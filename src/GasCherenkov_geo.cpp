// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Chao Peng

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>
#include <fmt/core.h>

using namespace dd4hep;
using namespace dd4hep::rec;

// A helper function to build stacked cone segments
// set solid and material into <vol>
// return the distance between (in z) the starting point and the volume center
double stackConeSegments(Volume &vol, xml_comp_t x_comp, const Material &mat)
{
    // store solids and their lengths
    std::vector<double> lengths;
    std::vector<ConeSegment> conesegs;

    // collect cone segments
    for (xml_coll_t il(x_comp, _Unicode(segment)); il; ++il) {
        xml_dim_t idim = il;
        conesegs.emplace_back(idim.length()/2., idim.rmin1(), idim.rmax1(), idim.rmin2(), idim.rmax2());
        lengths.push_back(idim.length());
    }
    // make a union solid out of the segments
    // 0 size also include and will raise errors
    if (conesegs.size() <= 1) {
        vol.setSolid(conesegs[0]);
    } else {
        UnionSolid segs_union(conesegs[0], conesegs[1], Position(0., 0., (lengths[0] + lengths[1])/2.));
        double mid_length = 0.;
        for (size_t i = 2; i < conesegs.size(); ++i) {
            mid_length += lengths[i - 1];
            segs_union = UnionSolid(segs_union, conesegs[i], Position(0., 0., (lengths[0] + lengths[i])/2. + mid_length));
        }
        vol.setSolid(segs_union);
    }
    vol.setMaterial(mat);
    // we use the first segment as the center solid
    return lengths[0]/2.;
}

// main geometry builder
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement x_det    = handle;
    std::string     det_name = x_det.nameStr();
    int             det_id   = x_det.id();
    DetElement      det(det_name, det_id);
    // sensitive detector type
    sens.setType("tracker");

    int  nsec    = x_det.numsides();

    xml_dim_t x_place = x_det.child(_U(placement));
    auto pos_x        = x_place.x();
    auto pos_y        = x_place.y();
    auto pos_z0       = x_place.z0();

    // --------------
    // Main tank
    // build from stacking conesegments
    auto x_tank = x_det.child(_Unicode(tank)); 
    Volume v_tank("v_tank");
    // using tank center (z) as the center point for all following volumes
    double shift_z = stackConeSegments(v_tank, x_tank, desc.material(x_tank.attr<std::string>(_Unicode(material))));
    v_tank.setVisAttributes(desc, x_tank.attr<std::string>(_Unicode(vis)));

    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(v_tank, Position(pos_x, pos_y, pos_z0 + shift_z));
    envPV.addPhysVolID("system", det_id);
    det.setPlacement(envPV);

    // build the radiator inside tank
    auto x_rad = x_tank.child(_Unicode(radiator));
    Volume v_rad("v_gas_radiator");
    // build radiator volume
    stackConeSegments(v_rad, x_rad, desc.material(x_rad.attr<std::string>(_Unicode(material))));
    v_rad.setVisAttributes(desc, x_rad.attr<std::string>(_Unicode(vis)));
    xml_dim_t rpl = x_rad.child(_U(placement));
    v_tank.placeVolume(v_rad, Position(rpl.x(), rpl.y(), rpl.z()));
    // shift from radiator center to the detector z0
    shift_z += rpl.z();

    // ---------------
    // Sectors
    // Everything that goes in the tank will be copies of the sector assembly volume
    Assembly v_sector("cherenkov_sector");
    DetElement de_sector("de_sector", 1);
    double sector_angle = 2.*M_PI / nsec;
    for (int isec = 1; isec <= nsec; isec++) {
        auto pv = v_rad.placeVolume(v_sector, Transform3D(RotationZ((isec - 1) * sector_angle)));
        pv.addPhysVolID("sector", isec);
        auto amod = (isec == 1 ? de_sector : de_sector.clone("de_sector" + std::to_string(isec), isec));
        amod.setPlacement(pv);
        det.add(amod);
    }
    // optical surface manager
    OpticalSurfaceManager surfMgr = desc.surfaceManager();

    // ---------------
    // Mirrors: intersection between a spherical shell and a trapezoid
    auto x_mirs = x_det.child(_Unicode(mirrors));
    for (xml_coll_t il(x_mirs, _Unicode(piece)); il; ++il) {
        xml_comp_t x_mir = il;
        auto       mid   = x_mir.id();
        xml_dim_t  mpl   = x_mir.child(_U(placement));
        auto       mmat  = desc.material(x_mir.materialStr());

        // build spherical shell
        xml_dim_t x_msph = x_mir.child(_Unicode(shell));
        // front half, this is to help debugging visualization
        Sphere mshell(x_msph.rmin(), x_msph.rmax(),
            x_msph.attr<double>(_Unicode(theta0)), x_msph.attr<double>(_Unicode(dtheta)),
            x_msph.attr<double>(_Unicode(phi0)), x_msph.attr<double>(_Unicode(dphi)));
        // build trapezoid
        xml_dim_t x_mtrd = x_mir.child(_Unicode(wedge));
        Trd1 mwedge(x_mtrd.attr<double>(_Unicode(dx1)), x_mtrd.attr<double>(_Unicode(dx2)), x_mtrd.dy(), x_mtrd.dz());
        auto mwedge_trans = Transform3D(Position(x_mtrd.x(), x_mtrd.y(), x_mtrd.z()))\
                          * RotationZYX(x_mtrd.attr<double>(_Unicode(rotz)),
                                        x_mtrd.attr<double>(_Unicode(roty)),
                                        x_mtrd.attr<double>(_Unicode(rotz)));

        // mirror volume
        // use union to debug if intersection does not exist
        // UnionSolid mir_solid(mshell, mwedge, mwedge_trans);
        IntersectionSolid mir_solid(mshell, mwedge, mwedge_trans);
        Volume            v_mir("v_mirror_" + std::to_string(mid), mir_solid, mmat);
        PlacedVolume      pv_mir = v_sector.placeVolume(v_mir, Position(mpl.x(), mpl.y(), mpl.z() - shift_z));
        DetElement        de_mir(det, "de_mirror_" + std::to_string(mid), mid);
        de_mir.setPlacement(pv_mir);
        v_mir.setVisAttributes(desc, x_mir.attr<std::string>(_Unicode(vis)));

        // optical surface
        auto msurf = surfMgr.opticalSurface(x_mir.attr<std::string>(_Unicode(surface)));
        SkinSurface mirror_skin(desc, de_mir, "mirror_surface_" + std::to_string(mid), msurf, v_mir);
        mirror_skin.isValid();
    }

    // ---------------
    // Winston Cone 
    auto x_winston   = x_det.child(_Unicode(winston_cone));
    auto winston_mat = desc.material(x_winston.attr<std::string>(_Unicode(material)));

    // build an assembly as its envelope
    Assembly     winston_assem("winston_assembly");
    xml_dim_t    wpl        = x_winston.child(_U(placement));
    xml_dim_t    wrot       = x_winston.child(_U(rotation));
    auto         wtrans     = Transform3D(Position(wpl.x(), wpl.y(), wpl.z() - shift_z))\
                            * RotationZYX(wrot.z(), wrot.y(), wrot.x());
    PlacedVolume pv_winston = v_sector.placeVolume(winston_assem, wtrans);
    DetElement   de_winston(det, "de_winston_assembly", 1);
    de_winston.setPlacement(pv_winston);

    // build pmt array (dummy sensitive surface at the assembly's center)
    // TODO: implement realistic material layers
    xml_comp_t x_pmt  = x_winston.child(_Unicode(pmt_array));
    auto       pmt_dx = x_pmt.attr<double>(_Unicode(dx));
    auto       pmt_dy = x_pmt.attr<double>(_Unicode(dy));
    double     pmt_dz = 2.*mm;
    DetElement de_pmt_array(det, "PMT_DE", 1);
    Box        pmt_array(pmt_dx/2., pmt_dy/2., pmt_dz/2.);
    Volume     v_pmt_array("v_pmt_array", pmt_array, winston_mat);
    if (x_pmt.isSensitive()) {
        v_pmt_array.setSensitiveDetector(sens);
    }
    PlacedVolume pv_pmt_array = winston_assem.placeVolume(v_pmt_array, Position(0., 0., 0.));
    pv_pmt_array.addPhysVolID("module", 1);
    de_pmt_array.setPlacement(pv_pmt_array);
    v_pmt_array.setVisAttributes(desc, x_pmt.attr<std::string>(_Unicode(vis)));

    // build cone (its end touches the PMT surface)
    xml_dim_t x_cone            = x_winston.child(_Unicode(cone));
    double    cone_thickness    = x_cone.thickness();
    double    cone_length       = x_cone.length();
    double    cone_radius1      = x_cone.rmin();
    double    cone_radius2      = x_cone.rmax();
    auto      cone_shape        = x_cone.attr<std::string>(_Unicode(shape));

    Volume    v_winston_cone("v_winston_cone");
    if (cone_shape == "paraboloid") {
        Paraboloid       winston_cone1(cone_radius1 + cone_thickness, cone_radius2 + cone_thickness, cone_length / 2.0 );
        Paraboloid       winston_cone2(cone_radius1, cone_radius2, cone_length / 2.0 + 0.1*mm );
        SubtractionSolid winston_cone_solid(winston_cone1, winston_cone2);
        v_winston_cone.setSolid(winston_cone_solid);
    } else if (cone_shape == "cone") {
        Cone winston_cone_solid(cone_length/2., cone_radius1, cone_radius1 + cone_thickness, cone_radius2, cone_radius2 + cone_thickness);
        v_winston_cone.setSolid(winston_cone_solid);
    } else {
        printout(ERROR, "SoLID_GasCherenkov",
                 fmt::format("Unknown shape {} for the winston cone, please use (Paraboloid, Cone).", cone_shape));
        throw std::runtime_error("Failed to build winston cone solid for SoLID_GasCherenkov");
    }
    v_winston_cone.setMaterial(winston_mat);
    DetElement       de_winston_cone(det, "de_winston_cone", 1);
    PlacedVolume     pv_winston_cone = winston_assem.placeVolume(v_winston_cone, Position(0., 0., (cone_length + pmt_dz)/2.));
    de_winston_cone.setPlacement(pv_winston_cone);
    v_winston_cone.setVisAttributes(desc, x_cone.attr<std::string>(_Unicode(vis)));

    // optical surface
    auto wsurf = surfMgr.opticalSurface(x_winston.attr<std::string>(_Unicode(surface)));
    SkinSurface winston_skin(desc, de_winston_cone, "winston_surface", wsurf, v_winston_cone);
    winston_skin.isValid();

    // build shield
    // double    cone_inset_length = cdims.attr<double>(_Unicode(inset_length));
    xml_dim_t x_shield      = x_winston.child(_Unicode(shield));
    double    shield_rmin   = x_shield.radius();
    double    shield_rmax   = shield_rmin + x_shield.thickness();
    double    shield_length = x_shield.length();
    auto      shield_mat    = desc.material(x_shield.attr<std::string>(_Unicode(material)));
    Tube      winston_shield_solid(shield_rmin, shield_rmax, shield_length / 2.0);
    Volume    v_winston_shield("v_winston_shield", winston_shield_solid, shield_mat);
    // wrapping around PMT
    winston_assem.placeVolume(v_winston_shield, Position(0., 0., x_shield.attr<double>(_Unicode(shift_z))));
    v_winston_shield.setVisAttributes(desc, x_shield.attr<std::string>(_Unicode(vis)));

    return det;
}
DECLARE_DETELEMENT(SoLID_GasCherenkov, createDetector)

