#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>
#include "TMath.h"

using namespace dd4hep;
using namespace dd4hep::rec;

/** \addtogroup PID Particle ID Detectors 
 */
/** \addtogroup ThresholdGasCherenkov Light Gas (threshold) Cherenkov detector.
 * \brief Type: **ThresholdGasCherenkov**.
 * \ingroup PID
 *
 * \code
 *   <detector>
 *   </detector>
 * \endcode
 *
 * @{
 */
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement x_det    = handle;
    std::string     det_name = x_det.nameStr();
    int             det_id   = x_det.id();
    DetElement      det(det_name, det_id);
    // sensitive detector type
    sens.setType("tracker");

    auto dims = x_det.dimensions();
    int  nsec = dims.numsides();
    xml_dim_t dims_tank = dims.child(_Unicode(main));
    xml_dim_t dims_snout = dims.child(_Unicode(snout));

    xml_dim_t x_place = x_det.child(_U(placement));
    auto pos_x = x_place.x();
    auto pos_y = x_place.y();
    auto pos_z = x_place.z();

    // Everything that goes in the tank will be copies of the sector assembly volume
    Assembly v_sector("cherenkov_sector_1");
    DetElement de_sector("de_sector" + std::to_string(1), 1);

    // gas tank
    auto        x_rad   = x_det.child(_U(radiator));
    auto        rad_mat = desc.material(x_rad.attr<std::string>(_U(material)));
    ConeSegment tank_main(dims_tank.length()/2., dims_tank.rmin1(), dims_tank.rmax1(), dims_tank.rmin2(), dims_tank.rmax2());
    ConeSegment tank_snout(dims_snout.length()/2., dims_snout.rmin1(), dims_snout.rmax1(), dims_snout.rmin2(), dims_snout.rmax2());
    UnionSolid  tank_solid(tank_main,tank_snout,Position(0, 0, -(dims_tank.length() + dims_snout.length())/2.));
    Volume      v_tank("vol_gas_tank", tank_solid, rad_mat);
    v_tank.setVisAttributes(desc, dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BlueVis"));
    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(v_tank, Position(pos_x, pos_y, pos_z));
    envPV.addPhysVolID("system", det_id);
    det.setPlacement(envPV);

    // optical surface manager
    OpticalSurfaceManager surfMgr = desc.surfaceManager();

    // mirrors
    auto x_mirrors = x_det.child(_Unicode(mirrors));
    int i = 1;
    for (xml_coll_t il(x_mirrors, Unicode("piece")); il; ++il) {
        xml_comp_t x_mir = il;
        xml_dim_t mdim = x_mir.child(_U(dimensions));
        xml_dim_t mloc = x_mir.child(_U(placement));
        xml_dim_t mrot = x_mir.child(_U(rotation));
        auto      mmat = desc.material(x_mir.materialStr());

        Sphere mir_shell(mdim.radius(), mdim.radius() + mdim.thickness(), 0., M_PI/2.);
        Trd1   mir_cutout(mdim.attr<double>(_Unicode(width1))/2., mdim.attr<double>(_Unicode(width2))/2.,
                          mdim.length()/2., mdim.length()/2.);
        auto   mir_trans = RotationX(M_PI/2.)*Transform3D(Position(0., 0., -mdim.radius()));

        Volume v_mir("vol_mirror_" + std::to_string(i), IntersectionSolid(mir_cutout, mir_shell, mir_trans), mmat);
        auto   mir_trans2 = Transform3D(Position(mloc.x(), mloc.y(), mloc.z()))*RotationZYX(mrot.z(), mrot.y(), mrot.x());
        PlacedVolume pv_mir = v_sector.placeVolume(v_mir, mir_trans2);
        DetElement de_mir(det, "de_mirror" + std::to_string(i), i);
        pv_mir.addPhysVolID("mirror", i);
        de_mir.setPlacement(pv_mir);
        // v_mir.setSensitiveDetector(sens);

        auto surface = surfMgr.opticalSurface(x_mir.attr<std::string>(_Unicode(surface)));

        // optical surface
        SkinSurface mirror_skin(desc, de_mir, "mirror_surface_" + std::to_string(i), surface, v_mir);
        mirror_skin.isValid();
        i++;
    }

    // sectors
    double sector_angle = 2.*M_PI / nsec;
    for (int isec = 1; isec <= nsec; isec++) {
        auto pv = v_tank.placeVolume(v_sector, Transform3D(RotationZ((isec - 1) * sector_angle)));
        pv.addPhysVolID("sector", isec);
        auto amod = (isec == 1 ? de_sector : de_sector.clone("de_sector" + std::to_string(isec), isec));
        amod.setPlacement(pv);
        det.add(amod);
    }

    // ---------------
    // Winston Cone 
    auto x_winston   = x_det.child(_Unicode(winston_cone));
    auto winston_mat = desc.material(x_winston.attr<std::string>(_Unicode(material)));

    xml_dim_t wpl  = x_winston.child(_U(placement));
    xml_dim_t wrot = x_winston.child(_U(rotation));

    xml_dim_t cdims             = x_winston.child(_Unicode(cone_dimensions));
    double    cone_thickness    = cdims.thickness();
    double    cone_length       = cdims.attr<double>(_Unicode(length));
    double    cone_radius1      = cdims.attr<double>(_Unicode(radius1));
    double    cone_radius2      = cdims.attr<double>(_Unicode(radius2));
    // double    cone_inset_length = cdims.attr<double>(_Unicode(inset_length));

    xml_dim_t tdims       = x_winston.child(_Unicode(tube_dimensions));
    double    tube_radius = tdims.radius();
    double    tube_length = tdims.length();

    DetElement       de_winston_cone(det, "de_winston_cone1", 1);
    Tube             winston_tube(tube_radius, tube_radius + cone_thickness, tube_length / 2.0);
    Paraboloid       winston_cone1(cone_radius1 + cone_thickness, cone_radius2 + cone_thickness, cone_length / 2.0 );
    Paraboloid       winston_cone2(cone_radius1, cone_radius2, cone_length / 2.0 );
    SubtractionSolid winston_cone(winston_cone1, winston_cone2);

    Volume v_winston_cone_solid("v_winston_cone_solid", winston_cone, winston_mat);
    PlacedVolume pv_winston_cone_solid = v_sector.placeVolume(
        v_winston_cone_solid, Transform3D(Position(wpl.x(), wpl.y(), wpl.z())) *
                       RotationZYX(wrot.z(), wrot.y(), wrot.x()) *
                       Transform3D(Position(0, 0, tube_length / 2.0 + 5.0 * mm)));
    de_winston_cone.setPlacement(pv_winston_cone_solid);
    // optical surface
    auto surface = surfMgr.opticalSurface(x_winston.attr<std::string>(_Unicode(surface)));
    SkinSurface winston_skin(desc, de_winston_cone, "winston_surface", surface, v_winston_cone_solid);
    winston_skin.isValid();

    // ---------------
    // Dummy PMT surface
    auto x_pmt = x_det.child(_Unicode(pmt_array));
    xml_dim_t dims_pmt = x_pmt.child(_Unicode(dimensions));
    auto      pmt_x    = dims_pmt.x();
    auto      pmt_y    = dims_pmt.y();
    auto      pmt_surf = surfMgr.opticalSurface(x_pmt.attr<std::string>(_Unicode(surface)));

    DetElement   de_pmt_array(det, "PMT_DE", 1);
    Box          pmt_array(pmt_x/2., pmt_y/2., 5 * mm / 2.0);
    Volume       v_pmt_array("v_pmt_array", pmt_array, rad_mat);
    PlacedVolume pv_pmt_array =
        v_sector.placeVolume(v_pmt_array, Transform3D(Position(wpl.x(), wpl.y(), wpl.z())) *
                                          RotationZYX(wrot.z(), wrot.y(), wrot.x()));

    pv_pmt_array.addPhysVolID("mirror", 3);
    de_pmt_array.setPlacement(pv_pmt_array);
    v_pmt_array.setSensitiveDetector(sens);

    // optical surface
    SkinSurface pmt_skin(desc, de_pmt_array, "LGCPMTsurface", pmt_surf, v_pmt_array);
    pmt_skin.isValid();

    // copper layer inside to stop photons
    Box  pmt_array_backing(pmt_x/2., pmt_y/2., 1*mm/2.0);
    auto Copper = desc.material("Copper");
    Volume v_pmt_array_backing("v_pmt_array_backing", pmt_array_backing, Copper);
    // PlacedVolume pv_pmt_array_backing = 
    v_pmt_array.placeVolume(v_pmt_array_backing, Position(0,0,0));

    return det;
}
//@}
// clang-format off
DECLARE_DETELEMENT(SoLID_GasCherenkov, createDetector)

