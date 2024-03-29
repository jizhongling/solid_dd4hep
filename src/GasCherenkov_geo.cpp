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
    sens.setType("photoncounter");

    auto dims = x_det.dimensions();
    auto r0   = dims.rmin();
    auto r1   = dims.rmax1();
    auto r2   = dims.rmax2();
    auto zmin = dims.zmin();
    auto zmax = dims.zmax();
    int  nsec = dims.numsides();

    xml_dim_t x_place = x_det.child(_U(placement));
    auto pos_x = x_place.x();
    auto pos_y = x_place.y();
    auto pos_z = x_place.z();

    double LGC_inner_radius1 = 71.0*cm;
    double LGC_inner_radius2 = 85.0*cm;
    double LGC_outer_radius1 = 265.0*cm;
    double LGC_main_length   = 105.0*cm;
    double LGC_snout_length   = 107.0*cm;
    double LGC_snout_inner_radius1 = 58.0*cm; 
    double LGC_snout_inner_radius2 = LGC_inner_radius1; 
    double LGC_snout_outer_radius1 = 127.0*cm;
    double LGC_snout_outer_radius2 = 144.0*cm;
    double LGC_entrance_window_thickness = 0.05*mm; // something tells this might be 5 mil, not mm
    double LGC_exit_window_thickness = 0.1*mm; // same here
    double LGC_pmt_array_size = 20.0*cm;

    // Everything that goes in the tank will be copies of the sector assembly volume
    Assembly v_sector("cherenkov_sector_1");
    DetElement de_sector("de_sector" + std::to_string(1), 1);

    // gas tank
    auto        x_rad   = x_det.child(_U(radiator));
    auto        rad_mat = desc.material(x_rad.attr<std::string>(_U(material)));
    ConeSegment tank_main(0.5 * LGC_main_length, LGC_inner_radius1, LGC_outer_radius1,
                          LGC_inner_radius2, LGC_outer_radius1);
    ConeSegment tank_snout(0.5 * LGC_snout_length, LGC_snout_inner_radius1, LGC_snout_outer_radius1,
                           LGC_snout_inner_radius2, LGC_snout_outer_radius2);
    UnionSolid  tank_solid(tank_main,tank_snout,Position(0, 0, -0.5 * LGC_main_length - 0.5 * LGC_snout_length));
    Volume      v_tank("vol_gas_tank", tank_solid, rad_mat);
    v_tank.setVisAttributes(desc, dd4hep::getAttrOrDefault<std::string>(x_det, _Unicode(vis), "BlueVis"));
    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(v_tank, Position(pos_x, pos_y, pos_z));
    envPV.addPhysVolID("system", det_id);
    det.setPlacement(envPV);

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
        auto   mir_trans2 = Transform3D(Position(0., mloc.y(), mloc.z()))*RotationZYX(mrot.z(), mrot.y(), mrot.x());
        PlacedVolume pv_mir = v_sector.placeVolume(v_mir, mir_trans2);
    }

    // sectors
    double sector_angle = 2.*M_PI / nsec;
    for (int isec = 1; isec <= nsec; isec++) {
        auto pv = v_tank.placeVolume(v_sector, Transform3D(RotationZ((isec - 1) * sector_angle)));
        pv.addPhysVolID("sector" + std::to_string(isec), isec);
        auto amod = (isec == 1 ? de_sector : de_sector.clone("de_sector" + std::to_string(isec), isec));
        amod.setPlacement(pv);
        det.add(amod);
    }

    // ---------------
    // Winston Cone 
    auto   x_winston    = x_det.child(_Unicode(winston_cone));

    xml_dim_t cdims       = x_winston.child(_Unicode(cone_dimensions));
    double cone_thickness = cdims.thickness();
    double cone_length    = cdims.length1();
    double cone_radius1   = cdims.radius1();
    double cone_radius2   = cdims.radius2();
    double cone_inset_length = cdims.inset_length();

    xml_dim_t tdims    = x_winston.child(_Unicode(tube_dimensions));
    double tube_radius = tdims.radius();
    double tube_length = tdims.length();
    auto   winston_mat = desc.material(x_winstone.materialStr());

    DetElement   de_wcone(det, "de_winston_cone1", 1);
    Tube         winston_tube(LGC_winston_tube_inner_radius,
                      LGC_winston_tube_inner_radius + LGC_winston_cone_thickness,
                      LGC_winston_tube_length / 2.0);
    Paraboloid winston_cone1(LGC_winston_cone_inner_radius1 + LGC_winston_cone_thickness,
                       LGC_winston_cone_inner_radius2 + LGC_winston_cone_thickness,
                       LGC_winston_cone_length / 2.0 );
    Paraboloid winston_cone2(LGC_winston_cone_inner_radius1,
                       LGC_winston_cone_inner_radius2,
                       LGC_winston_cone_length / 2.0 );
    SubtractionSolid  winston_cone(winston_cone1, winston_cone2);

    Volume v_winston_cone_solid("v_winston_cone_solid", winston_cone, winston_mat);
    PlacedVolume pv_winston_cone_solid = v_sector.placeVolume(
        v_winston_cone_solid, Transform3D(Position(0, LGC_pmt_y_pos, LGC_pmt_z_pos)) *
                       RotationX(LGC_pmt_tilt_angle) *
                       Transform3D(Position(0, 0, LGC_winston_tube_length / 2.0 + 5.0 * mm)));

    //// ---------------
    return det;
}
/*
    //std::cout << " LGC_pmt_y_pos/cm " << LGC_pmt_y_pos/cm  << "\n";
    //std::cout << " LGC_pmt_z_pos/cm " << LGC_pmt_z_pos/cm  << "\n";
    //mirrorPV.addPhysVolID("layer", 2).addPhysVolID("module", 1);
    //mirror_DE.setPlacement(mirrorPV);
    //sens.setType("photoncounter");
    //mirrorVol.setSensitiveDetector(sens);
  
    //// ---------------

    //// ---------------
    //// Dummy PMT surface
    DetElement   de_pmt_array(det, "PMT_DE", 1);
    Box          pmt_array(LGC_pmt_array_size / 2.0, LGC_pmt_array_size / 2.0, 5 * mm / 2.0);
    Volume       v_pmt_array("v_pmt_array", pmt_array, rad_mat);
    PlacedVolume pv_pmt_array =
        v_sector.placeVolume(v_pmt_array, Transform3D(Position(0, LGC_pmt_y_pos, LGC_pmt_z_pos)) *
                                                RotationX(LGC_pmt_tilt_angle));

    pv_pmt_array.addPhysVolID("mirror", 3);
    de_pmt_array.setPlacement(pv_pmt_array);
    sens.setType("photoncounter");
    v_pmt_array.setSensitiveDetector(sens);

    // copper layer inside to stop photons
    Box  pmt_array_backing(LGC_pmt_array_size/2.0, LGC_pmt_array_size/2.0, 1*mm/2.0);
    Volume v_pmt_array_backing("v_pmt_array_backing", pmt_array_backing, Copper);
    PlacedVolume pv_pmt_array_backing = v_pmt_array.placeVolume(v_pmt_array_backing, Position(0,0,0));

    // Optical Surfaces

    OpticalSurfaceManager surfMgr = desc.surfaceManager();
    OpticalSurface mirrorSurf  = surfMgr.opticalSurface("MirrorOpticalSurface");
    OpticalSurface pmtSurf    = surfMgr.opticalSurface("PMTOpticalSurface");
    //BorderSurface  mirrorBorder_Surf   = BorderSurface(desc, det, "RICHmirror", mirrorSurf, mirrorPV,   envPV);
    SkinSurface mirrorBorder_Surf(desc,de_mirror1_shape,"LGCmirror", mirrorSurf, v_mirror1_shape);
    SkinSurface winstonBorder_Surf(desc,de_winston_cone,"LGCWinstonCone", mirrorSurf, v_winston_cone_solid);
    SkinSurface pmtBorder_Surf(desc,de_pmt_array,"LGCPMTsurface", pmtSurf, v_pmt_array);
    //BorderSurface  bubbleSurf = BorderSurface(description, sdet, "TankBubble", airSurf,   bubblePlace, tankPlace);
    mirrorBorder_Surf.isValid();
    winstonBorder_Surf.isValid();
    pmtBorder_Surf.isValid();
    //tankSurf.isValid();

}
*/
//@}
// clang-format off
DECLARE_DETELEMENT(SoLID_GasCherenkov, createDetector)

