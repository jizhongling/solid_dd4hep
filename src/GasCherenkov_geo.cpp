#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <XML/Helper.h>
#include "TMath.h"
#include "DDRec/Surface.h"

#include "DD4hep/OpticalSurfaces.h"
#include "DDRec/DetectorData.h"

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
    xml_det_t x_det = handle;

    auto det_name = x_det.nameStr();
    auto dims = x_det.dimensions();
    auto r0 = dims.rmin();
    auto r1 = dims.rmax1();
    auto r2 = dims.rmax2();
    auto zmin = dims.zmin();
    auto zmax = dims.zmax();

    xml_dim_t x_place = x_det.child(_U(placement));
    auto pos_x = x_place.x();
    auto pox_y = x_place.y();
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
    double LGC_mirror1_radius = 277.51*cm;
    double LGC_mirror2_radius = 157.99*cm;
    double LGC_mirror1_length = 114.53*cm;
    double LGC_mirror2_length = 59.260*cm;
    double LGC_mirror1_width1 = 16.26*cm;
    double LGC_mirror1_width2 = 36.03*cm;
    double LGC_mirror2_width1 = 37.06*cm;
    double LGC_mirror2_width2 = 45.95*cm;
    double LGC_mirror1_thickness = 2.0*mm;
    double LGC_mirror2_thickness = 2.0*mm;

    double LGC_sector_angle       = M_PI * 15.0 / 180.0;
    //double LGC_scattering_angle   = 11.0 * M_PI / 180.0;

    /*
    double LGC_mirror1_tilt_angle = mirror_rotations[1][0];//25.0 * M_PI / 180.0;
    double LGC_mirror2_tilt_angle = mirror_rotations[2][0];//2.0 * M_PI / 180.0;
    double LGC_pmt_tilt_angle     = mirror_rotations[3][0];//45.0 * M_PI / 180.0;

    double LGC_pmt_z_pos = mirror_positions[3].z();//-30.0*cm;
    double LGC_pmt_y_pos = mirror_positions[3].y();//LGC_outer_radius1 - 20.0*cm;
    */
    double LGC_pmt_array_size = 20.0*cm;

    // Everything that goes in the tank will be copies of the sector assembly volume
    Assembly v_sector("cherenkov_sector_1");
    DetElement de_sector("de_sector"+std::to_string(1), 1);

    // build gas tank
    auto        x_rad   = x_det.child(_U(radiator));
    auto        rad_mat = desc.material(dd4hep::getAttrOrDefault<std::string>(x_rad, _U(material), "N2Optical"));
    ConeSegment tank_main(0.5 * LGC_main_length, LGC_inner_radius1, LGC_outer_radius1,
                          LGC_inner_radius2, LGC_outer_radius1);
    ConeSegment tank_snout(0.5 * LGC_snout_length, LGC_snout_inner_radius1, LGC_snout_outer_radius1,
                           LGC_snout_inner_radius2, LGC_snout_outer_radius2);
    UnionSolid  sidis_tank(tank_main,tank_snout,Position(0, 0, -0.5 * LGC_main_length - 0.5 * LGC_snout_length));
    Volume      v_lgc_tank("v_lgc_tank_gas", sidis_tank, rad_mat);
    v_lgc_tank.setVisAttributes(desc, dd4hep::getAttrOrDefault(x_det, _Unicode(vis), "BlueVis"));

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
        auto   mir_trans2 = Transform3D(Position(0., mloc.y(), mloc.z()))*RotationX(mrot.x())*RotationY(mrot.y())*RotationZ(mrot.z());
        PlacedVolume pv_mir = v_sector.placeVolume(v_mir, mir_trans2);
    }

}
/*

    Material air = desc.air();
    Material PyrexGlass = desc.material("PyrexGlass");
    Material Copper = desc.material("Copper");

    // the gas tank

    // Everything that goes in the tank will be copies of the sector assembly volume
    Assembly v_sector("cherenkov_sector_1");
    DetElement de_sector("de_sector"+std::to_string(1),1);

    // mirrors
    Sphere mirror1_shell(LGC_mirror1_radius, LGC_mirror1_radius + LGC_mirror1_thickness, 
                         0.0, M_PI / 2);
    Trd1   mirror1_cutout(LGC_mirror1_width1 / 2.0, LGC_mirror1_width2 / 2.0,
                        LGC_mirror1_length / 2.0, LGC_mirror1_length / 2.0);
    IntersectionSolid mirror1_shape(mirror1_cutout, mirror1_shell,
                                    RotationX(M_PI/2.0)*Transform3D(Position(0, 0, -LGC_mirror1_radius)));
    Sphere mirror2_shell(LGC_mirror2_radius, LGC_mirror2_radius + LGC_mirror2_thickness, 
                         0.0, M_PI / 2);
    Trd1   mirror2_cutout(LGC_mirror2_width1 / 2.0, LGC_mirror2_width2 / 2.0,
                        LGC_mirror2_length / 2.0, LGC_mirror2_length / 2.0);
    IntersectionSolid mirror2_shape(mirror2_cutout, mirror2_shell,
                                    RotationX(M_PI/2.0)*Transform3D(Position(0, 0, -LGC_mirror2_radius)));

    double z_mirror1 = mirror_positions[1].z();
    double z_mirror2 = mirror_positions[2].z();
    double y_mirror1 = mirror_positions[1].y();
    double y_mirror2 = mirror_positions[2].y();

    Volume     v_mirror1_shape("v_mirror1_shape", mirror1_shape, PyrexGlass);
    PlacedVolume pv_mirror1_shape = v_sector.placeVolume(
        v_mirror1_shape, Transform3D(Position(0, y_mirror1, z_mirror1)) *
                             RotationX(-M_PI / 2.0 + LGC_mirror1_tilt_angle));

    DetElement   de_mirror1_shape(det,"de_mirror1_shape"+std::to_string(1),1);
    pv_mirror1_shape.addPhysVolID("mirror", 1);
    de_mirror1_shape.setPlacement(pv_mirror1_shape);
    sens.setType("photoncounter");
    v_mirror1_shape.setSensitiveDetector(sens);

    Volume     v_mirror2_shape("v_mirror2_shape", mirror2_shape, PyrexGlass);
    PlacedVolume pv_mirror2_shape = v_sector.placeVolume(
        v_mirror2_shape, Transform3D(Position(0, y_mirror2, z_mirror2)) *
                             RotationX(-M_PI / 2.0 + LGC_mirror2_tilt_angle));

    DetElement   de_mirror2_shape(det,"de_mirror2_shape"+std::to_string(2),2);
    pv_mirror2_shape.addPhysVolID("mirror", 2);
    de_mirror2_shape.setPlacement(pv_mirror2_shape);
    sens.setType("photoncounter");
    v_mirror2_shape.setSensitiveDetector(sens);

    // ---------------
    // Winston Cone 
    double LGC_winston_cone_thickness = 4*mm;
    double LGC_winston_tube_inner_radius = 11.28*cm;
    double LGC_winston_tube_length = 30.0*cm;
    double LGC_winston_cone_length = 30.0*cm;
    double LGC_winston_cone_inner_radius1 = 7.8*cm;
    double LGC_winston_cone_inner_radius2 = 21.0*cm;
    double LGC_winston_cone_inset_length = 7.90909*cm;
    DetElement   de_winston_cone(det,"de_winston_cone1",1);
    Tube         winston_tube(LGC_winston_tube_inner_radius,
                      LGC_winston_tube_inner_radius + LGC_winston_cone_thickness,
                      LGC_winston_tube_length / 2.0);
    //Cone        winston_cone(LGC_winston_cone_length / 2.0, LGC_winston_cone_inner_radius1,
    //                  LGC_winston_cone_inner_radius1 + LGC_winston_cone_thickness,
    //                  LGC_winston_cone_inner_radius2,
    //                  LGC_winston_cone_inner_radius2 + LGC_winston_cone_thickness );
    //UnionSolid  winston_cone_solid(winston_tube,winston_cone,Position(0,0,LGC_winston_tube_length / 2.0 - LGC_winston_cone_inset_length));
    Paraboloid winston_cone1(LGC_winston_cone_inner_radius1 + LGC_winston_cone_thickness,
                       LGC_winston_cone_inner_radius2 + LGC_winston_cone_thickness,
                       LGC_winston_cone_length / 2.0 );
    Paraboloid winston_cone2(LGC_winston_cone_inner_radius1,
                       LGC_winston_cone_inner_radius2,
                       LGC_winston_cone_length / 2.0 );

    SubtractionSolid  winston_cone(winston_cone1, winston_cone2);

    Volume v_winston_cone_solid("v_winston_cone_solid", winston_cone, PyrexGlass);
    PlacedVolume pv_winston_cone_solid = v_sector.placeVolume(
        v_winston_cone_solid, Transform3D(Position(0, LGC_pmt_y_pos, LGC_pmt_z_pos)) *
                       RotationX(LGC_pmt_tilt_angle) *
                       Transform3D(Position(0, 0, LGC_winston_tube_length / 2.0 + 5.0 * mm)));

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


    // all sectors
    for (int i_sector = 1; i_sector <= 30; i_sector++) {
      //std::cout << i_sector  << " sector\n";
      PlacedVolume pv =
          v_lgc_tank.placeVolume(v_sector, Transform3D(RotationZ((i_sector - 1) * LGC_sector_angle)));
      pv.addPhysVolID("sector", i_sector);
      auto amod = (i_sector == 1 ? de_sector : de_sector.clone("de_sector" + std::to_string(i_sector), i_sector));
      amod.setPlacement(pv);
      det.add(amod);
    }


    //// ---------------


    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(v_lgc_tank, Position(pos_x, pos_y, pos_z));
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);

    return det;
}
*/
//@}
// clang-format off
DECLARE_DETELEMENT(SoLID_GasCherenkov, createDetector)

