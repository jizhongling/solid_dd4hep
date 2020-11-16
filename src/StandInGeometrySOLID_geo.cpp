#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"

using namespace std;
using namespace dd4hep;
//using namespace dd4hep::Geometry;


static Ref_t create_detector(Detector& lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = lcdd.air();
  DetElement sdet        (det_name,x_det.id());
  Assembly   assembly    ("SOLID_standin_assembly");

  Material steel = lcdd.material("Steel235");

  PlacedVolume pv;

  //assembly.setVisAttributes(lcdd, "OrangeVis");


  int n = 0;
  using namespace dd4hep;

  // ---------------------------------------------
  // Solenoid coil
  std::vector<double> rInner = {
    152.30*cm, 152.30*cm
  };
  std::vector<double> rOuter = { 
    154.30*cm, 154.30*cm 
  };
  std::vector<double> zPlane = { 
    -173.80*cm, 173.80*cm  
  };

  Polycone solenoid( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position solenoid_pos = {0*cm,0.0*cm,0.0*cm};
  Volume solenoid_vol("solenoid_vol", solenoid, air);
  solenoid_vol.setVisAttributes(lcdd, "GrayVis");
  pv = assembly.placeVolume(solenoid_vol, solenoid_pos);

  // ---------------------------------------------
  // Barrel yoke
  rInner = {
    177.91*cm,177.91*cm,176.60*cm,176.60*cm
  };
  rOuter = { 
    212.60*cm,212.60*cm,212.60*cm,212.60*cm
  };
  zPlane = { 
    -266.50*cm,-189.00*cm,-189.00*cm,189.00*cm
  };

  Polycone solenoid_barrel_yoke( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position solenoid_barrel_yoke_pos = {0*cm,0.0*cm,0.0*cm};
  Volume solenoid_barrel_yoke_vol("solenoid_barrel_yoke_vol", solenoid_barrel_yoke, steel);
  solenoid_barrel_yoke_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(solenoid_barrel_yoke_vol, solenoid_barrel_yoke_pos);

  // ---------------------------------------------
  // Barrel yoke_outer
  //
  rInner = {
    221.51*cm,221.51*cm
  };
  rOuter = { 
    257.50*cm,257.50*cm
  };
  zPlane = { 
    -266.50*cm,189.00*cm
  };

  Polycone solenoid_barrel_yoke_outer( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position solenoid_barrel_yoke_outer_pos = {0*cm,0.0*cm,0.0*cm};
  Volume solenoid_barrel_yoke_outer_vol("solenoid_barrel_yoke_outer_vol", solenoid_barrel_yoke_outer, steel);
  solenoid_barrel_yoke_outer_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(solenoid_barrel_yoke_outer_vol, solenoid_barrel_yoke_outer_pos);

  // ---------------------------------------------
  // solenoid_slab_spacer_upstream
  //
  rInner = {
    212.61*cm,212.61*cm
  };
  rOuter = { 
    221.50*cm,221.50*cm
  };
  zPlane = { 
    -266.50*cm,-235.90*cm
  };

  Polycone solenoid_slab_spacer_upstream( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position solenoid_slab_spacer_upstream_pos = {0*cm,0.0*cm,0.0*cm};
  Volume solenoid_slab_spacer_upstream_vol("solenoid_slab_spacer_upstream_vol", solenoid_slab_spacer_upstream, steel);
  solenoid_slab_spacer_upstream_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(solenoid_slab_spacer_upstream_vol, solenoid_slab_spacer_upstream_pos);

  // ---------------------------------------------
  // coil_collar_downstream
  //
  rInner = {
    144.00*cm,144.00*cm,156.00*cm
  };
  rOuter = { 
    285.00*cm,285.00*cm,285.00*cm
  };
  zPlane = { 
    189.01*cm,193.00*cm,209.00*cm
  };
  Polycone coil_collar_downstream( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position coil_collar_downstream_pos = {0*cm,0.0*cm,0.0*cm};
  Volume coil_collar_downstream_vol("coil_collar_downstream_vol", coil_collar_downstream, steel);
  coil_collar_downstream_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(coil_collar_downstream_vol, coil_collar_downstream_pos);


  // ---------------------------------------------
  // endcap_donut
  //
  rInner = {
    270.00*cm,270.00*cm
  };
  rOuter = { 
    285.00*cm,285.00*cm
  };

  // From Jay: 10ft from magnet to back of EC
  double z_endcap_donut0 = 209.01*cm + 10*12*2.54*cm;
  zPlane = { 
    209.01*cm, z_endcap_donut0 // 485.00*cm
  };
  Polycone endcap_donut( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position endcap_donut_pos = {0*cm,0.0*cm,0.0*cm};
  Volume endcap_donut_vol("endcap_donut_vol", endcap_donut, steel);
  endcap_donut_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(endcap_donut_vol, endcap_donut_pos);

  // ---------------------------------------------
  // endcap_bottom_inner
  //
  rInner = {
    30.00*cm,30.00*cm
  };
  rOuter = { 
    285.00*cm,285.00*cm
  };
  zPlane = { 
    z_endcap_donut0,z_endcap_donut0 + 15.0*cm
  };
  Polycone endcap_bottom_inner( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to inner surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position endcap_bottom_inner_pos = {0*cm,0.0*cm,0.0*cm};
  Volume endcap_bottom_inner_vol("endcap_bottom_inner_vol", endcap_bottom_inner, steel);
  endcap_bottom_inner_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(endcap_bottom_inner_vol, endcap_bottom_inner_pos);

  // ---------------------------------------------
  // endcap_bottom_outer
  //
  rInner = {
30.00*cm,45.00*cm
  };
  rOuter = { 
  185.00*cm,170.00*cm
  };
  zPlane = { 
  z_endcap_donut0 + 15.01*cm ,
  z_endcap_donut0 + 30.0*cm 
  //500.01*cm,515.00*cm
  };
  Polycone endcap_bottom_outer( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position endcap_bottom_outer_pos = {0*cm,0.0*cm,0.0*cm};
  Volume endcap_bottom_outer_vol("endcap_bottom_outer_vol", endcap_bottom_outer, steel);
  endcap_bottom_outer_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(endcap_bottom_outer_vol, endcap_bottom_outer_pos);

  // ---------------------------------------------
  // endcap_nose
  //
  rInner = {
    20.00*cm,30.00*cm,30.00*cm
  };
  rOuter = { 
    60.00*cm,90.00*cm,90.00*cm
  };
  zPlane = { 
    z_endcap_donut0 - 296.0*cm,z_endcap_donut0 - 80.0*cm ,z_endcap_donut0
  };
  Polycone endcap_nose( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position endcap_nose_pos = {0*cm,0.0*cm,0.0*cm};
  Volume endcap_nose_vol("endcap_nose_vol", endcap_nose, steel);
  endcap_nose_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(endcap_nose_vol, endcap_nose_pos);

  // ---------------------------------------------
  // front_piece
  //
  rInner = {
    55.60*cm,70.00*cm
  };
  rOuter = { 
    144.00*cm,144.00*cm
  };
  zPlane = { 
    -237.00*cm,-207.00*cm
  };
  Polycone front_piece( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position front_piece_pos = {0*cm,0.0*cm,0.0*cm};
  Volume front_piece_vol("front_piece_vol", front_piece, steel);
  front_piece_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(front_piece_vol, front_piece_pos);

  // ---------------------------------------------
  // upstream_shield
  //
  rInner = {
    48.50*cm,50.50*cm
  };
  rOuter = { 
    144.00*cm,144.00*cm
  };
  zPlane = { 
    -250.50*cm,-246.50*cm
  };
  Polycone upstream_shield( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position upstream_shield_pos = {0*cm,0.0*cm,0.0*cm};
  Volume upstream_shield_vol("upstream_shield_vol", upstream_shield, steel);
  upstream_shield_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(upstream_shield_vol, upstream_shield_pos);

  // ---------------------------------------------
  // upstream_shield2
  //
  rInner = {
    44.60*cm,46.50*cm
  };
  rOuter = { 
    144.00*cm,144.00*cm
  };
  zPlane = { 
    -258.50*cm,-254.50*cm
  };
  Polycone upstream_shield2( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position upstream_shield2_pos = {0*cm,0.0*cm,0.0*cm};
  Volume upstream_shield2_vol("upstream_shield2_vol", upstream_shield2, steel);
  upstream_shield2_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(upstream_shield2_vol, upstream_shield2_pos);

  // ---------------------------------------------
  // upstream_shield3
  //
  rInner = {
    40.10*cm,42.70*cm
  };
  rOuter = { 
    144.00*cm,144.00*cm
  };
  zPlane = { 
    -266.50*cm,-262.50*cm
  };
  Polycone upstream_shield3( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position upstream_shield3_pos = {0*cm,0.0*cm,0.0*cm};
  Volume upstream_shield3_vol("upstream_shield3_vol", upstream_shield3, steel);
  upstream_shield3_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(upstream_shield3_vol, upstream_shield3_pos);

  // ---------------------------------------------
  // upstream_shield4
  //
  rInner = {
    36.80*cm,38.80*cm
  };
  rOuter = { 
    257.50*cm,257.50*cm
  };
  zPlane = { 
    -274.50*cm,-270.50*cm
  };
  Polycone upstream_shield4( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position upstream_shield4_pos = {0*cm,0.0*cm,0.0*cm};
  Volume upstream_shield4_vol("upstream_shield4_vol", upstream_shield4, steel);
  upstream_shield4_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(upstream_shield4_vol, upstream_shield4_pos);

  // ---------------------------------------------
  // upstream_shield5
  //
  rInner = {
    33.90*cm,35.30*cm
  };
  rOuter = { 
    257.50*cm,257.50*cm
  };
  zPlane = { 
    -280.50*cm,-277.50*cm
  };
  Polycone upstream_shield5( 0,  ///< Initial Phi starting angle
      360*deg,  ///< Total Phi angle
      rInner,        ///< Tangent distance to outer surface
      rOuter,        ///< Tangent distance to outer surface
      zPlane);       ///< z coordinate of corners
  Position upstream_shield5_pos = {0*cm,0.0*cm,0.0*cm};
  Volume upstream_shield5_vol("upstream_shield5_vol", upstream_shield5, steel);
  upstream_shield5_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(upstream_shield5_vol, upstream_shield5_pos);

  // ---------------------------------------------
  // cryostat_inner
  //
  Tube cryostat_inner( 144*cm,144.1*cm,187*cm);
  Position cryostat_inner_pos = {0*cm,0.0*cm,0.0*cm};
  Volume cryostat_inner_vol("cryostat_inner_vol", cryostat_inner, air);
  cryostat_inner_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(cryostat_inner_vol, cryostat_inner_pos);

  // ---------------------------------------------
  // cryostat_outer
  //
  Tube cryostat_outer(176.48*cm, 176.59*cm, 187*cm);
  Position cryostat_outer_pos = {0*cm,0.0*cm,0.0*cm};
  Volume cryostat_outer_vol("cryostat_outer_vol", cryostat_outer, air);
  cryostat_outer_vol.setVisAttributes(lcdd, "SteelVis");
  pv = assembly.placeVolume(cryostat_outer_vol, cryostat_outer_pos);

  // ---------------------------------------------
  // flange_upstream
  //
  Tube flange_upstream(144*cm,176.59*cm,0.995*cm);
  Position flange_upstream_pos = {0*cm,0.0*cm,-188.0*cm};
  Volume flange_upstream_vol("flange_upstream_vol", flange_upstream, air);
  flange_upstream_vol.setVisAttributes(lcdd, "BlueVis");
  pv = assembly.placeVolume(flange_upstream_vol, flange_upstream_pos);


  // ---------------------------------------------
  // flange_downstream
  //
  Tube flange_downstream(144*cm,176.59*cm,0.995*cm);
  Position flange_downstream_pos = {0*cm,0.0*cm,188.0*cm};
  Volume flange_downstream_vol("flange_downstream_vol", flange_downstream, air);
  flange_downstream_vol.setVisAttributes(lcdd, "BlueVis");
  pv = assembly.placeVolume(flange_downstream_vol, flange_downstream_pos);

  //______________________________________________________________________________


  //// ---------------------------------------------
  //// Light Gas Chrenkov
  //// ---------------------------------------------
  //// cher_lg_Tank
  ////
  //rInner = {
  //  65*cm,67*cm,67*cm,85*cm
  //};
  //rOuter = { 
  //  144*cm,155*cm,265*cm,265*cm
  //};
  //zPlane = { 
  //  194*cm,209.01*cm,209.01*cm,301*cm
  //};
  //Polycone cher_lg_tank( 0,  ///< Initial Phi starting angle
  //    360*deg,  ///< Total Phi angle
  //    rInner,        ///< Tangent distance to outer surface
  //    rOuter,        ///< Tangent distance to outer surface
  //    zPlane);       ///< z coordinate of corners
  //Position cher_lg_tank_pos = {0*cm,0.0*cm,0.0*cm};
  //Volume cher_lg_tank_vol("cher_lg_tank_vol", cher_lg_tank, air);
  //cher_lg_tank_vol.setVisAttributes(lcdd, "GreenVis");
  //pv = assembly.placeVolume(cher_lg_tank_vol, cher_lg_tank_pos);

  ////  |  root |   tank |   0*cm 0*cm 0*cm |   0*deg 0*deg 0*deg | FF9900 | Polycone |0*deg 360*deg 4 65*cm 67*cm 67*cm 85*cm 144*cm 155*cm 265*cm 265*cm 194*cm 209.01*cm 209.01*cm 301*cm | SL_LGCCgas_SIDIS 
  ////cher_lg_Tank_window_back | cher_lg_Tank |  tank window back |   0*cm 0*cm 300.995*cm |   0*deg 0*deg 0*deg | FF9900 |  Tube |  85.1*cm 264.9*cm 0.005*cm 0*deg 360*deg |  SL_PVF

  //// ---------------------------------------------
  //// cher_lg_tank2
  ////
  //rInner = {
  //  58*cm,65*cm
  //};
  //rOuter = { 
  //  127*cm,144*cm
  //};
  //zPlane = { 
  //  97*cm,194*cm
  //};
  //Polycone cher_lg_tank2( 0,  ///< Initial Phi starting angle
  //    360*deg,  ///< Total Phi angle
  //    rInner,        ///< Tangent distance to outer surface
  //    rOuter,        ///< Tangent distance to outer surface
  //    zPlane);       ///< z coordinate of corners
  //Position cher_lg_tank2_pos = {0*cm,0.0*cm,0.0*cm};
  //Volume cher_lg_tank2_vol("cher_lg_tank2_vol", cher_lg_tank2, air);
  //cher_lg_tank2_vol.setVisAttributes(lcdd, "GreenVis");
  //pv = assembly.placeVolume(cher_lg_tank2_vol, cher_lg_tank2_pos);
  //
  // cher_lg_Tank2 |  root |   tank |   0*cm 0*cm 0*cm |   0*deg 0*deg 0*deg | FF9900 | Polycone |0*deg 360*deg 2*counts 58*cm 65*cm 127*cm 144*cm 97*cm 194*cm | SL_LGCCgas_SIDIS |  no | 1 | 1 | 1 | 1 | 1 |  no |  no |    no 
  //
  //cher_lg_Tank2_window_front | cher_lg_Tank2 |  tank window front |   0*cm 0*cm 97.0025*cm |   0*deg 0*deg 0*deg | FF9900 |  Tube |  58.1*cm 126.9*cm 0.0025*cm 0*deg 360*deg |  SL_PVF |  no | 1 | 1 | 1 | 1 | 1 |  no |  no |    no 




  //// ---------------------------------------------
  //// Heavy Gas Chrenkov
  //// ---------------------------------------------
  //// cher_hg_chamber
  ////
  //rInner = {
  //  80*cm,94*cm
  //};
  //rOuter = { 
  //  265*cm,265*cm
  //};
  //zPlane = { 
  //  306*cm,406*cm
  //};
  //Polycone cher_hg_chamber( 0,  ///< Initial Phi starting angle
  //    360*deg,  ///< Total Phi angle
  //    rInner,        ///< Tangent distance to outer surface
  //    rOuter,        ///< Tangent distance to outer surface
  //    zPlane);       ///< z coordinate of corners
  //Position cher_hg_chamber_pos = {0*cm,0.0*cm,0.0*cm};
  //Volume cher_hg_chamber_vol("cher_hg_chamber_vol", cher_hg_chamber, air);
  //cher_hg_chamber_vol.setVisAttributes(lcdd, "PurpleVis");
  //pv = assembly.placeVolume(cher_hg_chamber_vol, cher_hg_chamber_pos);

  //// ---------------------------------------------
  //// cher_hg_gas
  ////
  //rInner = {
  //  80.5*cm,94.5*cm
  //};
  //rOuter = { 
  //  264.5*cm,264.5*cm
  //};
  //zPlane = { 
  //  306.056*cm,405.9*cm
  //};
  //Polycone cher_hg_gas( 0,  ///< Initial Phi starting angle
  //    360*deg,  ///< Total Phi angle
  //    rInner,        ///< Tangent distance to outer surface
  //    rOuter,        ///< Tangent distance to outer surface
  //    zPlane);       ///< z coordinate of corners
  //Position cher_hg_gas_pos = {0*cm,0.0*cm,0.0*cm};
  //Volume cher_hg_gas_vol("cher_hg_gas_vol", cher_hg_gas, air);
  //cher_hg_gas_vol.setVisAttributes(lcdd, "PurpleVis");
  //pv = assembly.placeVolume(cher_hg_gas_vol, cher_hg_gas_pos);

  //______________________________________________________________________________

  pv = lcdd.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system",sdet.id()).addPhysVolID("barrel",0);
  sdet.setPlacement(pv);
  return sdet;

}

DECLARE_DETELEMENT(StandInGeometrySOLID,create_detector)
