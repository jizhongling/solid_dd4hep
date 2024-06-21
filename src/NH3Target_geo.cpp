// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Whitney Armstrong

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"

using namespace std;
using namespace dd4hep;


static Ref_t create_detector(Detector& lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t  x_det     = e;
  string     det_name  = x_det.nameStr();
  Material   air       = lcdd.air();

  DetElement sdet(det_name,x_det.id());
  Assembly   assembly("NH3Target");

  double z_center = -300.0*cm;

  // ---------------------------------------------
  // Target tube
  //
  Tube target_tube( 0*cm,2.5*cm,5.0*cm);
  Volume v_target_tube("v_target_tube", target_tube, air);
  v_target_tube.setVisAttributes(lcdd, "RedVis");
  auto pv_target_tube = assembly.placeVolume(v_target_tube);

  //______________________________________________________________________________

  auto pv = lcdd.pickMotherVolume(sdet).placeVolume(assembly, Position(0.0,0.0,z_center));
  pv.addPhysVolID("system",sdet.id()).addPhysVolID("layer",1);
  sdet.setPlacement(pv);
  return sdet;

}

DECLARE_DETELEMENT(NH3Target,create_detector)
