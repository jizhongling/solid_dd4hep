#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TMath.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
//using namespace DDSurfaces;

static Ref_t create_detector(Detector& lcdd, xml_h e, SensitiveDetector sens)
{
  typedef vector<PlacedVolume> Placements;

  xml_det_t   x_det     = e;
  Material    air       = lcdd.air();
  Material    carbon    = lcdd.material("CarbonFiber");
  Material    silicon   = lcdd.material("SiliconOxide");
  int         det_id    = x_det.id();
  string      det_name  = x_det.nameStr();
  PlacedVolume             pv;

  DetElement  sdet(det_name, det_id);
  Assembly    assembly(det_name+"_assembly");

  sens.setType("tracker");
  string module_name = "GEM";
  
  double thickness = 1.0*dd4hep::cm;

  int N_layers      = 0;

  for(xml_coll_t lay( x_det, _U(layer) ); lay; ++lay, ++N_layers)  {

    xml_comp_t x_layer  = lay;
    double     inner_r     = x_layer.attr<double>(  _Unicode(inner_r) ) ;
    double     outer_r     = x_layer.attr<double>(  _Unicode(outer_r) ) ;
    double     phi0_offset = x_layer.attr<double>(  _Unicode(phi0_offset) ) ;
    double     z           = x_layer.attr<double>(  _Unicode(z) ) ;
    int        layer_id    = x_layer.id();//attr<double>(  _Unicode(z) ) ;
    
    string  layer_name = std::string("gem_layer") + std::to_string(layer_id) ;

    Tube    gem_layer(inner_r, outer_r, thickness/2.0);
    Volume  gem_layer_vol("gem_layer_vol", gem_layer, carbon);


    // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
    Vector3D u( 1. , 0. , 0. ) ;
    Vector3D v( 0. , 1. , 0. ) ;
    Vector3D n( 0. , 0. , 1. ) ;
    Vector3D o( 0. , 0. , 0. ) ;
    double inner_thickness = thickness/2.0;
    double outer_thickness = thickness/2.0;
    SurfaceType type( SurfaceType::Sensitive ) ;
    VolPlane    surf( gem_layer_vol, type, inner_thickness , outer_thickness , u,v,n,o ) ;

    gem_layer_vol.setSensitiveDetector(sens);

    DetElement layer_DE( sdet, _toString(layer_id,"layer%d"), layer_id );

    //Assembly   layer_assembly( layer_name+"_assembly" );
    pv = assembly.placeVolume( gem_layer_vol, Transform3D(RotationZ(phi0_offset),Position(0.0,0.0,z)) );
    pv.addPhysVolID( "layer", layer_id );
    layer_DE.setPlacement(pv);
    //layer_DE.setAttributes(lcdd, layer_assembly, "", "", "SiVertexLayerVis");

  }

  sdet.setAttributes(lcdd, assembly,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  //assembly.setVisAttributes(lcdd.invisible());

  pv = lcdd.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);      // Set the subdetector system ID.
  sdet.setPlacement(pv);

  assembly->GetShape()->ComputeBBox() ;
  return sdet;
}

DECLARE_DETELEMENT(GEMTrackerDiscSOLID,create_detector)
