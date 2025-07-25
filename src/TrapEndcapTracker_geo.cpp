// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Whitney Armstrong, Zhongling Ji

/** \addtogroup Trackers Trackers
 * \brief Type: **BarrelTrackerWithFrame**.
 * \author W. Armstrong
 *
 * \ingroup trackers
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DD4hepDetectorHelper.h"
#include <array>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

/** Endcap Trapezoidal Tracker.
 *
 * @author Whitney Armstrong
 *
 */
static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t x_det = e;
  Material vacuum = description.vacuum();
  int det_id      = x_det.id();
  string det_name = x_det.nameStr();
  bool reflect    = x_det.reflect(false);
  DetElement sdet(det_name, det_id);
  Assembly assembly(det_name);

  Material air     = description.material("Air");
  Volume motherVol = description.pickMotherVolume(sdet);
  int m_id = 0, c_id = 0, n_sensor = 0;
  map<string, Volume> modules;
  map<string, Placements> sensitives;
  map<string, std::vector<VolPlane>> volplane_surfaces;
  map<string, std::array<double, 2>> module_thicknesses;
  PlacedVolume pv;

  // Set detector type flag
  dd4hep::xml::setDetectorTypeFlag(x_det, sdet);
  auto& params = DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(sdet);

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
        "boundary_material");
  }

  assembly.setVisAttributes(description.invisible());
  sens.setType("tracker");

  for (xml_coll_t su(x_det, _U(support)); su; ++su) {
    xml_comp_t x_support     = su;
    double support_thickness = getAttrOrDefault(x_support, _U(thickness), 2.0 * mm);
    double support_length    = getAttrOrDefault(x_support, _U(length), 2.0 * mm);
    double support_rmin      = getAttrOrDefault(x_support, _U(rmin), 2.0 * mm);
    double support_zstart    = getAttrOrDefault(x_support, _U(zstart), 2.0 * mm);
    std::string support_name =
      getAttrOrDefault<std::string>(x_support, _Unicode(name), "support_tube");
    std::string support_vis = getAttrOrDefault<std::string>(x_support, _Unicode(vis), "AnlRed");
    xml_dim_t pos(x_support.child(_U(position), false));
    xml_dim_t rot(x_support.child(_U(rotation), false));
    Solid support_solid;
    if (x_support.hasChild(_U(shape))) {
      xml_comp_t shape(x_support.child(_U(shape)));
      string shape_type = shape.typeStr();
      support_solid     = xml::createShape(description, shape_type, shape);
    } else {
      support_solid = Tube(support_rmin, support_rmin + support_thickness, support_length / 2);
    }
    Transform3D tr =
      Transform3D(Rotation3D(),
          Position(0, 0, (reflect ? -1.0 : 1.0) * (support_zstart + support_length / 2)));
    if (pos.ptr() && rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      Position pos3D(pos.x(0), pos.y(0), pos.z(0));
      tr = Transform3D(rot3D, pos3D);
    } else if (pos.ptr()) {
      tr = Transform3D(Rotation3D(), Position(pos.x(0), pos.y(0), pos.z(0)));
    } else if (rot.ptr()) {
      Rotation3D rot3D(RotationZYX(rot.z(0), rot.y(0), rot.x(0)));
      tr = Transform3D(rot3D, Position());
    }
    Material support_mat = description.material(x_support.materialStr());
    Volume support_vol(support_name, support_solid, support_mat);
    support_vol.setVisAttributes(description.visAttributes(support_vis));
    pv = assembly.placeVolume(support_vol, tr);
  }

  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer(li);
    int l_id    = x_layer.id();
    int mod_num = 1;

    string layer_name     = det_name + std::string("_layer") + std::to_string(l_id);
    string layer_modnam   = x_layer.child(_U(module)).attr<string>("name");

    std::string layer_vis = x_layer.attr<std::string>(_Unicode(vis));
    double layer_rmin     = x_layer.attr<double>(_Unicode(rmin));
    double layer_rmax     = x_layer.attr<double>(_Unicode(rmax));
    int nmodules          = x_layer.attr<int>(_Unicode(nmodules));
    double layer_phi0     = x_layer.phi0(0);
    double layer_length   = 0;
    double layer_zstart   = x_layer.attr<double>(_Unicode(zstart));

    for (xml_coll_t mi(x_det, _U(module)); mi; ++mi, ++m_id) {
      xml_comp_t x_mod = mi;
      if(layer_modnam != x_mod.nameStr())
        continue;

      string m_nam           = layer_modnam + string("_layer") + to_string(l_id);
      double ang             = M_PI / nmodules;
      double x1              = layer_rmin * tan(ang);
      double x2              = layer_rmax * sin(ang);
      double z               = (layer_rmax * cos(ang) - layer_rmin) / 2;
      double posY;
      double total_thickness = 0.;
      xml_coll_t ci(x_mod, _U(module_component));
      for (ci.reset(), total_thickness = 0.0; ci; ++ci)
        total_thickness += xml_comp_t(ci).thickness();
      layer_length += total_thickness;

      double thickness_so_far = 0.0;
      double y1               = total_thickness / 2;
      double y2               = total_thickness / 2;
      Trapezoid m_solid(x1, x2, y1, y2, z);
      Volume m_volume(m_nam, m_solid, vacuum);
      m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));

      Solid frame_s;
      if (x_mod.hasChild(_U(frame))) {
        // build frame from trd (assumed to be smaller)
        xml_comp_t m_frame     = x_mod.child(_U(frame));
        xml_comp_t f_pos       = m_frame.child(_U(position));
        xml_comp_t frame_trd   = m_frame.trd();
        double frame_thickness = getAttrOrDefault(m_frame, _U(thickness), total_thickness);
        double frame_x1        = frame_trd.x1();
        double frame_x2        = frame_trd.x2();
        double frame_z         = frame_trd.z();
        // make the frame match the total thickness if thickness attribute is not given
        Trapezoid f_solid1(x1, x2, frame_thickness / 2.0, frame_thickness / 2.0, z);
        Trapezoid f_solid(frame_x1, frame_x2, frame_thickness / 2.0, frame_thickness / 2.0, frame_z);
        SubtractionSolid frame_shape(f_solid1, f_solid);
        frame_s = frame_shape;

        Material f_mat = description.material(m_frame.materialStr());
        Volume f_vol(m_nam + "_frame", frame_shape, f_mat);
        f_vol.setVisAttributes(description.visAttributes(m_frame.visStr()));

        // figure out how to best place
        pv = m_volume.placeVolume(f_vol, Position(f_pos.x(), f_pos.y(), f_pos.z()));
      }

      for (ci.reset(), n_sensor = 1, c_id = 0, posY = -y1; ci; ++ci, ++c_id) {
        xml_comp_t c     = ci;
        double c_thick   = c.thickness();
        auto comp_x1     = getAttrOrDefault(c, _Unicode(x1), x1);
        auto comp_x2     = getAttrOrDefault(c, _Unicode(x2), x2);
        auto comp_height = getAttrOrDefault(c, _Unicode(height), z);

        Material c_mat = description.material(c.materialStr());
        string c_name  = _toString(c_id, "component%d");

        Trapezoid comp_s1(comp_x1, comp_x2, c_thick / 2e0, c_thick / 2e0, comp_height);
        Solid comp_shape = comp_s1;
        if (frame_s.isValid()) {
          comp_shape = SubtractionSolid(comp_s1, frame_s);
        }
        Volume c_vol(c_name, comp_shape, c_mat);

        c_vol.setVisAttributes(description.visAttributes(c.visStr()));
        pv = m_volume.placeVolume(c_vol, Position(0, posY + c_thick / 2, 0));
        if (c.isSensitive()) {
          module_thicknesses[m_nam] = {thickness_so_far + c_thick / 2.0,
            total_thickness - thickness_so_far - c_thick / 2.0};
          sdet.check(n_sensor > 2,
              "TrapEndcapTracker::fromCompact: " + c_name + " Max of 2 modules allowed!");
          pv.addPhysVolID("sensor", n_sensor);
          c_vol.setSensitiveDetector(sens);
          sensitives[m_nam].push_back(pv);
          ++n_sensor;
          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u(0., 0., -1.);
          Vector3D v(-1., 0., 0.);
          Vector3D n(0., 1., 0.);
          // Vector3D o( 0. , 0. , 0. ) ;

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = module_thicknesses[m_nam][0];
          double outer_thickness = module_thicknesses[m_nam][1];

          SurfaceType type(SurfaceType::Sensitive);

          VolPlane surf(c_vol, type, inner_thickness, outer_thickness, u, v, n); //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

          //--------------------------------------------
        }
        posY += c_thick;
        thickness_so_far += c_thick;
      }
      modules[m_nam] = m_volume;
    }

    double layer_center_z = layer_zstart + layer_length / 2.0;
    Tube layer_tub(layer_rmin, layer_rmax, layer_length / 2);
    Volume layer_vol(layer_name, layer_tub, air); // Create the layer envelope volume.
    layer_vol.setVisAttributes(description.visAttributes(layer_vis));

    PlacedVolume layer_pv;
    if (reflect) {
      layer_pv = assembly.placeVolume(
          layer_vol, Transform3D(RotationZYX(0.0, -M_PI, 0.0), Position(0, 0, -layer_center_z)));
      layer_pv.addPhysVolID("layer", l_id);
    } else {
      layer_pv = assembly.placeVolume(layer_vol, Position(0, 0, layer_center_z));
      layer_pv.addPhysVolID("layer", l_id);
    }
    DetElement layer_element(sdet, layer_name, l_id);
    layer_element.setPlacement(layer_pv);

    auto& layerParams =
      DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(layer_element);

    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      DD4hepDetectorHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
          "layer_material");
    }

    { // ring
      double ang           = M_PI / nmodules;
      double r             = layer_rmin + (layer_rmax * cos(ang) - layer_rmin) / 2;
      double phi0          = layer_phi0;
      double zstart        = 0;
      double dz            = 0;
      string m_nam         = layer_modnam + string("_layer") + to_string(l_id);
      Volume m_vol         = modules[m_nam];
      double iphi          = ang * 2;
      double phi           = phi0;
      Placements& sensVols = sensitives[m_nam];

      for (int k = 0; k < nmodules; ++k) {
        string m_base = _toString(l_id, "layer%d") + _toString(mod_num, "_module%d");
        double x      = -r * std::cos(phi);
        double y      = -r * std::sin(phi);

        if (!reflect) {
          DetElement module(layer_element, m_base + "_pos", det_id);
          pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
                Position(x, y, zstart + dz)));
          pv.addPhysVolID("module", mod_num);
          module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        } else {
          pv = layer_vol.placeVolume(m_vol, Transform3D(RotationZYX(0, -M_PI / 2 - phi, -M_PI / 2),
                Position(x, y, -zstart - dz)));
          pv.addPhysVolID("module", mod_num);
          DetElement r_module(layer_element, m_base + "_neg", det_id);
          r_module.setPlacement(pv);
          for (size_t ic = 0; ic < sensVols.size(); ++ic) {
            PlacedVolume sens_pv = sensVols[ic];
            DetElement comp_elt(r_module, sens_pv.volume().name(), mod_num);
            auto& comp_elt_params =
              DD4hepDetectorHelper::ensureExtension<dd4hep::rec::VariantParameters>(comp_elt);
            comp_elt_params.set<string>("axis_definitions", "XZY");
            comp_elt.setPlacement(sens_pv);
            volSurfaceList(comp_elt)->push_back(volplane_surfaces[m_nam][ic]);
          }
        }
        dz = -dz;
        phi += iphi;
        ++mod_num;
      }
    }
  }

  pv = motherVol.placeVolume(assembly, Position(0, 0, (reflect ? -1.0e-9 : 1.0e-9)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);
  return sdet;
}

DECLARE_DETELEMENT(SoLID_TrapEndcapTracker, create_detector)
