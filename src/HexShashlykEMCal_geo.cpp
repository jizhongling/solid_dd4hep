//==========================================================================
//  SoLID Electromagnetic Calorimeter Implementation
//--------------------------------------------------------------------------
// Author     : C. Peng
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "TMath.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "Math/Point2D.h"
#include <algorithm>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

typedef ROOT::Math::XYPoint Point;


// checks if a point is inside a unit hexagon
// lside = 1, center at (0, 0), and vertices at (1, 0), (-1, 0)
template<class T>
bool in_hex(const T &p)
{
    // Check length (squared) against inner and outer radius
    double l2 = p.x() * p.x() + p.y() * p.y();
    if (l2 > 1.0) return false;
    if (l2 < 0.75) return true; // (sqrt(3)/2)^2 = 3/4

    // Check against borders
    double py = p.y() * 1.15470053838; // 2/sqrt(3)
    if (py > 1.0 || py < -1.0) return false;

    double px = 0.5 * py + p.x();
    if (px > 1.0 || px < -1.0) return false;

    if ((py - px) > 1.0 || (py - px) < -1.0) return false;

    return true;
}


// recursively fill hexagon in a ring
void add_hex(Point p, std::vector<Point> &res, double lside, double rmin, double rmax, double tol = 1e-6)
{
    // center is outside of the boundaries
    if (p.r()  > (rmax + tol) || p.r() < (rmin - tol)) { return; }

    // already exist
    for (auto &pt : res) {
        if (in_hex((p - pt)/lside)) { return; }
    }

    // add into container
    res.emplace_back(p);

    // recursively check neigbors, 2.*sqrt(3)/2.*lside
    double d = lside*2.0;//1.732050807568877;
    // pi/3
    static double sext = 1.0471975512;
    for (int i = 0; i < 6; ++i) {
        add_hex(Point(p.x() + d*sin(i*sext), p.y() + d*cos(i*sext)), res, lside, rmin, rmax, tol);
    }

}

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    static double tolerance = 0e0;
    Layering layering(handle);
    xml_det_t x_det = handle;
    Material air = desc.air();
    int det_id = x_det.id();
    string det_name = x_det.nameStr();

    xml_comp_t x_dim = x_det.dimensions();
    double rmin = x_dim.rmin();
    double rmax = x_dim.rmax();
    double rmod = x_dim.attr<double>(_Unicode(rmod));
    double rtol = x_dim.attr<double>(_Unicode(rtol));
    double det_z0 = x_dim.z0();
    double len_z = layering.totalThickness();

    double module_gap  = 10.0*2.54*0.001; // 10 mil gap
    double side_length = rmod*2.0/std::sqrt(3.0) - module_gap;

    // detector volume
    DetElement sdet(det_name, det_id);
    Volume motherVol = desc.pickMotherVolume(sdet);
    // Tube tube(rmin - 2.0*rmod, rmax + 2.0*rmod, len_z/2.);
    // Volume envelope(det_name, tube, air);
    Assembly envelope(det_name);
    PlacedVolume env_phv = motherVol.placeVolume(envelope, Position(0., 0., det_z0 + len_z/2.));

    env_phv.addPhysVolID("system", det_id);
    sdet.setPlacement(env_phv);

    sens.setType("calorimeter");

    // a modular volume
    PolyhedraRegular m_hex(6, 0., side_length, len_z);
    Volume mod_vol("module", m_hex, air);
    mod_vol.setVisAttributes(desc.visAttributes(x_det.visStr()));
    DetElement mod_DE("module0", 0);

    // layer start point
    double l_pos_z  = -len_z/2.;
    int l_num = 1;
    // Loop over the sets of layer elements in the detector.
    for (xml_coll_t li(x_det,_U(layer)); li; ++li) {
        xml_comp_t x_layer = li;
        int repeat = x_layer.repeat();
        // Loop over number of repeats for this layer.
        for (int j = 0; j < repeat; j++) {
            string l_name = _toString(l_num,"layer%d");
            double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness.

            Position l_pos(0, 0, l_pos_z + l_thickness/2);      // Position of the layer.
            PolyhedraRegular l_hex(6, 0., side_length, l_thickness);
            Volume l_vol(l_name,l_hex,air);
            DetElement layer_DE(mod_DE, l_name, det_id);

            // Loop over the sublayers or slices for this layer.
            int s_num = 1;
            double s_pos_z = -(l_thickness / 2);
            for (xml_coll_t si(x_layer,_U(slice)); si; ++si) {
                xml_comp_t x_slice = si;
                string s_name = _toString(s_num,"slice%d");
                double s_thick = x_slice.thickness();
                PolyhedraRegular s_hex(6, 0., side_length, s_thick);
                Volume s_vol(s_name, s_hex, desc.material(x_slice.materialStr()));
                DetElement slice(layer_DE, s_name, det_id);

                if (x_slice.isSensitive()) {
                    s_vol.setSensitiveDetector(sens);
                }
                slice.setAttributes(desc, s_vol, x_slice.regionStr(), x_slice.limitsStr(), "InvisibleNoDaughters");
                // s_vol.setVisAttributes(desc.invisible());

                // Slice placement.
                PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick/2));
                slice_phv.addPhysVolID("slice", s_num);
                slice.setPlacement(slice_phv);
                // Increment Z position of slice.
                s_pos_z += s_thick;

                // Increment slice number.
                ++s_num;
            }

            // Set region, limitset, and vis of layer.
            layer_DE.setAttributes(desc, l_vol, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
            // l_vol.setVisAttributes(desc.invisible());

            PlacedVolume layer_phv = mod_vol.placeVolume(l_vol, l_pos);
            layer_phv.addPhysVolID("layer", l_num);
            layer_DE.setPlacement(layer_phv);
            // Increment to next layer Z position.
            l_pos_z += l_thickness;
            ++l_num;
        }
    }


    // automatically fill the ring, start with a seed module centered at (0., 0.)
    std::vector<Point> res;
    int nm = int(rmin/(rmod*sqrt(3.))) + 1;
    add_hex(Point(0., nm*rmod*sqrt(3.)), res, rmod*2.0/std::sqrt(3.0), rmin, rmax, rtol);
    // std::cout << det_name << ": " << res.size() << " modules." << std::endl;

    // sort for a better organized id
    sort(res.begin(), res.end(),
        [](const Point &p1, const Point &p2) {
            if (p1.y() == p2.y()) { return p1.x() < p2.x(); }
            return p1.y() < p2.y();
        });

    int nmod = 0;
    for (auto &p : res) {
        PlacedVolume pv = envelope.placeVolume(mod_vol, Position(p.x(), p.y(), 0.));
        // cout << p.x() << ", " << p.y() << endl;
        pv.addPhysVolID("system", det_id);
        pv.addPhysVolID("module", nmod+1);
        auto amod = (nmod == 0 ? mod_DE : mod_DE.clone("module" + std::to_string(nmod+1), nmod+1));
        amod.setPlacement(pv);
        sdet.add(amod);
        nmod++;
    }

    // Set envelope volume attributes.
    envelope.setAttributes(desc, x_det.regionStr(), x_det.limitsStr(), "InvisibleWithDaughters");
    return sdet;
}

DECLARE_DETELEMENT(SoLID_HexShashlykEMCal, createDetector)

