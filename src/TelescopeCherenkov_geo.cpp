#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <XML/Helper.h>
#include "TMath.h"

using namespace dd4hep;
using namespace dd4hep::rec;

static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement x_det    = handle;
    std::string     det_name = x_det.nameStr();
    int             det_id   = x_det.id();
    DetElement      det(det_name, det_id);
    sens.setType("photoncounter");

    xml_dim_t x_pl = x_det.child(_U(placement));

    // build tank
    xml_dim_t x_tank = x_det.child(_Unicode(tank));
    auto tank_mat     = desc.material(x_tank.attr<std::string>(_Unicode(material)));
    auto rmin         = x_tank.rmin();
    auto rmax1        = x_tank.rmax1();
    auto rmax2        = x_tank.rmax2();
    auto length1      = x_tank.attr<double>(_Unicode(length1));
    auto length2      = x_tank.attr<double>(_Unicode(length2));
    auto vert_height  = x_tank.attr<double>(_Unicode(vert_height));
    auto flange_width = x_tank.attr<double>(_Unicode(flange_width));
    auto inter_shift  = x_tank.attr<double>(_Unicode(intersect_shift));

    // main tube
    Tube tank_hori(0., rmax2, length1/2.);
    Tube tank_vert(0., rmax2, vert_height/2.);
    UnionSolid tank_part1(tank_hori, tank_vert, Transform3D(Position(0., vert_height/2., inter_shift))*RotationX(M_PI/2.));

    // subtract the inner part
    Tube tank_hori2(0., rmin, length1/2.);
    Tube tank_vert2(0., rmin, vert_height/2.);
    UnionSolid tank_part2(tank_hori2, tank_vert2, Transform3D(Position(0., vert_height/2., inter_shift))*RotationX(M_PI/2.));
    SubtractionSolid tank_union(tank_part1, tank_part2);

    // subtract the outer tube
    double out1_length = (length1 - length2)/2. - inter_shift - flange_width;
    double out2_length = (length1 - length2)/2. + inter_shift - flange_width;
    Tube tank_out1(rmax1, rmax2, out1_length/2.);
    Tube tank_out2(rmax1, rmax2, out2_length/2.);
    SubtractionSolid tank_sub1(tank_union, tank_out1, Position(0., 0., (length1 - out1_length - flange_width)/2.));
    SubtractionSolid tank_solid(tank_sub1, tank_out2, Position(0., 0., -(length1 - out2_length - flange_width)/2.));
    // tank
    Volume v_tank("vol_tcd_tank", tank_solid, tank_mat);
    v_tank.setVisAttributes(desc, x_tank.attr<std::string>(_Unicode(vis)));
    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(v_tank, Position(x_pl.x(), x_pl.y(), x_pl.z()));
    envPV.addPhysVolID("system", det_id);
    det.setPlacement(envPV);

    return det;
}
//@}
// clang-format off
DECLARE_DETELEMENT(SoLID_TelescopeCherenkov, createDetector)

