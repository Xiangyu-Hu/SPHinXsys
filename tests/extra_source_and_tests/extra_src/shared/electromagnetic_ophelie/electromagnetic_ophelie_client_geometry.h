/**
 * @file electromagnetic_ophelie_client_geometry.h
 * @brief Client cold-crucible assembly metadata (Doc software plan + Doc228/229).
 *
 * CAD STLs live in the corresponding relax case data/ folders (mm, +Y up):
 *   test_3d_ophelie_client_nc_glass_relax/data/
 *   test_3d_ophelie_client_stirring_glass_relax/data/
 * Recommended solver frame: translate crucible outer-bottom axis to origin, then
 * (x,y,z)_cad → (x, -z, y)_solver / 1000 → metres, +Z up.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_CLIENT_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_CLIENT_GEOMETRY_H

#include "sphinxsys.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OphelieClientAssemblyKind
{
    NaturalConvection = 0,
    MechanicallyStirred = 1
};

struct OphelieClientVec3
{
    Real x = 0;
    Real y = 0;
    Real z = 0;
};

struct OphelieClientAabbMm
{
    OphelieClientVec3 min_mm{};
    OphelieClientVec3 max_mm{};
};

inline OphelieClientVec3 clientCadMmToSolverM(const OphelieClientVec3 &cad_mm_translated)
{
    // (x,y,z)_cad → (x, -z, y) then mm→m
    OphelieClientVec3 s;
    s.x = cad_mm_translated.x * Real(1.0e-3);
    s.y = -cad_mm_translated.z * Real(1.0e-3);
    s.z = cad_mm_translated.y * Real(1.0e-3);
    return s;
}

struct OphelieClientCoilRingMm
{
    int index = 0;
    Real z_lo_mm = 0; // axial interval along CAD +Y (outer-bottom origin)
    Real z_hi_mm = 0;
    Real center_height_mm = 0;
    Real r_inner_mm = Real(168.5);
    Real r_outer_mm = Real(171.5);
    Real r_center_mm = Real(170);
};

inline std::array<OphelieClientCoilRingMm, 4> clientCoilRingsRelativeToOuterBottom()
{
    // Software plan §4.1; relative to crucible outer bottom = 0 on CAD vertical.
    return {{
        {1, Real(39.078), Real(119.078), Real(79.078), Real(168.5), Real(171.5), Real(170)},
        {2, Real(129.078), Real(209.078), Real(169.078), Real(168.5), Real(171.5), Real(170)},
        {3, Real(220), Real(300), Real(260), Real(168.5), Real(171.5), Real(170)},
        {4, Real(310), Real(390), Real(350), Real(168.5), Real(171.5), Real(170)},
    }};
}

struct OphelieClientPartRole
{
    const char *role = "";
    const char *relative_stl = "";
    bool enters_solver_bodies = false;
    bool probe_or_display_only = true;
};

struct OphelieClientAssemblyCatalog
{
    OphelieClientAssemblyKind kind = OphelieClientAssemblyKind::NaturalConvection;
    const char *assembly_id = "";
    const char *stl_subdir = "";
    /** Measured glass AABB in raw CAD mm (from binary STL audit 2026-11-20). */
    OphelieClientAabbMm glass_aabb_raw_mm{};
    /** Thermocouple AABB in raw CAD mm. */
    OphelieClientAabbMm thermocouple_aabb_raw_mm{};
    /** Paddle AABB in raw CAD mm; NC unused. */
    OphelieClientAabbMm paddle_aabb_raw_mm{};
    bool has_paddle = false;
    /** Historical eccentricity of paddle axis vs crucible axis in CAD X (mm). */
    Real paddle_eccentricity_x_mm = 0;
    std::vector<OphelieClientPartRole> parts;
};

inline OphelieClientAssemblyCatalog makeClientNaturalConvectionAssemblyCatalog()
{
    OphelieClientAssemblyCatalog c;
    c.kind = OphelieClientAssemblyKind::NaturalConvection;
    c.assembly_id = "client_nc_2026-11";
    c.stl_subdir =
        "tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/data";
    c.glass_aabb_raw_mm = {{Real(30.715), Real(0), Real(30.715)}, {Real(315.715), Real(250), Real(315.715)}};
    c.thermocouple_aabb_raw_mm = {{Real(166.227), Real(130), Real(166.232)},
                                  {Real(180.212), Real(400), Real(180.214)}};
    c.has_paddle = false;
    c.parts = {
        {"crucible_shell", "glass_NC.STL", false, true},
        {"coil", "coil_NC.STL", false, true},
        {"cooling_pipes", "cooling-pipes_NC.STL", false, true},
        {"thermocouple", "thermal-couple_NC.STL", false, true},
    };
    return c;
}

inline OphelieClientAssemblyCatalog makeClientStirringAssemblyCatalog()
{
    OphelieClientAssemblyCatalog c;
    c.kind = OphelieClientAssemblyKind::MechanicallyStirred;
    c.assembly_id = "client_stir_2026-11";
    c.stl_subdir =
        "tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/data";
    c.glass_aabb_raw_mm = {{Real(32.739), Real(3.732), Real(1029)}, {Real(317.739), Real(253.732), Real(1314)}};
    c.thermocouple_aabb_raw_mm = {{Real(238.739), Real(133.732), Real(1164.506)},
                                  {Real(252.714), Real(403.732), Real(1178.494)}};
    c.paddle_aabb_raw_mm = {{Real(43.265), Real(183.732), Real(1126.661)},
                            {Real(173.213), Real(683.732), Real(1216.339)}};
    c.has_paddle = true;
    c.paddle_eccentricity_x_mm = Real(-67); // audit / plan; verify before hard-coding drive axis
    c.parts = {
        {"crucible_shell", "glass.STL", false, true},
        {"coil", "coil.STL", false, true},
        {"cooling_pipes", "cooling-pipes-1.STL", false, true},
        {"thermocouple", "thermal-couple-1.STL", false, true},
        {"stirring_paddle", "stirring paddle.STL", true, false},
    };
    return c;
}

struct OphelieClientCrucibleMetricsM
{
    Real outer_diameter_m = Real(0.285);
    Real outer_height_m = Real(0.250);
    Real inner_diameter_m = Real(0.283); // plan estimate; confirm before production
    Real inner_bottom_above_outer_m = Real(0.001);
    Real melt_radius_m = Real(0.1415); // 0.5 * inner diameter estimate
};

inline OphelieClientCrucibleMetricsM clientCrucibleMetricsFromPlan()
{
    return {};
}

/**
 * Solver-frame (+Z up) parametric melt used for particle generation.
 * Fluid is NOT filled from the shell STL; shell remains display-only.
 *
 * melt_height_m is a DEVELOPMENT default until the client confirms liquid depth.
 */
struct OphelieClientMeltGeometryM
{
    Real melt_radius_m = Real(0.1415);
    Real melt_height_m = Real(0.180); // provisional; client liquid depth TBD
    Real inner_bottom_z_m = Real(0.001);
    Real glass_center_z_m = Real(0.001 + 0.5 * 0.180);
    Vecd glass_center = Vecd(0, 0, 0.091);
    Real glass_half_height = Real(0.09);
    Real coil_loop_radius_m = Real(0.170);
    int coil_num_loops = 4;
    Real coil_z_min_m = Real(0.079078); // ring-1 centre
    Real coil_z_max_m = Real(0.350);    // ring-4 centre
    const char *melt_height_status = "provisional_development_liquid_depth";
};

inline OphelieClientMeltGeometryM makeClientNcMeltGeometryDevelopment()
{
    OphelieClientMeltGeometryM g;
    const auto m = clientCrucibleMetricsFromPlan();
    g.melt_radius_m = m.melt_radius_m;
    g.melt_height_m = Real(0.180);
    g.inner_bottom_z_m = m.inner_bottom_above_outer_m;
    g.glass_half_height = Real(0.5) * g.melt_height_m;
    g.glass_center_z_m = g.inner_bottom_z_m + g.glass_half_height;
    g.glass_center = Vecd(0, 0, g.glass_center_z_m);
    const auto rings = clientCoilRingsRelativeToOuterBottom();
    g.coil_loop_radius_m = rings[0].r_center_mm * Real(1.0e-3);
    g.coil_num_loops = 4;
    g.coil_z_min_m = rings[0].center_height_mm * Real(1.0e-3);
    g.coil_z_max_m = rings[3].center_height_mm * Real(1.0e-3);
    return g;
}

/** Parametric melt cylinder (lattice / relax), matching French natural style. */
class OphelieClientNcGlassCylinderShape : public ComplexShape
{
  public:
    OphelieClientNcGlassCylinderShape(const std::string &shape_name, const OphelieClientMeltGeometryM &melt,
                                      int mesh_resolution = 20)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), melt.melt_radius_m, melt.glass_half_height, mesh_resolution,
                                       melt.glass_center, "ClientNcMelt");
    }
};

/** Open-top numerical wall: annular side + bottom (display shell STL is separate). */
class OphelieClientNcCrucibleWallShape : public ComplexShape
{
  public:
    OphelieClientNcCrucibleWallShape(const std::string &shape_name, const OphelieClientMeltGeometryM &melt,
                                     Real wall_thickness, int mesh_resolution = 20, Real rim_height = 0.0)
        : ComplexShape(shape_name)
    {
        const Vecd &c = melt.glass_center;
        const Real r_in = melt.melt_radius_m;
        const Real r_out = melt.melt_radius_m + wall_thickness;
        const Real half_h = melt.glass_half_height;
        const Real rim = std::max(rim_height, Real(0));
        const Vecd outer_center(c[0], c[1], c[2] + Real(0.5) * rim - Real(0.5) * wall_thickness);
        const Real outer_half_h = half_h + Real(0.5) * rim + Real(0.5) * wall_thickness;
        const Vecd inner_center(c[0], c[1], c[2] + Real(0.5) * rim);
        const Real inner_half_h = half_h + Real(0.5) * rim;
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_out, outer_half_h, mesh_resolution, outer_center,
                                       "OuterBoundary");
        subtract<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_in, inner_half_h, mesh_resolution, inner_center,
                                            "InnerBoundary");
    }
};

inline BoundingBoxd clientNcDomainBounds(const OphelieClientMeltGeometryM &melt, Real boundary_width)
{
    const Real radial = std::max(melt.melt_radius_m, melt.coil_loop_radius_m) + boundary_width;
    const Real z_lo = std::min(melt.inner_bottom_z_m, melt.coil_z_min_m) - boundary_width;
    const Real z_hi = std::max(melt.inner_bottom_z_m + melt.melt_height_m, melt.coil_z_max_m) + boundary_width;
    return BoundingBoxd(Vecd(-radial, -radial, z_lo), Vecd(radial, radial, z_hi));
}

struct OphelieClientProbeDefinition
{
    const char *probe_id = "thermocouple_nominal_tip";
    /** Status: unconfirmed_sensor_location until client confirms active junction. */
    const char *status = "unconfirmed_sensor_location";
    OphelieClientVec3 cad_tip_mm_raw{};
    OphelieClientVec3 solver_m{};
    Real target_temperature_k = 0;
    const char *target_meaning = "";
};

/**
 * Nominal tip = lower end of thermocouple AABB (CAD vertical = Y).
 * NOT confirmed as the true RTD/junction location.
 */
inline OphelieClientProbeDefinition makeClientNominalThermocoupleProbe(const OphelieClientAssemblyCatalog &cat)
{
    OphelieClientProbeDefinition p;
    const auto &bb = cat.thermocouple_aabb_raw_mm;
    p.cad_tip_mm_raw = {Real(0.5) * (bb.min_mm.x + bb.max_mm.x), bb.min_mm.y,
                        Real(0.5) * (bb.min_mm.z + bb.max_mm.z)};

    // Translate so crucible outer-bottom axis → origin (horizontal mid of glass AABB, y = glass ymin).
    const auto &g = cat.glass_aabb_raw_mm;
    const OphelieClientVec3 origin_cad = {Real(0.5) * (g.min_mm.x + g.max_mm.x), g.min_mm.y,
                                          Real(0.5) * (g.min_mm.z + g.max_mm.z)};
    OphelieClientVec3 tip_t = {p.cad_tip_mm_raw.x - origin_cad.x, p.cad_tip_mm_raw.y - origin_cad.y,
                               p.cad_tip_mm_raw.z - origin_cad.z};
    p.solver_m = clientCadMmToSolverM(tip_t);

    if (cat.kind == OphelieClientAssemblyKind::NaturalConvection)
    {
        p.target_temperature_k = Real(1273.15);
        p.target_meaning = "NC thermocouple ~1000C target (not domain max / not IC)";
    }
    else
    {
        p.target_temperature_k = Real(1673.15);
        p.target_meaning = "stirring ~1400C qualitative peak trend only; time/location TBD";
    }
    return p;
}

inline bool writeClientNcMeltGeometryJson(const OphelieClientMeltGeometryM &melt, const OphelieClientProbeDefinition &probe,
                                          const std::string &path)
{
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "{\n";
    out << "  \"melt_radius_m\": " << melt.melt_radius_m << ",\n";
    out << "  \"melt_height_m\": " << melt.melt_height_m << ",\n";
    out << "  \"melt_height_status\": \"" << melt.melt_height_status << "\",\n";
    out << "  \"inner_bottom_z_m\": " << melt.inner_bottom_z_m << ",\n";
    out << "  \"glass_center_m\": [" << melt.glass_center[0] << "," << melt.glass_center[1] << ","
        << melt.glass_center[2] << "],\n";
    out << "  \"coil_loop_radius_m\": " << melt.coil_loop_radius_m << ",\n";
    out << "  \"coil_num_loops\": " << melt.coil_num_loops << ",\n";
    out << "  \"coil_z_min_m\": " << melt.coil_z_min_m << ",\n";
    out << "  \"coil_z_max_m\": " << melt.coil_z_max_m << ",\n";
    out << "  \"probe_id\": \"" << probe.probe_id << "\",\n";
    out << "  \"probe_status\": \"" << probe.status << "\",\n";
    out << "  \"probe_solver_m\": [" << probe.solver_m.x << "," << probe.solver_m.y << "," << probe.solver_m.z
        << "],\n";
    out << "  \"probe_target_k\": " << probe.target_temperature_k << ",\n";
    out << "  \"note\": \"parametric melt; shell STL display-only; not French R=0.25 geometry\"\n";
    out << "}\n";
    return true;
}

inline bool writeClientAssemblyMetadataJson(const OphelieClientAssemblyCatalog &cat, const std::string &path)
{
    const auto probe = makeClientNominalThermocoupleProbe(cat);
    const auto rings = clientCoilRingsRelativeToOuterBottom();
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "{\n";
    out << "  \"assembly_id\": \"" << cat.assembly_id << "\",\n";
    out << "  \"stl_subdir\": \"" << cat.stl_subdir << "\",\n";
    out << "  \"cad_vertical\": \"+Y\",\n";
    out << "  \"stl_unit\": \"mm\",\n";
    out << "  \"solver_up\": \"+Z\",\n";
    out << "  \"transform\": \"translate crucible outer-bottom axis to origin; (x,y,z)_cad -> (x,-z,y)/1000\",\n";
    out << "  \"has_paddle\": " << (cat.has_paddle ? "true" : "false") << ",\n";
    out << "  \"paddle_eccentricity_x_mm_historical\": " << cat.paddle_eccentricity_x_mm << ",\n";
    out << "  \"glass_aabb_raw_mm\": {"
        << "\"min\":[" << cat.glass_aabb_raw_mm.min_mm.x << "," << cat.glass_aabb_raw_mm.min_mm.y << ","
        << cat.glass_aabb_raw_mm.min_mm.z << "],"
        << "\"max\":[" << cat.glass_aabb_raw_mm.max_mm.x << "," << cat.glass_aabb_raw_mm.max_mm.y << ","
        << cat.glass_aabb_raw_mm.max_mm.z << "]},\n";
    out << "  \"probe\": {"
        << "\"id\":\"" << probe.probe_id << "\","
        << "\"status\":\"" << probe.status << "\","
        << "\"cad_tip_mm_raw\":[" << probe.cad_tip_mm_raw.x << "," << probe.cad_tip_mm_raw.y << ","
        << probe.cad_tip_mm_raw.z << "],"
        << "\"solver_m\":[" << probe.solver_m.x << "," << probe.solver_m.y << "," << probe.solver_m.z << "],"
        << "\"target_temperature_k\":" << probe.target_temperature_k << ","
        << "\"target_meaning\":\"" << probe.target_meaning << "\"},\n";
    out << "  \"coil_rings_mm_relative_outer_bottom\": [\n";
    for (size_t i = 0; i < rings.size(); ++i)
    {
        const auto &r = rings[i];
        out << "    {\"index\":" << r.index << ",\"z_lo\":" << r.z_lo_mm << ",\"z_hi\":" << r.z_hi_mm
            << ",\"center\":" << r.center_height_mm << ",\"r_center\":" << r.r_center_mm << "}"
            << (i + 1 < rings.size() ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"parts\": [\n";
    for (size_t i = 0; i < cat.parts.size(); ++i)
    {
        const auto &p = cat.parts[i];
        out << "    {\"role\":\"" << p.role << "\",\"file\":\"" << p.relative_stl
            << "\",\"enters_solver\":" << (p.enters_solver_bodies ? "true" : "false")
            << ",\"display_or_probe_only\":" << (p.probe_or_display_only ? "true" : "false") << "}"
            << (i + 1 < cat.parts.size() ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"notes\": ["
        << "\"cooling pipes are display-only\","
        << "\"thermocouple tip is nominal AABB lower end until client confirms\","
        << "\"coil/cooling STL may contain zero-area triangles; keep originals\""
        << "]\n";
    out << "}\n";
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLIENT_GEOMETRY_H
