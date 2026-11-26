/**
 * @file electromagnetic_ophelie_client_stirring_geometry.h
 * @brief Client stirring assembly knobs — same *structure* as French stirring
 *        (GlassBody / Rotor / WallBoundary + Reload), but client STL + materials.
 *
 * Architecture target (Doc228 + user 2026-11-21):
 *   mirror test_3d_ophelie_french_stirring_{glass_relax,em}
 *   Rotor mechanics: elastic_driven (stainless 304) instead of host-spin / Simbody rigid
 *   Glass: parametric melt minus paddle (shell STL display-only, like NC)
 *   Material / power: Doc229 client laws; 30 kW generator input; current unknown
 *
 * CAD STLs + preprocess: test_3d_ophelie_client_stirring_glass_relax/data/
 *   preprocess_client_paddle_stl.py → client_stir_paddle_solver_m.stl (mm +Y → m +Z).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_CLIENT_STIRRING_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_CLIENT_STIRRING_GEOMETRY_H

#include "electromagnetic_ophelie_client_case_params.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_stainless_steel.h"
#include "electromagnetic_ophelie_french_stirring_geometry.h"

#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OphelieClientPaddleMechanics
{
    RigidPrescribed = 0, // French-style host spin / Simbody (regression only)
    ElasticDriven = 1    // Doc228 production path for client stirring
};

struct OphelieClientStirringCaseParams
{
    /** Reuse French stirring scaffolding (coil/relax knobs); override paths below. */
    OphelieFrenchStirringCaseParams french_stirring;

    OphelieClientCaseParams client = makeClientStirringDevelopmentParams();
    OphelieClientStainless304Development steel = makeClientStainless304Development();
    OphelieClientPaddleMechanics paddle_mechanics = OphelieClientPaddleMechanics::ElasticDriven;

    /** Solver-frame paddle STL (metres, +Z up). Empty → must preprocess. */
    std::string rotor_stl_solver_m = "./data/client_stir_paddle_solver_m.stl";

    /** Melt fill (parametric cylinder); shell STL not used as fluid.
     *  Client CAD: paddle blade lowest ≈ 180 mm; blade disk ~40–50 mm thick.
     *  Default free surface 230 mm so blades are immersed (was wrongly 180 mm = surface skim). */
    Real melt_radius_m = Real(0.1415);
    Real melt_height_m = Real(0.230);
    Real melt_bottom_z_m = Real(0.001);

    /** Historical CAD eccentricity; verify before locking drive axis. */
    Real paddle_axis_x_m = Real(-0.067);
    Vecd rotation_axis = Vecd(0, 0, 1);
};

inline void applyClientStirringDevelopmentDefaults(OphelieClientStirringCaseParams &c)
{
    applyFrenchReducedStirringDefaults(c.french_stirring);
    c.french_stirring.rotor_stl_path = c.rotor_stl_solver_m;
    c.french_stirring.geometry_scale = Real(1); // client paddle STL is already metres
    // French glass-z.stl is millimetres (scale 0.001). Leave path empty so summaries / domain
    // helpers cannot accidentally inflate the box to ~500 m (Doc248 C5 bad_alloc).
    c.french_stirring.glass_stl_path.clear();
    c.french_stirring.rotor_translation = Vecd::Zero();
    c.french_stirring.rotation_axis = c.rotation_axis;
    c.french_stirring.rotation_center = Vecd(c.paddle_axis_x_m, 0, Real(0.5) * c.melt_height_m);
    c.french_stirring.rotation_speed_rad_s =
        c.client.stirrer_rpm * Real(2) * Real(3.14159265358979323846) / Real(60); // 10 rpm default
    c.french_stirring.rho_rotor = c.steel.rho_kg_m3;
    c.french_stirring.french.dp = c.client.particle_spacing_m > TinyReal ? c.client.particle_spacing_m
                                                                        : Real(0.005);
    c.french_stirring.french.frequency_hz = c.client.frequency_hz;
    c.french_stirring.french.glass_radius = c.melt_radius_m;
    c.french_stirring.french.glass_half_height = Real(0.5) * c.melt_height_m;
    c.french_stirring.french.glass_center = Vecd(0, 0, c.melt_bottom_z_m + c.french_stirring.french.glass_half_height);
    c.french_stirring.french.require_canonical_bottom_at_zero = false;
    // Coil stack = NC development rings remapped into stir melt band (H=0.230 m).
    {
        OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
        melt.melt_height_m = c.melt_height_m;
        melt.glass_half_height = Real(0.5) * c.melt_height_m;
        melt.glass_center_z_m = c.melt_bottom_z_m + melt.glass_half_height;
        melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);
        adaptClientNcCoilStackToMeltHeight(melt);
        applyClientNcMeltToFrenchReduced(c.french_stirring.french, melt, c.french_stirring.french.dp);
        // applyClientNcMeltToFrenchReduced overwrites glass_* from melt; restore stir fill.
        c.french_stirring.french.glass_radius = c.melt_radius_m;
        c.french_stirring.french.glass_half_height = Real(0.5) * c.melt_height_m;
        c.french_stirring.french.glass_center =
            Vecd(0, 0, c.melt_bottom_z_m + c.french_stirring.french.glass_half_height);
        c.french_stirring.french.frequency_hz = c.client.frequency_hz;
        c.french_stirring.french.require_canonical_bottom_at_zero = false;
    }
}

inline const char *clientPaddleMechanicsName(OphelieClientPaddleMechanics m)
{
    return m == OphelieClientPaddleMechanics::ElasticDriven ? "elastic_driven" : "rigid_prescribed";
}

/**
 * Client stirring thermocouple tip in solver metres (from thermal-couple-1.STL AABB).
 * Same convention as makeClientNominalThermocoupleProbe — tip = CAD AABB lower end;
 * sensor junction unconfirmed until client says otherwise.
 */
inline OphelieClientProbeDefinition clientStirringThermocoupleProbe()
{
    return makeClientNominalThermocoupleProbe(makeClientStirringAssemblyCatalog());
}

inline Vecd clientStirringThermocoupleProbePosition()
{
    const OphelieClientProbeDefinition p = clientStirringThermocoupleProbe();
    return Vecd(p.solver_m.x, p.solver_m.y, p.solver_m.z);
}

/** CAD STL filename under stirring relax data/ (display / audit; probe uses tip AABB). */
inline const char *clientStirringThermocoupleStlName() { return "thermal-couple-1.STL"; }

/** Melt = parametric cylinder minus solver-frame paddle STL (French stirring structure).
 *  SPHinXsys GeometricCylinder is X-axis by default — must use TriangleMeshShapeCylinder(+Z)
 *  like NC / French reduced, otherwise melt appears lying on its side (90° flip vs wall/CAD). */
class OphelieClientStirringGlassFluidShape : public ComplexShape
{
  public:
    OphelieClientStirringGlassFluidShape(const std::string &shape_name, const OphelieClientStirringCaseParams &c,
                                         int mesh_resolution = 20)
        : ComplexShape(shape_name)
    {
        const Vecd center = c.french_stirring.french.glass_center;
        const Real R = c.french_stirring.french.glass_radius;
        const Real half_h = c.french_stirring.french.glass_half_height;
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), R, half_h, mesh_resolution, center, "ClientStirMelt");
        subtract<TriangleMeshShapeSTL>(c.rotor_stl_solver_m, Vecd::Zero(), Real(1), "RotorCutout");
    }
};

class OphelieClientStirringRotorShape : public ComplexShape
{
  public:
    OphelieClientStirringRotorShape(const std::string &shape_name, const OphelieClientStirringCaseParams &c)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(c.rotor_stl_solver_m, Vecd::Zero(), Real(1), "Rotor");
    }
};

inline BoundingBoxd clientStirringRelaxDomainBounds(const OphelieClientStirringCaseParams &c)
{
    const Real pad = Real(8) * c.french_stirring.french.dp;
    TriangleMeshShapeSTL rotor(c.rotor_stl_solver_m, Vecd::Zero(), Real(1), "RotorBounds");
    BoundingBoxd b = rotor.getBounds();
    const Vecd gc = c.french_stirring.french.glass_center;
    const Real R = c.french_stirring.french.glass_radius;
    const Real hh = c.french_stirring.french.glass_half_height;
    b.lower_ = b.lower_.cwiseMin(gc - Vecd(R, R, hh));
    b.upper_ = b.upper_.cwiseMax(gc + Vecd(R, R, hh));
    return b.expand(pad);
}

/**
 * Flow / SPHSystem domain for client stir on french_stirring_em.
 * Must NOT call frenchStirringFlowDomainBounds: that loads French glass-z.stl (mm) with
 * geometry_scale=1 → ~500 m box → cell-list bad_alloc.
 */
inline BoundingBoxd clientStirringFlowDomainBounds(const OphelieFrenchStirringCaseParams &stirring)
{
    const Real dp = stirring.french.dp;
    const Real pad = Real(8) * dp;
    const Real wall = frenchStirringWallThickness(stirring);
    TriangleMeshShapeSTL rotor(stirring.rotor_stl_path, stirring.rotor_translation, Real(1), "ClientRotorBounds");
    BoundingBoxd b = rotor.getBounds();
    const Vecd gc = stirring.french.glass_center;
    const Real R = stirring.french.glass_radius;
    const Real hh = stirring.french.glass_half_height;
    b.lower_ = b.lower_.cwiseMin(gc - Vecd(R, R, hh));
    b.upper_ = b.upper_.cwiseMax(gc + Vecd(R, R, hh));
    b.lower_ -= Vecd(wall + pad, wall + pad, wall + pad);
    b.upper_ += Vecd(wall + pad, wall + pad, pad);
    return b;
}

/** Advection reference speed from parametric melt radius (not French mm glass STL). */
inline Real clientStirringReferenceSpeed(const OphelieFrenchStirringCaseParams &stirring)
{
    return Real(1.2) * stirring.rotation_speed_rad_s * stirring.french.glass_radius;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLIENT_STIRRING_GEOMETRY_H
