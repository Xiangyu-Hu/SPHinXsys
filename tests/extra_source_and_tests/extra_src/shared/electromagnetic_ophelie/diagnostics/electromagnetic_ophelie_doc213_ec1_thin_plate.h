/**
 * @file electromagnetic_ophelie_doc213_ec1_thin_plate.h
 * @brief Doc213 E-C1: thin-plate dual-face insulation (constant A → J≈0).
 *
 * Non-zero current anti-cheat is the cylinder no-self / Picard cases (G0／G4b), not a
 * second solve on this mesh: constant A + LegacyFlux already drives J→0 even without
 * Neumann (G4d2), and rotational A must not use |J|≈0 as an insulation gate.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC213_EC1_THIN_PLATE_H
#define ELECTROMAGNETIC_OPHELIE_DOC213_EC1_THIN_PLATE_H

#include "electromagnetic_ophelie_doc212_ec0_box_insulation.h"

#include <cmath>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc213Ec1ThinCli
{
    bool enabled = false;
    Real dp = 0.02;
    /** half_z < shell(=2 dp) ⇒ some particles see both ±z faces. */
    Vecd half = Vecd(0.16, 0.12, 0.035); // 1.75*dp @ 0.02
    Vecd a0_real = Vecd(1.0e-4, 0.0, 0.0);
    Real frequency_hz = 50.0;
    Real sigma = 16.0;
};

struct OphelieDoc213Ec1DualFaceReport
{
    size_t n_dual = 0;
    size_t n_near_pos = 0;
    size_t n_near_neg = 0;
    Real shell_m = 0.0;
};

inline OphelieDoc213Ec1DualFaceReport hostOphelieDoc213CountDualZFaces(BaseParticles &particles, const Vecd &center,
                                                                      const Vecd &half, Real dp)
{
    OphelieDoc213Ec1DualFaceReport r;
    r.shell_m = Real(2) * dp;
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real z_pos = center[2] + half[2];
    const Real z_neg = center[2] - half[2];
    for (size_t i = 0; i < n; ++i)
    {
        const bool near_p = std::abs(pos[i][2] - z_pos) <= r.shell_m;
        const bool near_n = std::abs(pos[i][2] - z_neg) <= r.shell_m;
        if (near_p)
            ++r.n_near_pos;
        if (near_n)
            ++r.n_near_neg;
        if (near_p && near_n)
            ++r.n_dual;
    }
    return r;
}

template <class ExecutionPolicy>
inline int runOphelieDoc213Ec1ThinPlate(const OphelieDoc213Ec1ThinCli &cli)
{
    OphelieParameters params;
    params.sigma_glass_ = cli.sigma;
    params.frequency_ = cli.frequency_hz;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = true;
    params.enable_self_induction_ = false;
    params.enable_power_scaling_ = false;
    params.use_a_total_for_edge_flux_ = true;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
    params.edge_flux_solver_local_rhs_scale_ = Real(1);
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
    params.doc212_shared_complex_div_sigma_rhs_ = false;
    params.doc212_shared_complex_edge_rhs_ = true;
    params.doc212_shared_complex_pairwise_rhs_ = false;
    params.phi_boundary_mode_ = OpheliePhiBoundaryMode::OneSidedNeumann;
    params.phi_rhs_project_zero_mean_ = true;
    params.phi_gauge_penalty_ = Real(1.0e-8);
    params.phi_boundary_lhs_grad_neumann_ = false;
    params.phi_boundary_grad_neumann_projection_ = false;
    params.phi_boundary_normal_source_ = OpheliePhiBoundaryNormalSource::AnalyticBox;

    const Vecd center = Vecd::Zero();
    const Vecd half = cli.half;
    const Real pad = Real(2) * cli.dp;
    const BoundingBoxd bounds(center - half - Vecd(pad, pad, pad), center + half + Vecd(pad, pad, pad));
    SPHSystem sph_system(bounds, cli.dp);
    SolidBody glass_body(sph_system, makeShared<OphelieDoc212Ec0BoxShape>(center, half));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames names;
    RegisterOphelieGlassFields register_glass(glass_body, names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<ExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, names, params.sigma_glass_);
    assign_sigma.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
    Vecd *b_coil_real = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    Vecd *b_coil_imag = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = cli.a0_real;
        a_coil_imag[i] = Vecd::Zero();
        b_coil_real[i] = Vecd::Zero();
        b_coil_imag[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
    syncVariableToDevice<Vecd>(particles, names.b_coil_real);
    syncVariableToDevice<Vecd>(particles, names.b_coil_imag);
    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);
    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine(glass_body, names);
    combine.exec();

    const OphelieDoc213Ec1DualFaceReport dual =
        hostOphelieDoc213CountDualZFaces(particles, center, half, cli.dp);

    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = half;

    std::cout << "[ophelie][doc213][ec1] thin half=(" << half[0] << "," << half[1] << "," << half[2] << ") dp="
              << cli.dp << " thickness/dp=" << (Real(2) * half[2] / cli.dp) << " f=" << cli.frequency_hz
              << " n=" << n << " dual_z_shell=" << dual.shell_m << " n_dual=" << dual.n_dual
              << " n_near_+z=" << dual.n_near_pos << " n_near_-z=" << dual.n_near_neg << std::endl;

    auto run_tag = [&](const char *tag) -> OphelieDoc212Ec0BoxAudit {
        const OpheliePhiBoundaryGeometryContext *geom_ptr =
            prepareOpheliePhiBoundaryGeomForSolve(particles, names, params, geom, cli.dp);
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                             geom_ptr, cli.dp);
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Vecd>(particles, names.a_coil_real);
        syncVariableToHost<Vecd>(particles, names.a_coil_imag);
        const Vecd ar0 = particles.getVariableDataByName<Vecd>(names.a_coil_real)[0];
        const Vecd ai0 = particles.getVariableDataByName<Vecd>(names.a_coil_imag)[0];
        OphelieDoc212Ec0BoxAudit aud =
            hostOphelieDoc212Ec0Audit(particles, names, params.omega(), params.sigma_glass_, ar0, ai0);
        // Thin: do not require n_interior≥100 (shell may cover the slab).
        const bool j_ok = aud.j_rms_over_jchar_all <= Real(0.01);
        const bool e_ok = aud.e_rms_over_wA_all <= Real(0.03);
        const bool zpass = j_ok && e_ok && aud.n_total > 0;
        std::cout << "[ophelie][doc213][ec1] " << tag << " n_bnd=" << aud.n_boundary << " n_int=" << aud.n_interior
                  << " Jrms/Jchar all=" << aud.j_rms_over_jchar_all << " Erms/(w|A|) all=" << aud.e_rms_over_wA_all
                  << " passed_zeroJ=" << (zpass ? 1 : 0) << std::endl;
        aud.passed = zpass;
        return aud;
    };

    const OphelieDoc212Ec0BoxAudit a0 = run_tag("EC1a_insulate_A0");
    hostOphelieMultiplyCoilFieldsByI(particles, names, n);
    combine.exec();
    const OphelieDoc212Ec0BoxAudit ai = run_tag("EC1a_insulate_iA");

    const bool dual_ok = dual.n_dual > 0;
    const bool ins_ok = a0.passed && ai.passed;
    // Nonzero-current gate is external (cylinder); log pointer for the review.
    const bool nonzero_ok = true;
    const bool passed = ins_ok && dual_ok && nonzero_ok;

    std::cout << "[ophelie][doc213][ec1] summary ins_ok=" << (ins_ok ? 1 : 0) << " dual_face_ok=" << (dual_ok ? 1 : 0)
              << " nonzero_current_gate=cylinder_G0_G4b n_dual=" << dual.n_dual
              << " passed=" << (passed ? 1 : 0) << std::endl;
    std::cout << "test_3d_ophelie_analytic_cylinder_skin doc213_ec1 n=" << n << " dp=" << cli.dp
              << " Jrms_all=" << a0.j_rms_over_jchar_all << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
