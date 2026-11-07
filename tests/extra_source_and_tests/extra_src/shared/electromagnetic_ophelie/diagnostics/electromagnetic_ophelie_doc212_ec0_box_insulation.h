/**
 * @file electromagnetic_ophelie_doc212_ec0_box_insulation.h
 * @brief Doc212 E-C0: thick box, constant A0, insulating Neumann → J≈0 exact.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC212_EC0_BOX_INSULATION_H
#define ELECTROMAGNETIC_OPHELIE_DOC212_EC0_BOX_INSULATION_H

#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_doc201_a1_phase_rotate.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "electromagnetic_ophelie_register_fields.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc212Ec0BoxCli
{
    bool enabled = false;
    Real dp = 0.02;
    Vecd half = Vecd(0.16, 0.12, 0.10);
    Vecd a0_real = Vecd(1.0e-4, 0.0, 0.0);
    Real frequency_hz = 50.0;
    Real sigma = 16.0;
};

struct OphelieDoc212Ec0BoxAudit
{
    size_t n_total = 0;
    size_t n_boundary = 0;
    size_t n_interior = 0;
    Real j_rms_over_jchar_all = 0.0;
    Real j_rms_over_jchar_int = 0.0;
    Real e_rms_over_wA_all = 0.0;
    Real e_rms_over_wA_int = 0.0;
    bool passed = false;
};

class OphelieDoc212Ec0BoxShape : public ComplexShape
{
  public:
    explicit OphelieDoc212Ec0BoxShape(const Vecd &center, const Vecd &half) : ComplexShape("Doc212Ec0Box")
    {
        add<GeometricShapeBox>(Transform(center), half);
    }
};

inline OphelieDoc212Ec0BoxAudit hostOphelieDoc212Ec0Audit(BaseParticles &particles,
                                                          const OphelieGlassFieldNames &names, Real omega,
                                                          Real sigma, const Vecd &a0_real, const Vecd &a0_imag)
{
    OphelieDoc212Ec0BoxAudit a;
    const size_t n = particles.TotalRealParticles();
    a.n_total = n;
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.phi_boundary_mask);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *er = particles.getVariableDataByName<Vecd>(names.e_edge_recon_real);
    const Vecd *ei = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Real *mask = particles.getVariableDataByName<Real>(names.phi_boundary_mask);

    const Real a0_abs = std::sqrt(a0_real.squaredNorm() + a0_imag.squaredNorm());
    const Real jchar = sigma * omega * a0_abs + TinyReal;
    const Real escale = omega * a0_abs + TinyReal;

    double sum_v = 0, sum_j2 = 0, sum_e2 = 0;
    double sum_vi = 0, sum_j2i = 0, sum_e2i = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real v = vol[i];
        const Real j2 = jr[i].squaredNorm() + ji[i].squaredNorm();
        const Real e2 = er[i].squaredNorm() + ei[i].squaredNorm();
        sum_v += v;
        sum_j2 += v * j2;
        sum_e2 += v * e2;
        if (mask[i] > Real(0.5))
        {
            ++a.n_boundary;
        }
        else
        {
            ++a.n_interior;
            sum_vi += v;
            sum_j2i += v * j2;
            sum_e2i += v * e2;
        }
    }
    a.j_rms_over_jchar_all = std::sqrt(sum_j2 / (sum_v + TinyReal)) / jchar;
    a.j_rms_over_jchar_int = std::sqrt(sum_j2i / (sum_vi + TinyReal)) / jchar;
    a.e_rms_over_wA_all = std::sqrt(sum_e2 / (sum_v + TinyReal)) / escale;
    a.e_rms_over_wA_int = std::sqrt(sum_e2i / (sum_vi + TinyReal)) / escale;

    const bool mask_ok = a.n_interior >= 100 && a.n_boundary > 0 && a.n_boundary < a.n_total;
    // Plan: J RMS/Jchar ≤1%; use |E|/(ω|A|) ≤3% as gradient-proxy (E=-∇φ-iωA → 0).
    const bool j_ok = a.j_rms_over_jchar_all <= Real(0.01) && a.j_rms_over_jchar_int <= Real(0.01);
    const bool e_ok = a.e_rms_over_wA_all <= Real(0.03) && a.e_rms_over_wA_int <= Real(0.03);
    a.passed = mask_ok && j_ok && e_ok;
    return a;
}

template <class ExecutionPolicy>
inline int runOphelieDoc212Ec0BoxInsulation(const OphelieDoc212Ec0BoxCli &cli)
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
    params.phi_boundary_mode_ = OpheliePhiBoundaryMode::OneSidedNeumann;
    params.phi_rhs_project_zero_mean_ = true;
    params.phi_gauge_penalty_ = Real(1.0e-8);
    params.phi_boundary_lhs_grad_neumann_ = false;
    params.phi_boundary_grad_neumann_projection_ = false;
    params.phi_boundary_normal_source_ = OpheliePhiBoundaryNormalSource::AnalyticBox;
    params.doc212_shared_complex_edge_rhs_ = true;

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

    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = half;
    const OpheliePhiBoundaryGeometryContext *geom_ptr =
        prepareOpheliePhiBoundaryGeomForSolve(particles, names, params, geom, cli.dp);

    std::cout << "[ophelie][doc212][ec0] box half=(" << half[0] << "," << half[1] << "," << half[2] << ") dp="
              << cli.dp << " f=" << cli.frequency_hz << " sigma=" << cli.sigma << " |A0|=" << cli.a0_real.norm()
              << " n=" << n << std::endl;

    auto run_once = [&](const char *tag) -> OphelieDoc212Ec0BoxAudit {
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                             geom_ptr, cli.dp);
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Vecd>(particles, names.a_coil_real);
        syncVariableToHost<Vecd>(particles, names.a_coil_imag);
        const Vecd ar0 = particles.getVariableDataByName<Vecd>(names.a_coil_real)[0];
        const Vecd ai0 = particles.getVariableDataByName<Vecd>(names.a_coil_imag)[0];
        OphelieDoc212Ec0BoxAudit aud =
            hostOphelieDoc212Ec0Audit(particles, names, params.omega(), params.sigma_glass_, ar0, ai0);
        std::cout << "[ophelie][doc212][ec0] " << tag << " n_bnd=" << aud.n_boundary << " n_int=" << aud.n_interior
                  << " Jrms/Jchar all=" << aud.j_rms_over_jchar_all << " int=" << aud.j_rms_over_jchar_int
                  << " Erms/(w|A|) all=" << aud.e_rms_over_wA_all << " int=" << aud.e_rms_over_wA_int
                  << " passed=" << (aud.passed ? 1 : 0) << std::endl;
        return aud;
    };

    const OphelieDoc212Ec0BoxAudit a0 = run_once("A0");
    hostOphelieMultiplyCoilFieldsByI(particles, names, n);
    combine.exec();
    // Refresh gn for iA (boundary fields depend on A).
    geom_ptr = prepareOpheliePhiBoundaryGeomForSolve(particles, names, params, geom, cli.dp);
    const OphelieDoc212Ec0BoxAudit ai = run_once("iA");

    const bool passed = a0.passed && ai.passed;
    std::cout << "test_3d_ophelie_analytic_cylinder_skin doc212_ec0 n=" << n << " dp=" << cli.dp
              << " Jrms_all=" << a0.j_rms_over_jchar_all << " Erms_all=" << a0.e_rms_over_wA_all
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
