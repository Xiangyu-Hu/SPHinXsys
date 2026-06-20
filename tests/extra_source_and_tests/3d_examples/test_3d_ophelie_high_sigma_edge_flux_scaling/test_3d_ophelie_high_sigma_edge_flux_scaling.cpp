/**
 * @file test_3d_ophelie_high_sigma_edge_flux_scaling.cpp
 * @brief High-sigma edge-flux benchmarks: uniform E, constant-A gauge cancellation, rotational A.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_aind_lenz_audit.h"
#include "electromagnetic_ophelie_p6b_box_boundary.h"
#include "electromagnetic_ophelie_phi_boundary_diagnostics.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

enum class HighSigmaCaseKind
{
    UniformEReal,
    /** Spatially constant A (pure gauge): edge-only reconstruction diagnostic. */
    ConstantAGaugeCancellationEdgeOnly,
    /** Full φ solve cancels constant A; J≈0 is expected (gauge cancellation). */
    ConstantAGaugeCancellationPhiSolve,
    ConstantAGaugeCancellationPhiSolveSolverLocal,
    /** A=0.5*B0×r, curl A=B0 z — non-conservative drive; edge-only diagnostic. */
    RotationalAUniformBEdgeOnly,
    /** Full φ solve on rotational A; J must NOT gauge-cancel. */
    RotationalAUniformBPhiSolve
};

struct HighSigmaBenchmarkRecord
{
    std::string case_name;
    Real sigma = 0.0;
    Real frequency_hz = 0.0;
    std::string normalization_mode;
    Real input_scale_applied = 1.0;
    Real input_scale_measured = 1.0;
    size_t n_particles = 0;
    Real e_recon_over_exact = 0.0;
    Real j_recon_over_sigma_e = 0.0;
    Real e_edge_recon_over_exact = 0.0;
    Real j_edge_recon_over_sigma_e = 0.0;
    Real max_j_imag = 0.0;
    Real max_j_edge_imag = 0.0;
    Real p_recon_over_exact = 0.0;
    Real e_edge_em_mismatch = 0.0;
    Real p_graph_over_p_recon = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real phi_eq_res_imag = 0.0;
    Real restore_invariance_error_j = 0.0;
    bool gauge_cancellation_passed = false;
    Real phi_linear_fit_error = 0.0;
    Real edge_drop_after_phi_l2 = 0.0;
    Real edge_drop_after_phi_linf = 0.0;
    Real j_after_phi_norm = 0.0;
    Real j_after_phi_over_edge_only = 0.0;
    Real j_edge_only_norm = 0.0;
    Real radial_e_leak = 0.0;
    Real normal_j_boundary = 0.0;
    Real q_antisym_rel_l2 = 0.0;
    Real edge_res_red = 0.0;
    Real corr_b_ind_z_b_drive_z = 0.0;
    Real b_ind_vol_norm_lenz = 0.0;
    bool lenz_opposes_uniform_b_drive = false;
    bool finite_fields = false;
    bool row_passed = false;
    std::string boundary_mode = "none";
    Real boundary_jn_rel = 0.0;
};

inline Real hostTotalVolume(BaseParticles &particles, size_t n)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real total = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        total += vol[i];
    }
    return total;
}

inline Real hostVecdVolWeightedRelErrorToConstant(BaseParticles &particles, const std::string &field_name,
                                                 const Vecd &exact, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real err_sq = 0.0;
    Real ref_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd delta = values[i] - exact;
        err_sq += vol[i] * delta.squaredNorm();
        ref_sq += vol[i] * exact.squaredNorm();
    }
    return std::sqrt(err_sq) / (std::sqrt(ref_sq) + TinyReal);
}

inline Real hostVecdVolWeightedRelErrorToField(BaseParticles &particles, const std::string &field_name,
                                               const StdVec<Vecd> &exact, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real err_sq = 0.0;
    Real ref_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd delta = values[i] - exact[i];
        err_sq += vol[i] * delta.squaredNorm();
        ref_sq += vol[i] * exact[i].squaredNorm();
    }
    return std::sqrt(err_sq) / (std::sqrt(ref_sq) + TinyReal);
}

inline Real hostScalarVolWeightedRelError(BaseParticles &particles, const std::string &field_name,
                                        const StdVec<Real> &exact, size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real err_sq = 0.0;
    Real ref_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real delta = values[i] - exact[i];
        err_sq += vol[i] * delta * delta;
        ref_sq += vol[i] * exact[i] * exact[i];
    }
    return std::sqrt(err_sq) / (std::sqrt(ref_sq) + TinyReal);
}

inline bool hostFieldsFinite(BaseParticles &particles, const OphelieGlassFieldNames &names, size_t n)
{
    syncVariableToHost<Vecd>(particles, names.e_imag);
    syncVariableToHost<Vecd>(particles, names.j_imag);
    const Vecd *e = particles.getVariableDataByName<Vecd>(names.e_imag);
    const Vecd *j = particles.getVariableDataByName<Vecd>(names.j_imag);
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(e[i].sum()) || !std::isfinite(j[i].sum()))
        {
            return false;
        }
    }
    return true;
}

inline Vecd rotationalAUniformB(const Vecd &pos, Real b0)
{
    return Vecd(-Real(0.5) * b0 * pos[1], Real(0.5) * b0 * pos[0], 0.0);
}

inline Real hostBoundaryNormalJNorm(BaseParticles &particles, const std::string &j_field, const Vecd &box_lower,
                                    const Vecd &box_upper, Real skin_h, size_t n)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, j_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j = particles.getVariableDataByName<Vecd>(j_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real jn_sq = 0.0;
    Real vol_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx_lo = pos[i][0] - box_lower[0];
        const Real dx_hi = box_upper[0] - pos[i][0];
        const Real dy_lo = pos[i][1] - box_lower[1];
        const Real dy_hi = box_upper[1] - pos[i][1];
        const Real dz_lo = pos[i][2] - box_lower[2];
        const Real dz_hi = box_upper[2] - pos[i][2];
        const Real min_dist = std::min({dx_lo, dx_hi, dy_lo, dy_hi, dz_lo, dz_hi});
        if (min_dist > skin_h)
        {
            continue;
        }
        Vecd n_hat = Vecd::Zero();
        if (dx_lo == min_dist)
        {
            n_hat[0] = -1.0;
        }
        else if (dx_hi == min_dist)
        {
            n_hat[0] = 1.0;
        }
        else if (dy_lo == min_dist)
        {
            n_hat[1] = -1.0;
        }
        else if (dy_hi == min_dist)
        {
            n_hat[1] = 1.0;
        }
        else if (dz_lo == min_dist)
        {
            n_hat[2] = -1.0;
        }
        else
        {
            n_hat[2] = 1.0;
        }
        const Real jn = j[i].dot(n_hat);
        jn_sq += vol[i] * jn * jn;
        vol_sum += vol[i];
    }
    (void)vol_sum;
    return std::sqrt(jn_sq);
}

inline Real hostRadialELeak(BaseParticles &particles, const std::string &e_field, const Vecd &center, size_t n)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, e_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *e = particles.getVariableDataByName<Vecd>(e_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real leak_sq = 0.0;
    Real e_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd r = pos[i] - center;
        const Real r_norm = r.head<2>().norm() + TinyReal;
        const Vecd r_hat(r[0] / r_norm, r[1] / r_norm, 0.0);
        const Real er = e[i].dot(r_hat);
        leak_sq += vol[i] * er * er;
        e_sq += vol[i] * e[i].squaredNorm();
    }
    return std::sqrt(leak_sq) / (std::sqrt(e_sq) + TinyReal);
}

inline bool evaluateHighSigmaRowPass(const HighSigmaBenchmarkRecord &r)
{
    if (!r.finite_fields)
    {
        return false;
    }
    if (r.case_name == "uniform_E_real")
    {
        return r.e_recon_over_exact < Real(1e-4) && r.j_recon_over_sigma_e < Real(1e-4) && r.p_recon_over_exact > Real(0.5) &&
               r.p_recon_over_exact < Real(2.0);
    }
    if (r.case_name == "constant_A_edge_only")
    {
        return r.e_recon_over_exact < Real(1e-4) && r.j_recon_over_sigma_e < Real(1e-4);
    }
    if (r.case_name == "constant_A_gauge_cancellation" ||
        r.case_name == "constant_A_gauge_cancellation_solver_local")
    {
        return r.gauge_cancellation_passed;
    }
    if (r.case_name == "rotational_A_uniform_B_edge_only")
    {
        return r.e_recon_over_exact < Real(0.15) && r.j_recon_over_sigma_e < Real(0.15);
    }
    if (r.case_name == "rotational_A_uniform_B_phi_solve")
    {
        return r.e_recon_over_exact < Real(0.5) && r.j_recon_over_sigma_e < Real(0.5) &&
               r.j_after_phi_over_edge_only > Real(0.05);
    }
    return false;
}

inline HighSigmaBenchmarkRecord runHighSigmaBenchmarkCase(HighSigmaCaseKind case_kind, Real sigma, Real frequency_hz,
                                                          Real dp, const Vecd &center, const Vecd &halfsize,
                                                          OphelieParameters &params, Real b0_tesla = Real(0.01))
{
    HighSigmaBenchmarkRecord record;
    record.sigma = sigma;
    record.frequency_hz = frequency_hz;

    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));
    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    params.sigma_glass_ = sigma;
    params.frequency_ = frequency_hz;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = true;

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, sigma);
    assign_sigma.exec();

    if (ophelieParamsNeedBoxEdgeFluxNormals(params))
    {
        initBoxEdgeFluxBoundaryNormalsFromShape(sph_system, glass_body);
    }

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    record.n_particles = n;

    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    Real *phi_real = particles.getVariableDataByName<Real>(glass_names.phi_real);

    const Real omega = params.omega();
    const Vecd box_lower = center - halfsize;
    const Vecd box_upper = center + halfsize;

    const bool is_constant_a =
        case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationEdgeOnly ||
        case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolve ||
        case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolveSolverLocal;
    const bool is_rotational_a = case_kind == HighSigmaCaseKind::RotationalAUniformBEdgeOnly ||
                                 case_kind == HighSigmaCaseKind::RotationalAUniformBPhiSolve;
    const bool is_phi_solve = case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolve ||
                              case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolveSolverLocal ||
                              case_kind == HighSigmaCaseKind::RotationalAUniformBPhiSolve;

    const OphelieEdgeFluxNormalizationMode normalization_mode_saved = params.edge_flux_normalization_mode_;
    if (case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolveSolverLocal ||
        case_kind == HighSigmaCaseKind::RotationalAUniformBPhiSolve)
    {
        params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::SolverLocal;
    }
    record.normalization_mode = ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_);
    record.boundary_mode = p6bBoundarySweepLabel(params);

    StdVec<Vecd> e_exact_field(n, Vecd::Zero());
    StdVec<Vecd> j_exact_field(n, Vecd::Zero());
    StdVec<Real> phi_exact_imag(n, 0.0);
    Vecd e_exact_const = Vecd::Zero();
    Vecd j_exact_const = Vecd::Zero();
    Real p_exact = 0.0;

    std::string e_primary_field = glass_names.e_imag;
    std::string j_primary_field = glass_names.j_imag;
    std::string e_edge_field = glass_names.e_edge_recon_imag;
    std::string j_edge_field = glass_names.j_edge_recon_imag;

    if (case_kind == HighSigmaCaseKind::UniformEReal)
    {
        record.case_name = "uniform_E_real";
        const Vecd e0(100.0, 0.0, 0.0);
        e_exact_const = e0;
        j_exact_const = sigma * e0;
        e_primary_field = glass_names.e_real;
        j_primary_field = glass_names.j_real;
        e_edge_field = glass_names.e_edge_recon_real;
        j_edge_field = glass_names.j_edge_recon_real;
        for (size_t i = 0; i < n; ++i)
        {
            a_coil_real[i] = Vecd::Zero();
            a_coil_imag[i] = Vecd::Zero();
            phi_imag[i] = Real(0);
            phi_real[i] = -e0.dot(pos[i]);
        }
    }
    else if (is_constant_a)
    {
        if (case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationEdgeOnly)
        {
            record.case_name = "constant_A_edge_only";
        }
        else if (case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolveSolverLocal)
        {
            record.case_name = "constant_A_gauge_cancellation_solver_local";
        }
        else
        {
            record.case_name = "constant_A_gauge_cancellation";
        }
        const Vecd a0(0.001, 0.0, 0.0);
        e_exact_const = -omega * a0;
        j_exact_const = sigma * e_exact_const;
        for (size_t i = 0; i < n; ++i)
        {
            a_coil_real[i] = a0;
            a_coil_imag[i] = Vecd::Zero();
            phi_imag[i] = Real(0);
            phi_real[i] = Real(0);
            phi_exact_imag[i] = -omega * a0.dot(pos[i]);
            e_exact_field[i] = e_exact_const;
            j_exact_field[i] = j_exact_const;
        }
        record.j_edge_only_norm = j_exact_const.norm() * std::sqrt(hostTotalVolume(particles, n));
    }
    else if (is_rotational_a)
    {
        record.case_name = case_kind == HighSigmaCaseKind::RotationalAUniformBEdgeOnly ? "rotational_A_uniform_B_edge_only"
                                                                                       : "rotational_A_uniform_B_phi_solve";
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        Real j_edge_sq = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            a_coil_real[i] = rotationalAUniformB(pos[i], b0_tesla);
            a_coil_imag[i] = Vecd::Zero();
            phi_imag[i] = Real(0);
            phi_real[i] = Real(0);
            e_exact_field[i] = -omega * a_coil_real[i];
            j_exact_field[i] = sigma * e_exact_field[i];
            j_edge_sq += vol[i] * j_exact_field[i].squaredNorm();
        }
        record.j_edge_only_norm = std::sqrt(j_edge_sq);
    }

    const Real total_volume = hostTotalVolume(particles, n);
    if (case_kind == HighSigmaCaseKind::UniformEReal || is_constant_a)
    {
        p_exact = Real(0.5) * sigma * e_exact_const.squaredNorm() * total_volume;
    }
    else if (is_rotational_a)
    {
        Real p_sum = 0.0;
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        for (size_t i = 0; i < n; ++i)
        {
            p_sum += Real(0.5) * sigma * e_exact_field[i].squaredNorm() * vol[i];
        }
        p_exact = p_sum;
    }

    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);

    Real input_scale = 1.0;
    if (is_phi_solve && ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_))
    {
        Real rhs_l2_pre = 0.0;
        input_scale = applyOphelieEdgeFluxInputNormalization<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, n, params.edge_flux_safe_rhs_l2_,
            params.edge_flux_safe_rhs_max_abs_, &rhs_l2_pre, &record.input_scale_measured);
    }
    else
    {
        record.input_scale_measured = 1.0;
    }
    record.input_scale_applied = input_scale;

    Real p_recon = 0.0;
    if (case_kind == HighSigmaCaseKind::UniformEReal)
    {
        const OphelieEdgeFluxPowerMetrics power_metrics =
            execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
        syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, glass_names, params);
        p_recon = power_metrics.p_total_recon;
        record.phi_eq_res_imag = 0.0;
    }
    else if (is_phi_solve)
    {
        OphelieComplexEdgeFluxSolveReport solve_report;
        const OpheliePhiBoundaryGeometryContext box_geom = makeAnalyticBoxPhiBoundaryGeometry(center, halfsize);
        const OpheliePhiBoundaryGeometryContext *boundary_geom =
            params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann ? &box_geom : nullptr;
        p_recon = execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, boundary_geom, dp, &solve_report);
        record.phi_eq_res_imag = solve_report.phi_eq_res_vol_imag;
        record.edge_res_red = solve_report.phi_eq_res_vol_imag;
    }
    else
    {
        const OphelieEdgeFluxPowerMetrics power_metrics =
            execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
        syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, glass_names, params);
        p_recon = power_metrics.p_total_recon;
        record.phi_eq_res_imag = 0.0;
    }

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(*glass_inner,
                                                                                                            glass_names);
    compute_grad_phi.exec();

    record.e_edge_em_mismatch =
        (is_constant_a || is_rotational_a)
            ? hostOphelieImagEdgeEmfMismatchVolRatio(particles, glass_names, params, n)
            : Real(0);

    if (case_kind == HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolve &&
        ophelieEdgeFluxUsesFieldScaleRestore(params.edge_flux_normalization_mode_) &&
        std::abs(input_scale - 1.0) > TinyReal)
    {
        restoreOphelieEdgeFluxInputNormalization<MainExecutionPolicy>(glass_body, glass_names, params, input_scale, n);
    }

    const OphelieEdgeFluxEdgeDropMetrics edge_drop_metrics =
        evaluateOphelieEdgeFluxImagEdgeDropMetrics<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    record.edge_drop_after_phi_l2 = edge_drop_metrics.edge_drop_l2;
    record.edge_drop_after_phi_linf = edge_drop_metrics.edge_drop_linf;

    const OphelieEdgeFluxQAntisymMetrics q_metrics =
        evaluateOphelieEdgeFluxQAntisymmetry<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    record.q_antisym_rel_l2 = q_metrics.q_antisym_rel_l2;

    if (is_rotational_a && is_phi_solve)
    {
        const OphelieAIndLenzAuditRecord lenz = execOphelieGlassSelfInducedBiotAndLenzAudit<MainExecutionPolicy>(
            glass_body, glass_names, params, glass_names.j_real, j_primary_field, b0_tesla);
        record.corr_b_ind_z_b_drive_z = lenz.corr_b_ind_z_b_drive_z;
        record.b_ind_vol_norm_lenz = lenz.b_ind_vol_norm;
        record.lenz_opposes_uniform_b_drive = lenz.lenz_opposes_uniform_b_drive;
        printOphelieAIndLenzAudit(lenz, record.case_name);
    }

    record.j_imag_vol_norm = hostVecdVolWeightedNorm(particles, j_primary_field, n);
    record.j_after_phi_norm = record.j_imag_vol_norm;
    record.j_after_phi_over_edge_only = record.j_imag_vol_norm / (record.j_edge_only_norm + TinyReal);

    if (is_constant_a)
    {
        record.phi_linear_fit_error = hostScalarVolWeightedRelError(particles, glass_names.phi_imag, phi_exact_imag, n);
    }

    if (is_rotational_a)
    {
        record.radial_e_leak = hostRadialELeak(particles, e_primary_field, center, n);
        record.normal_j_boundary =
            hostBoundaryNormalJNorm(particles, j_primary_field, box_lower, box_upper, dp, n);
        record.e_recon_over_exact = hostVecdVolWeightedRelErrorToField(particles, e_primary_field, e_exact_field, n);
        record.j_recon_over_sigma_e = hostVecdVolWeightedRelErrorToField(particles, j_primary_field, j_exact_field, n);
        record.e_edge_recon_over_exact = hostVecdVolWeightedRelErrorToField(particles, e_edge_field, e_exact_field, n);
        record.j_edge_recon_over_sigma_e = hostVecdVolWeightedRelErrorToField(particles, j_edge_field, j_exact_field, n);
    }
    else if (case_kind != HighSigmaCaseKind::UniformEReal)
    {
        record.e_recon_over_exact = hostVecdVolWeightedRelErrorToConstant(particles, e_primary_field, e_exact_const, n);
        record.j_recon_over_sigma_e = hostVecdVolWeightedRelErrorToConstant(particles, j_primary_field, j_exact_const, n);
        record.e_edge_recon_over_exact = hostVecdVolWeightedRelErrorToConstant(particles, e_edge_field, e_exact_const, n);
        record.j_edge_recon_over_sigma_e = hostVecdVolWeightedRelErrorToConstant(particles, j_edge_field, j_exact_const, n);
    }
    else
    {
        record.e_recon_over_exact = hostVecdVolWeightedRelErrorToConstant(particles, e_primary_field, e_exact_const, n);
        record.j_recon_over_sigma_e = hostVecdVolWeightedRelErrorToConstant(particles, j_primary_field, j_exact_const, n);
        record.e_edge_recon_over_exact = hostVecdVolWeightedRelErrorToConstant(particles, e_edge_field, e_exact_const, n);
        record.j_edge_recon_over_sigma_e = hostVecdVolWeightedRelErrorToConstant(particles, j_edge_field, j_exact_const, n);
    }

    record.max_j_imag = hostVecdFieldMax(particles, j_primary_field, n);
    record.max_j_edge_imag = hostVecdFieldMax(particles, j_edge_field, n);
    record.p_recon_over_exact = p_recon / (p_exact + TinyReal);

    const OphelieEdgeFluxPowerAuditDetail power_detail = hostOphelieEdgeFluxPowerAuditDetail(particles, glass_names, n, params);
    record.p_graph_over_p_recon = power_detail.p_graph_over_recon;
    record.finite_fields = hostFieldsFinite(particles, glass_names, n);
    if (is_constant_a && is_phi_solve)
    {
        record.gauge_cancellation_passed =
            record.finite_fields && record.j_after_phi_over_edge_only < Real(0.01) &&
            record.edge_drop_after_phi_l2 < Real(1e-2) * (std::abs(omega * Real(0.001)) + Real(1.0)) &&
            record.phi_linear_fit_error < Real(0.05);
    }
    record.row_passed = evaluateHighSigmaRowPass(record);

    if (is_phi_solve)
    {
        const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
        const OpheliePhiBoundaryJnMetrics jn_metrics =
            computeBoxBoundaryJnMetrics(particles, glass_names, n, params, center, halfsize, boundary_width);
        record.boundary_jn_rel = jn_metrics.jn_boundary_rel;
    }

    std::cout << "[p1a] case=" << record.case_name << " boundary=" << record.boundary_mode << " sigma=" << sigma << " f_hz=" << frequency_hz
              << " norm=" << record.normalization_mode << " scale=" << record.input_scale_applied << " n=" << n
              << " E_rel=" << record.e_recon_over_exact << " J_rel=" << record.j_recon_over_sigma_e
              << " gauge_pass=" << (record.gauge_cancellation_passed ? 1 : 0)
              << " J_phi/J_edge=" << record.j_after_phi_over_edge_only
              << " edge_drop_l2=" << record.edge_drop_after_phi_l2
              << " phi_fit_err=" << record.phi_linear_fit_error << " radial_E_leak=" << record.radial_e_leak
              << " normal_J_bnd=" << record.normal_j_boundary << " q_antisym=" << record.q_antisym_rel_l2
              << " corr_Bind_z_Bdrive=" << record.corr_b_ind_z_b_drive_z
              << " lenz_vs_drive=" << (record.lenz_opposes_uniform_b_drive ? 1 : 0)
              << " boundary_Jn_rel=" << record.boundary_jn_rel
              << " P_ratio=" << record.p_recon_over_exact << " pass=" << (record.row_passed ? 1 : 0) << std::endl;
    params.edge_flux_normalization_mode_ = normalization_mode_saved;
    return record;
}

inline void writeHighSigmaBenchmarkCsv(const std::string &path, const StdVec<HighSigmaBenchmarkRecord> &records)
{
    const fs::path parent = fs::path(path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(path);
    if (!out)
    {
        std::cerr << "error: could not write " << path << std::endl;
        return;
    }
    out << "case_name,sigma,frequency_hz,normalization_mode,input_scale_applied,input_scale_measured,n_particles,"
           "e_recon_over_exact,j_recon_over_sigma_e,e_edge_recon_over_exact,j_edge_recon_over_sigma_e,"
           "p_recon_over_exact,e_edge_em_mismatch,p_graph_over_p_recon,j_imag_vol_norm,max_j_primary,max_j_edge,"
           "phi_eq_res_imag,restore_invariance_error_j,gauge_cancellation_passed,phi_linear_fit_error,"
           "edge_drop_after_phi_l2,edge_drop_after_phi_linf,j_after_phi_norm,j_after_phi_over_edge_only,"
           "j_edge_only_norm,radial_e_leak,normal_j_boundary,boundary_mode,boundary_jn_rel,q_antisym_rel_l2,edge_res_red,corr_b_ind_z_b_drive_z,"
           "b_ind_vol_norm_lenz,lenz_opposes_uniform_b_drive,finite_fields,row_passed\n";
    out << std::setprecision(10);
    for (const HighSigmaBenchmarkRecord &r : records)
    {
        out << r.case_name << "," << r.sigma << "," << r.frequency_hz << "," << r.normalization_mode << ","
            << r.input_scale_applied << "," << r.input_scale_measured << "," << r.n_particles << ","
            << r.e_recon_over_exact << "," << r.j_recon_over_sigma_e << "," << r.e_edge_recon_over_exact << ","
            << r.j_edge_recon_over_sigma_e << "," << r.p_recon_over_exact << "," << r.e_edge_em_mismatch << ","
            << r.p_graph_over_p_recon << "," << r.j_imag_vol_norm << "," << r.max_j_imag << "," << r.max_j_edge_imag
            << "," << r.phi_eq_res_imag << "," << r.restore_invariance_error_j << ","
            << (r.gauge_cancellation_passed ? 1 : 0) << "," << r.phi_linear_fit_error << ","
            << r.edge_drop_after_phi_l2 << "," << r.edge_drop_after_phi_linf << "," << r.j_after_phi_norm << ","
            << r.j_after_phi_over_edge_only << "," << r.j_edge_only_norm << "," << r.radial_e_leak << ","
            << r.normal_j_boundary << "," << r.boundary_mode << "," << r.boundary_jn_rel << "," << r.q_antisym_rel_l2
            << "," << r.edge_res_red << "," << r.corr_b_ind_z_b_drive_z << "," << r.b_ind_vol_norm_lenz << ","
            << (r.lenz_opposes_uniform_b_drive ? 1 : 0) << ","
            << (r.finite_fields ? 1 : 0) << "," << (r.row_passed ? 1 : 0) << "\n";
    }
    std::cout << "[p1a] wrote CSV: " << path << " rows=" << records.size() << std::endl;
}

} // namespace

int main(int ac, char *av[])
{
    OphelieParameters params;
    OphelieTestCliOptions cli_options;
    const StdVec<std::string> filtered_arguments = filterOphelieTestCommandLine(ac, av, params, cli_options);
    (void)filtered_arguments;

    Real dp = 0.04;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            dp = static_cast<Real>(std::atof(av[i] + 5));
        }
    }

    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const StdVec<Real> sigma_values = {Real(16), Real(1.0e3), Real(1.0e5), Real(1.0e7), Real(3.526e7)};
    const StdVec<Real> frequency_values = {Real(50), Real(200)};

    std::string csv_path = "./output/high_sigma_edge_flux_scaling.csv";
    bool p6b_boundary_sweep = false;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--output-csv=", 13) == 0)
        {
            csv_path = std::string(av[i] + 13);
        }
        else if (std::strcmp(av[i], "--p6b-boundary-sweep=1") == 0 || std::strcmp(av[i], "--p6b-boundary-sweep") == 0)
        {
            p6b_boundary_sweep = true;
        }
    }

    if (p6b_boundary_sweep)
    {
        csv_path = "./output/p6b_box_boundary_sweep.csv";
        const Real sigma = Real(1.0e7);
        const Real frequency_hz = Real(50);
        struct BoundarySweepConfig
        {
            const char *label;
            OphelieEdgeReconBoundaryMode edge_mode;
            OpheliePhiBoundaryMode phi_mode;
        };
        const BoundarySweepConfig configs[] = {
            {"none", OphelieEdgeReconBoundaryMode::None, OpheliePhiBoundaryMode::None},
            {"no-flux-ghost-edge", OphelieEdgeReconBoundaryMode::NoFluxGhostEdge, OpheliePhiBoundaryMode::None},
            {"no-flux-full", OphelieEdgeReconBoundaryMode::NoFluxFull, OpheliePhiBoundaryMode::None},
            {"phi-neumann", OphelieEdgeReconBoundaryMode::None, OpheliePhiBoundaryMode::OneSidedNeumann},
        };
        const HighSigmaCaseKind cases[] = {HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolve,
                                           HighSigmaCaseKind::RotationalAUniformBPhiSolve};
        StdVec<HighSigmaBenchmarkRecord> records;
        std::cout << "[p6b] box/slab boundary sweep sigma=" << sigma << " f_hz=" << frequency_hz << " dp=" << dp
                  << std::endl;
        for (const BoundarySweepConfig &config : configs)
        {
            for (HighSigmaCaseKind case_kind : cases)
            {
                OphelieParameters sweep_params = params;
                sweep_params.edge_recon_boundary_mode_ = config.edge_mode;
                sweep_params.phi_boundary_mode_ = config.phi_mode;
                sweep_params.edge_recon_boundary_width_factor_ = Real(2);
                HighSigmaBenchmarkRecord record =
                    runHighSigmaBenchmarkCase(case_kind, sigma, frequency_hz, dp, center, halfsize, sweep_params);
                record.case_name += std::string("_") + config.label;
                records.push_back(record);
            }
        }
        writeHighSigmaBenchmarkCsv(csv_path, records);
        size_t n_pass = 0;
        for (const HighSigmaBenchmarkRecord &r : records)
        {
            if (r.row_passed)
            {
                ++n_pass;
            }
        }
        std::cout << "[p6b] sweep summary: pass=" << n_pass << " fail=" << (records.size() - n_pass)
                  << " csv=" << csv_path << std::endl;
        return n_pass == records.size() ? 0 : 1;
    }

    StdVec<HighSigmaBenchmarkRecord> records;
    records.reserve(sigma_values.size() * frequency_values.size() * 6);

    std::cout << "[p1a] high-sigma edge-flux benchmark dp=" << dp
              << " normalization_mode=" << ophelieEdgeFluxNormalizationModeName(params.edge_flux_normalization_mode_)
              << std::endl;

    for (Real frequency_hz : frequency_values)
    {
        for (Real sigma : sigma_values)
        {
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::UniformEReal, sigma, frequency_hz, dp,
                                                        center, halfsize, params));
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::ConstantAGaugeCancellationEdgeOnly, sigma,
                                                        frequency_hz, dp, center, halfsize, params));
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolve, sigma,
                                                        frequency_hz, dp, center, halfsize, params));
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::ConstantAGaugeCancellationPhiSolveSolverLocal,
                                                        sigma, frequency_hz, dp, center, halfsize, params));
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::RotationalAUniformBEdgeOnly, sigma,
                                                        frequency_hz, dp, center, halfsize, params));
            records.push_back(runHighSigmaBenchmarkCase(HighSigmaCaseKind::RotationalAUniformBPhiSolve, sigma,
                                                        frequency_hz, dp, center, halfsize, params));
        }
    }

    writeHighSigmaBenchmarkCsv(csv_path, records);

    size_t n_pass = 0;
    for (const HighSigmaBenchmarkRecord &r : records)
    {
        if (r.row_passed)
        {
            ++n_pass;
        }
    }

    const bool passed = n_pass == records.size();
    std::cout << "test_3d_ophelie_high_sigma_edge_flux_scaling rows=" << records.size() << " pass_rows=" << n_pass
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
