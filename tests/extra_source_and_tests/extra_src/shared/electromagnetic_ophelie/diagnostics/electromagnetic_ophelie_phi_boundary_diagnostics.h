#ifndef ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_DIAGNOSTICS_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OpheliePhiRhsCompatibilityMetrics
{
    Real rhs_volume_integral = 0.0;
    Real rhs_volume_mean = 0.0;
    Real rhs_vol_l2 = 0.0;
    Real total_volume = 0.0;
};

inline OpheliePhiRhsCompatibilityMetrics computeOpheliePhiRhsCompatibilityMetrics(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OpheliePhiRhsCompatibilityMetrics metrics;
    Real rhs_l2_sq = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        metrics.rhs_volume_integral += vol[i] * rhs[i];
        metrics.total_volume += vol[i];
        rhs_l2_sq += vol[i] * rhs[i] * rhs[i];
    }
    metrics.rhs_volume_mean = metrics.rhs_volume_integral / (metrics.total_volume + TinyReal);
    metrics.rhs_vol_l2 = std::sqrt(rhs_l2_sq);
    return metrics;
}

/** b_i <- b_i - sum(b_j V_j) / sum(V_j). Returns mean removed. */
inline Real applyOpheliePhiRhsZeroMeanProjection(BaseParticles &particles, const std::string &rhs_field_name,
                                                 size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, rhs_field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    Real *rhs = particles.getVariableDataByName<Real>(rhs_field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    Real volume_integral = 0.0;
    Real total_volume = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        volume_integral += vol[i] * rhs[i];
        total_volume += vol[i];
    }
    const Real mean = volume_integral / (total_volume + TinyReal);
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        rhs[i] -= mean;
    }
    syncVariableToDevice<Real>(particles, rhs_field_name);
    return mean;
}

inline void finalizeOpheliePhiImagRhsHost(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                          OphelieParameters &params, const OpheliePhiBoundaryGeometryContext *geom,
                                          Real dp, OpheliePhiNeumannRhsCorrectionStats *neumann_stats = nullptr)
{
    if (geom != nullptr && params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
    {
        const OpheliePhiNeumannRhsCorrectionStats stats =
            applyOpheliePhiOneSidedNeumannRhsCorrection(particles, names, params, *geom, dp);
        if (neumann_stats != nullptr)
        {
            *neumann_stats = stats;
        }
    }
    if (params.phi_rhs_project_zero_mean_)
    {
        applyOpheliePhiRhsZeroMeanProjection(particles, names.phi_rhs_imag, particles.TotalRealParticles());
    }
}

struct OpheliePhiBoundaryJnMetrics
{
    size_t n_boundary = 0;
    Real boundary_vol_fraction = 0.0;
    Real jn_boundary_weighted_l2 = 0.0;
    Real j_boundary_weighted_l2 = 0.0;
    Real jn_boundary_rel = 0.0;
    Real jn_boundary_max = 0.0;
    /** Mean |n·(σA)| on boundary (Neumann target magnitude reference). */
    Real n_sigma_a_boundary_mean = 0.0;
    Real n_sigma_a_boundary_max = 0.0;
};

inline OpheliePhiBoundaryJnMetrics computeOpheliePhiBoundaryJnMetrics(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles,
    const OpheliePhiBoundaryGeometryContext &geom, const OphelieParameters &params, Real boundary_width)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.j_imag);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(names.j_imag);
    const Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real omega = params.omega();

    OpheliePhiBoundaryJnMetrics metrics;
    Real jn_l2_sq = 0.0;
    Real j_l2_sq = 0.0;
    Real total_vol = 0.0;
    Real boundary_vol = 0.0;
    Real n_sigma_a_sum = 0.0;
    Real n_sigma_a_weight = 0.0;

    for (size_t i = 0; i < total_real_particles; ++i)
    {
        total_vol += vol[i];
        const Real dist = boundaryDistanceFromContext(pos[i], geom);
        if (dist > boundary_width)
        {
            continue;
        }
        ++metrics.n_boundary;
        boundary_vol += vol[i];
        const Vecd n_out = boundaryOutwardNormalFromContext(pos[i], geom);
        const Real jn = j_imag[i].dot(n_out);
        const Real n_sigma_a = n_out.dot(sigma[i] * a_src[i]);
        metrics.jn_boundary_max = std::max(metrics.jn_boundary_max, std::abs(jn));
        metrics.n_sigma_a_boundary_max = std::max(metrics.n_sigma_a_boundary_max, std::abs(n_sigma_a));
        jn_l2_sq += vol[i] * jn * jn;
        j_l2_sq += vol[i] * j_imag[i].squaredNorm();
        n_sigma_a_sum += vol[i] * std::abs(n_sigma_a);
        n_sigma_a_weight += vol[i];
        (void)omega;
    }

    metrics.boundary_vol_fraction = boundary_vol / (total_vol + TinyReal);
    metrics.jn_boundary_weighted_l2 = std::sqrt(jn_l2_sq);
    metrics.j_boundary_weighted_l2 = std::sqrt(j_l2_sq);
    metrics.jn_boundary_rel =
        metrics.jn_boundary_weighted_l2 / (metrics.j_boundary_weighted_l2 + TinyReal);
    metrics.n_sigma_a_boundary_mean = n_sigma_a_sum / (n_sigma_a_weight + TinyReal);
    return metrics;
}

inline OpheliePhiBoundaryJnMetrics computeBoxBoundaryJnMetrics(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles,
    const OphelieParameters &params, const Vecd &center, const Vecd &halfsize, Real boundary_width)
{
    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = halfsize;
    return computeOpheliePhiBoundaryJnMetrics(particles, names, total_real_particles, geom, params, boundary_width);
}

inline OpheliePhiBoundaryJnMetrics computeFrenchCylinderBoundaryJnMetrics(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles,
    const OphelieFrenchReducedCaseParams &french, const OphelieParameters &params, Real boundary_width)
{
    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
    geom.french = french;
    return computeOpheliePhiBoundaryJnMetrics(particles, names, total_real_particles, geom, params, boundary_width);
}

struct OpheliePhiP0DiagnosticResult
{
    OpheliePhiRhsCompatibilityMetrics rhs;
    bool rhs_zero_mean_projected = false;
    Real rhs_mean_removed = 0.0;
    Real phi_eq_res_vol_pre_solve = 0.0;
    Real phi_eq_res_vol_pre_solve_after_rhs_proj = 0.0;
    OpheliePhiBoundaryJnMetrics boundary_jn_level0;
    OpheliePhiBoundaryJnMetrics boundary_jn_post_phi;
    Real phi_eq_res_vol_post_solve = 0.0;
    Real div_j_l2_reduction = 0.0;
    Real joule_power_raw = 0.0;
    Real neumann_correction_l2 = 0.0;
    size_t neumann_n_boundary = 0;
};

inline void logOpheliePhiP0Diagnostics(const OpheliePhiP0DiagnosticResult &diag, const OphelieParameters &params)
{
    std::cout << "[ophelie] phi_p0_rhs: integral=" << diag.rhs.rhs_volume_integral
              << " mean=" << diag.rhs.rhs_volume_mean << " l2=" << diag.rhs.rhs_vol_l2
              << " project_zero_mean=" << (params.phi_rhs_project_zero_mean_ ? 1 : 0)
              << " boundary_mode=" << phiBoundaryModeName(params.phi_boundary_mode_);
    if (params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
    {
        std::cout << " neumann_corr_l2=" << diag.neumann_correction_l2
                  << " neumann_n_boundary=" << diag.neumann_n_boundary;
    }
    std::cout << " eq_res_pre_solve=" << diag.phi_eq_res_vol_pre_solve << std::endl;

    if (params.phi_boundary_diagnostics_)
    {
        std::cout << "[ophelie] phi_p0_boundary_jn_level0: n_boundary=" << diag.boundary_jn_level0.n_boundary
                  << " vol_frac=" << diag.boundary_jn_level0.boundary_vol_fraction
                  << " Jn_l2=" << diag.boundary_jn_level0.jn_boundary_weighted_l2
                  << " J_l2=" << diag.boundary_jn_level0.j_boundary_weighted_l2
                  << " Jn_rel=" << diag.boundary_jn_level0.jn_boundary_rel
                  << " Jn_max=" << diag.boundary_jn_level0.jn_boundary_max
                  << " |n·σA|_mean=" << diag.boundary_jn_level0.n_sigma_a_boundary_mean << std::endl;
        std::cout << "[ophelie] phi_p0_boundary_jn_post_phi: n_boundary=" << diag.boundary_jn_post_phi.n_boundary
                  << " vol_frac=" << diag.boundary_jn_post_phi.boundary_vol_fraction
                  << " Jn_l2=" << diag.boundary_jn_post_phi.jn_boundary_weighted_l2
                  << " J_l2=" << diag.boundary_jn_post_phi.j_boundary_weighted_l2
                  << " Jn_rel=" << diag.boundary_jn_post_phi.jn_boundary_rel
                  << " Jn_max=" << diag.boundary_jn_post_phi.jn_boundary_max
                  << " |n·σA|_mean=" << diag.boundary_jn_post_phi.n_sigma_a_boundary_mean << std::endl;
    }
}

inline std::string opheliePhiP0CaseLabel(const OphelieParameters &params)
{
    const bool proj = params.phi_rhs_project_zero_mean_;
    const bool neumann = params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann;
    const bool lhs_grad = params.phi_boundary_lhs_grad_neumann_;
    const bool grad_corr = params.phi_gradient_correction_;
    const bool compatible = params.phi_compatible_correction_;
    if (compatible)
    {
        if (proj && neumann)
        {
            return "D_rhs_proj_neumann_compatible";
        }
        if (neumann)
        {
            return "C_neumann_compatible";
        }
        if (proj)
        {
            return "B_rhs_proj_compatible";
        }
        return "F_compatible";
    }
    if (proj && neumann)
    {
        if (grad_corr)
        {
            return lhs_grad ? "D_rhs_proj_neumann_lhs_gradcorr" : "D_rhs_proj_neumann_gradcorr";
        }
        return lhs_grad ? "D_rhs_proj_neumann_lhs" : "D_rhs_proj_neumann";
    }
    if (neumann)
    {
        if (grad_corr)
        {
            return lhs_grad ? "C_neumann_lhs_gradcorr" : "C_neumann_gradcorr";
        }
        return lhs_grad ? "C_neumann_lhs" : "C_neumann";
    }
    if (grad_corr)
    {
        return proj ? "B_rhs_proj_gradcorr" : "E_gradcorr";
    }
    if (proj)
    {
        return "B_rhs_proj";
    }
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::LegacyPairwise &&
        params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux)
    {
        if (params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux)
        {
            return "H_edge_flux_production";
        }
        return "G_edge_flux";
    }
    return "A_baseline";
}

inline void appendOpheliePhiP0DiagnosticCsv(const std::string &path, const std::string &case_label,
                                            const std::string &particle_source, size_t n, Real dp,
                                            const OphelieParameters &params, const OpheliePhiP0DiagnosticResult &diag)
{
    namespace fs = std::filesystem;
    if (!path.empty())
    {
        const fs::path parent = fs::path(path).parent_path();
        if (!parent.empty())
        {
            fs::create_directories(parent);
        }
    }

    std::ofstream csv(path, std::ios::app);
    if (!csv)
    {
        std::cerr << "[ophelie] phi_p0_csv: cannot open " << path << std::endl;
        return;
    }
    csv.seekp(0, std::ios::end);
    const bool write_header = csv.tellp() == 0;
    if (write_header)
    {
        csv << "case_label,particle_source,n,dp,phi_boundary_mode,rhs_volume_integral,rhs_volume_mean,rhs_vol_l2,"
               "rhs_zero_mean_projected,neumann_correction_l2,phi_eq_res_pre_solve,phi_eq_res_post_solve,"
               "div_j_l2_reduction,joule_power_raw,"
               "jn_level0_rel,jn_level0_max,jn_post_phi_rel,jn_post_phi_max,"
               "n_sigma_a_level0_mean,n_sigma_a_post_phi_mean\n";
    }
    csv << case_label << "," << particle_source << "," << n << "," << dp << ","
        << phiBoundaryModeName(params.phi_boundary_mode_) << "," << diag.rhs.rhs_volume_integral << ","
        << diag.rhs.rhs_volume_mean << "," << diag.rhs.rhs_vol_l2 << "," << (diag.rhs_zero_mean_projected ? 1 : 0)
        << "," << diag.neumann_correction_l2 << "," << diag.phi_eq_res_vol_pre_solve << ","
        << diag.phi_eq_res_vol_post_solve << "," << diag.div_j_l2_reduction << "," << diag.joule_power_raw << ","
        << diag.boundary_jn_level0.jn_boundary_rel << "," << diag.boundary_jn_level0.jn_boundary_max << ","
        << diag.boundary_jn_post_phi.jn_boundary_rel << "," << diag.boundary_jn_post_phi.jn_boundary_max << ","
        << diag.boundary_jn_level0.n_sigma_a_boundary_mean << ","
        << diag.boundary_jn_post_phi.n_sigma_a_boundary_mean << "\n";
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_DIAGNOSTICS_H
