#ifndef ELECTROMAGNETIC_OPHELIE_VECTOR_DIVERGENCE_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_VECTOR_DIVERGENCE_DIAGNOSTICS_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_gradient.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "update_body_relation.h"

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

enum class OphelieVectorDivergenceMmsCase
{
    ConstantX,
    LinearXYZ,
    RotationalXY,
    QuadraticXYZ,
    SingularToroidalLegacy,
    SmoothToroidal
};

inline const char *ophelieVectorDivergenceMmsCaseName(OphelieVectorDivergenceMmsCase kind)
{
    switch (kind)
    {
    case OphelieVectorDivergenceMmsCase::ConstantX:
        return "constant_x";
    case OphelieVectorDivergenceMmsCase::LinearXYZ:
        return "linear_xyz";
    case OphelieVectorDivergenceMmsCase::RotationalXY:
        return "rotational_xy";
    case OphelieVectorDivergenceMmsCase::QuadraticXYZ:
        return "quadratic_xyz";
    case OphelieVectorDivergenceMmsCase::SingularToroidalLegacy:
        return "singular_toroidal_legacy";
    default:
        return "smooth_toroidal";
    }
}

inline Real boxDistanceToBoundary(const Vecd &position, const Vecd &center, const Vecd &halfsize)
{
    const Vecd rel = (position - center).cwiseAbs();
    return std::min(halfsize[0] - rel[0], std::min(halfsize[1] - rel[1], halfsize[2] - rel[2]));
}

inline Vecd manufacturedVectorFieldExact(OphelieVectorDivergenceMmsCase kind, const Vecd &position, const Vecd &center)
{
    const Vecd r = position - center;
    switch (kind)
    {
    case OphelieVectorDivergenceMmsCase::ConstantX:
        return Vecd(1.0, 0.0, 0.0);
    case OphelieVectorDivergenceMmsCase::LinearXYZ:
        return r;
    case OphelieVectorDivergenceMmsCase::RotationalXY:
        return Vecd(-r[1], r[0], 0.0);
    case OphelieVectorDivergenceMmsCase::QuadraticXYZ:
        return Vecd(r[0] * r[0], r[1] * r[1], r[2] * r[2]);
    case OphelieVectorDivergenceMmsCase::SmoothToroidal:
    {
        const Real r_xy_sq = r[0] * r[0] + r[1] * r[1];
        const Real f = Real(1) + Real(0.1) * r_xy_sq + Real(0.2) * r[2] * r[2];
        return Vecd(-r[1] * f, r[0] * f, 0.0);
    }
    default:
    {
        const Real r_xy_sq = r[0] * r[0] + r[1] * r[1];
        const Real r_xy = std::sqrt(std::max(r_xy_sq, TinyReal));
        const Real f = Real(1) + Real(0.1) * r_xy_sq + Real(0.2) * r[2] * r[2];
        return Vecd(-r[1] * f / r_xy, r[0] * f / r_xy, 0.0);
    }
    }
}

inline Real manufacturedDivergenceExact(OphelieVectorDivergenceMmsCase kind, const Vecd &position, const Vecd &center)
{
    const Vecd r = position - center;
    switch (kind)
    {
    case OphelieVectorDivergenceMmsCase::ConstantX:
    case OphelieVectorDivergenceMmsCase::RotationalXY:
        return Real(0);
    case OphelieVectorDivergenceMmsCase::LinearXYZ:
        return Real(3);
    case OphelieVectorDivergenceMmsCase::QuadraticXYZ:
        return Real(2) * (r[0] + r[1] + r[2]);
    case OphelieVectorDivergenceMmsCase::SmoothToroidal:
    case OphelieVectorDivergenceMmsCase::SingularToroidalLegacy:
        return Real(0);
    default:
        return Real(0);
    }
}

struct OphelieVectorDivergenceErrorMetrics
{
    size_t n_total = 0;
    size_t n_interior = 0;
    size_t n_boundary = 0;
    Real l2_all = 0.0;
    Real l2_interior = 0.0;
    Real l2_boundary = 0.0;
    Real linf_all = 0.0;
    Real signed_l2_all = 0.0;
    Real signed_l2_interior = 0.0;
    Real signed_l2_boundary = 0.0;
    Real flipped_l2_all = 0.0;
    Real flipped_l2_interior = 0.0;
    Real flipped_l2_boundary = 0.0;
    Real best_sign_l2_all = 0.0;
    Real best_sign_l2_interior = 0.0;
    Real best_sign_l2_boundary = 0.0;
    Real sign_alpha = 0.0;
    Real mean_exact_div = 0.0;
    Real mean_discrete_div = 0.0;
    Real global_integral_exact = 0.0;
    Real global_integral_discrete = 0.0;
    Real cosine_with_exact = 0.0;
};

inline OphelieVectorDivergenceErrorMetrics computeHostVectorDivergenceErrorMetrics(
    BaseParticles &particles, const Real *discrete_div, const Real *exact_div, size_t total_real_particles,
    const Vecd &center, const Vecd &halfsize, Real boundary_shell_width)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");

    OphelieVectorDivergenceErrorMetrics metrics;
    metrics.n_total = total_real_particles;

    Real diff_all_l2 = 0.0;
    Real ref_all_l2 = 0.0;
    Real diff_interior_l2 = 0.0;
    Real ref_interior_l2 = 0.0;
    Real diff_boundary_l2 = 0.0;
    Real ref_boundary_l2 = 0.0;
    Real flipped_all_l2 = 0.0;
    Real flipped_interior_l2 = 0.0;
    Real flipped_boundary_l2 = 0.0;
    Real vol_all = 0.0;
    Real vol_interior = 0.0;
    Real vol_boundary = 0.0;
    Real dot_de = 0.0;
    Real dot_ee = 0.0;

    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real signed_d = discrete_div[i] - exact_div[i];
        const Real flipped_d = -discrete_div[i] - exact_div[i];
        const Real ref = exact_div[i];
        const Real v = vol[i];
        const bool is_boundary = boxDistanceToBoundary(pos[i], center, halfsize) < boundary_shell_width;
        if (is_boundary)
        {
            metrics.n_boundary++;
            diff_boundary_l2 += v * signed_d * signed_d;
            flipped_boundary_l2 += v * flipped_d * flipped_d;
            ref_boundary_l2 += v * ref * ref;
            vol_boundary += v;
        }
        else
        {
            metrics.n_interior++;
            diff_interior_l2 += v * signed_d * signed_d;
            flipped_interior_l2 += v * flipped_d * flipped_d;
            ref_interior_l2 += v * ref * ref;
            vol_interior += v;
        }
        diff_all_l2 += v * signed_d * signed_d;
        flipped_all_l2 += v * flipped_d * flipped_d;
        ref_all_l2 += v * ref * ref;
        dot_de += v * discrete_div[i] * exact_div[i];
        dot_ee += v * exact_div[i] * exact_div[i];
        metrics.linf_all = std::max(metrics.linf_all, std::abs(signed_d));
        metrics.mean_exact_div += v * exact_div[i];
        metrics.mean_discrete_div += v * discrete_div[i];
        metrics.global_integral_exact += v * exact_div[i];
        metrics.global_integral_discrete += v * discrete_div[i];
        vol_all += v;
    }

    const auto normalized_l2 = [](Real diff_l2, Real ref_l2, Real region_vol)
    {
        const Real err = std::sqrt(diff_l2);
        const Real ref = std::sqrt(ref_l2);
        const Real vol_scale = std::sqrt(std::max(region_vol, TinyReal));
        if (ref > Real(1.0e-6) * vol_scale)
        {
            return err / ref;
        }
        return err / vol_scale;
    };

    metrics.l2_all = normalized_l2(diff_all_l2, ref_all_l2, vol_all);
    metrics.l2_interior = normalized_l2(diff_interior_l2, ref_interior_l2, vol_interior);
    metrics.l2_boundary = normalized_l2(diff_boundary_l2, ref_boundary_l2, vol_boundary);
    metrics.signed_l2_all = metrics.l2_all;
    metrics.signed_l2_interior = metrics.l2_interior;
    metrics.signed_l2_boundary = metrics.l2_boundary;
    metrics.flipped_l2_all = normalized_l2(flipped_all_l2, ref_all_l2, vol_all);
    metrics.flipped_l2_interior = normalized_l2(flipped_interior_l2, ref_interior_l2, vol_interior);
    metrics.flipped_l2_boundary = normalized_l2(flipped_boundary_l2, ref_boundary_l2, vol_boundary);
    metrics.best_sign_l2_all = std::min(metrics.signed_l2_all, metrics.flipped_l2_all);
    metrics.best_sign_l2_interior = std::min(metrics.signed_l2_interior, metrics.flipped_l2_interior);
    metrics.best_sign_l2_boundary = std::min(metrics.signed_l2_boundary, metrics.flipped_l2_boundary);
    metrics.sign_alpha = dot_de / (dot_ee + TinyReal);
    metrics.mean_exact_div /= (vol_all + TinyReal);
    metrics.mean_discrete_div /= (vol_all + TinyReal);
    metrics.cosine_with_exact = hostVolWeightedCosine(particles, discrete_div, exact_div, total_real_particles);
    return metrics;
}

inline void assignManufacturedVectorFieldToParticles(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                     OphelieVectorDivergenceMmsCase kind, const Vecd &center)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_field = particles.getVariableDataByName<Vecd>(names.a_src_real);
    for (size_t i = 0; i < n; ++i)
    {
        a_field[i] = manufacturedVectorFieldExact(kind, pos[i], center);
    }
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
}

template <class ExecutionPolicy>
inline OphelieVectorDivergenceErrorMetrics evaluateOphelieVectorDivergenceMmsCase(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieVectorDivergenceMmsCase kind,
    const Vecd &center, const Vecd &halfsize, Real dp, bool use_corrected_divergence)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real boundary_shell_width = Real(2) * dp;

    assignManufacturedVectorFieldToParticles(particles, names, kind, center);

    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();

    if (use_corrected_divergence)
    {
        execOphelieVecdDivergenceCorrected<ExecutionPolicy>(glass_body, inner, names, names.a_src_real,
                                                            names.div_j_imag);
    }
    else
    {
        execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.a_src_real, names.div_j_imag);
    }

    syncVariableToHost<Real>(particles, names.div_j_imag);
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *discrete_div = particles.getVariableDataByName<Real>(names.div_j_imag);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");

    StdVec<Real> exact_div(n, Real(0));
    for (size_t i = 0; i < n; ++i)
    {
        exact_div[i] = manufacturedDivergenceExact(kind, pos[i], center);
    }
    return computeHostVectorDivergenceErrorMetrics(particles, discrete_div, exact_div.data(), n, center, halfsize,
                                                   boundary_shell_width);
}

struct OphelieVectorDivergenceMmsCaseReport
{
    OphelieVectorDivergenceMmsCase kind = OphelieVectorDivergenceMmsCase::ConstantX;
    OphelieVectorDivergenceErrorMetrics uncorrected;
    OphelieVectorDivergenceErrorMetrics corrected;
    Real corrected_vs_uncorrected_l2 = 0.0;
    Real corrected_vs_uncorrected_cosine = 0.0;
};

template <class ExecutionPolicy>
inline OphelieVectorDivergenceMmsCaseReport evaluateOphelieVectorDivergenceMmsCaseBothOperators(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieVectorDivergenceMmsCase kind,
    const Vecd &center, const Vecd &halfsize, Real dp)
{
    OphelieVectorDivergenceMmsCaseReport report;
    report.kind = kind;
    report.uncorrected =
        evaluateOphelieVectorDivergenceMmsCase<ExecutionPolicy>(glass_body, inner, names, kind, center, halfsize, dp,
                                                                false);

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> div_unc(n, Real(0));
    hostReadScalarField(particles, names.div_j_imag, div_unc.data(), n);

    report.corrected =
        evaluateOphelieVectorDivergenceMmsCase<ExecutionPolicy>(glass_body, inner, names, kind, center, halfsize, dp,
                                                                true);
    StdVec<Real> div_corr(n, Real(0));
    hostReadScalarField(particles, names.div_j_imag, div_corr.data(), n);

    report.corrected_vs_uncorrected_l2 =
        hostVolWeightedRelativeDifference(particles, div_unc.data(), div_corr.data(), n);
    report.corrected_vs_uncorrected_cosine = hostVolWeightedCosine(particles, div_unc.data(), div_corr.data(), n);
    return report;
}

inline void logOphelieVectorDivergenceOperatorMetrics(const char *operator_kind,
                                                      const OphelieVectorDivergenceErrorMetrics &metrics)
{
    std::cout << " operator=" << operator_kind << " l2_interior=" << metrics.l2_interior
              << " flipped_l2_interior=" << metrics.flipped_l2_interior
              << " best_sign_l2_interior=" << metrics.best_sign_l2_interior << " sign_alpha=" << metrics.sign_alpha
              << " mean_exact=" << metrics.mean_exact_div << " mean_discrete=" << metrics.mean_discrete_div
              << " cosine=" << metrics.cosine_with_exact;
}

inline void logOphelieVectorDivergenceMmsCaseReport(Real dp, const OphelieVectorDivergenceMmsCaseReport &report)
{
    const char *case_name = ophelieVectorDivergenceMmsCaseName(report.kind);
    std::cout << "[ophelie] vector_div_mms dp=" << dp << " case=" << case_name << " n_total=" << report.uncorrected.n_total
              << " n_interior=" << report.uncorrected.n_interior << " n_boundary=" << report.uncorrected.n_boundary;
    logOphelieVectorDivergenceOperatorMetrics("D_unc", report.uncorrected);
    logOphelieVectorDivergenceOperatorMetrics("D_c", report.corrected);
    std::cout << " D_c_vs_D_l2=" << report.corrected_vs_uncorrected_l2
              << " D_c_vs_D_cosine=" << report.corrected_vs_uncorrected_cosine << std::endl;

    if (report.kind == OphelieVectorDivergenceMmsCase::LinearXYZ &&
        report.uncorrected.l2_interior > Real(1.5) && report.uncorrected.flipped_l2_interior < Real(0.2) &&
        report.uncorrected.sign_alpha < Real(-0.5))
    {
        std::cout << "[ophelie][warning] divergence sign convention appears flipped for linear_xyz MMS case."
                  << std::endl;
    }
}

inline void appendOphelieVectorDivergenceMmsCsv(const std::string &path, Real dp,
                                                const OphelieVectorDivergenceMmsCaseReport &report)
{
    namespace fs = std::filesystem;
    const bool write_header = path.empty() || !fs::exists(path) || fs::file_size(path) == 0;
    std::ofstream csv(path, std::ios::app);
    if (!csv)
    {
        return;
    }
    if (write_header)
    {
        csv << "dp,case_name,operator_kind,n_total,n_interior,l2_all,l2_interior,l2_boundary,linf_all,"
               "flipped_l2_all,flipped_l2_interior,flipped_l2_boundary,best_sign_l2_all,best_sign_l2_interior,"
               "best_sign_l2_boundary,sign_alpha,mean_exact_div,mean_discrete_div,global_integral_exact,"
               "global_integral_discrete,cosine_with_exact\n";
    }
    const char *case_name = ophelieVectorDivergenceMmsCaseName(report.kind);
    auto write_row = [&](const char *operator_kind, const OphelieVectorDivergenceErrorMetrics &metrics)
    {
        csv << dp << "," << case_name << "," << operator_kind << "," << metrics.n_total << "," << metrics.n_interior
            << "," << metrics.l2_all << "," << metrics.l2_interior << "," << metrics.l2_boundary << "," << metrics.linf_all
            << "," << metrics.flipped_l2_all << "," << metrics.flipped_l2_interior << "," << metrics.flipped_l2_boundary
            << "," << metrics.best_sign_l2_all << "," << metrics.best_sign_l2_interior << "," << metrics.best_sign_l2_boundary
            << "," << metrics.sign_alpha << "," << metrics.mean_exact_div << "," << metrics.mean_discrete_div << ","
            << metrics.global_integral_exact << "," << metrics.global_integral_discrete << "," << metrics.cosine_with_exact
            << "\n";
    };
    write_row("D_unc", report.uncorrected);
    write_row("D_c", report.corrected);
}

inline Real cylinderDistanceToBoundary(const Vecd &position, const Vecd &center, Real radius, Real half_height)
{
    const Vecd rel = position - center;
    const Real r_xy = std::sqrt(rel[0] * rel[0] + rel[1] * rel[1]);
    const Real radial_gap = radius - r_xy;
    const Real axial_gap = half_height - std::abs(rel[2]);
    return std::min(radial_gap, axial_gap);
}

struct OphelieScalarFieldNormMetrics
{
    size_t n_total = 0;
    size_t n_interior = 0;
    size_t n_boundary = 0;
    Real l2_all = 0.0;
    Real l2_interior = 0.0;
    Real l2_boundary = 0.0;
    Real linf_all = 0.0;
    Real mean = 0.0;
    Real rms = 0.0;
    Real volume_integral = 0.0;
};

inline OphelieScalarFieldNormMetrics computeHostScalarFieldNormMetrics(
    BaseParticles &particles, const Real *field, size_t total_real_particles, const Vecd &center, Real radius,
    Real half_height, Real boundary_shell_width)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");

    OphelieScalarFieldNormMetrics metrics;
    metrics.n_total = total_real_particles;

    Real sum_all_l2 = 0.0;
    Real sum_interior_l2 = 0.0;
    Real sum_boundary_l2 = 0.0;
    Real total_vol = 0.0;

    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real value = field[i];
        const Real v = vol[i];
        const bool is_boundary =
            cylinderDistanceToBoundary(pos[i], center, radius, half_height) < boundary_shell_width;
        if (is_boundary)
        {
            metrics.n_boundary++;
            sum_boundary_l2 += v * value * value;
        }
        else
        {
            metrics.n_interior++;
            sum_interior_l2 += v * value * value;
        }
        sum_all_l2 += v * value * value;
        metrics.linf_all = std::max(metrics.linf_all, std::abs(value));
        metrics.mean += v * value;
        metrics.volume_integral += v * value;
        total_vol += v;
    }

    metrics.l2_all = std::sqrt(sum_all_l2);
    metrics.l2_interior = std::sqrt(sum_interior_l2);
    metrics.l2_boundary = std::sqrt(sum_boundary_l2);
    metrics.mean /= (total_vol + TinyReal);
    metrics.rms = std::sqrt(sum_all_l2 / (total_vol + TinyReal));
    return metrics;
}

struct OphelieBiotVectorFieldDecomposition
{
    Real a_r_rms = 0.0;
    Real a_theta_rms = 0.0;
    Real a_z_rms = 0.0;
    Real a_total_rms = 0.0;
    Real a_theta_fraction = 0.0;
};

inline OphelieBiotVectorFieldDecomposition computeHostCylindricalAFieldDecomposition(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles, const Vecd &center)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    Real sum_ar = 0.0;
    Real sum_ath = 0.0;
    Real sum_az = 0.0;
    Real sum_atot = 0.0;
    Real total_vol = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Vecd rel = pos[i] - center;
        const Real r_xy = std::sqrt(rel[0] * rel[0] + rel[1] * rel[1]);
        const Real inv_r = Real(1) / std::max(r_xy, TinyReal);
        const Real cos_t = rel[0] * inv_r;
        const Real sin_t = rel[1] * inv_r;
        const Vecd &a = a_src[i];
        const Real a_r = a[0] * cos_t + a[1] * sin_t;
        const Real a_theta = -a[0] * sin_t + a[1] * cos_t;
        const Real a_z = a[2];
        const Real v = vol[i];
        sum_ar += v * a_r * a_r;
        sum_ath += v * a_theta * a_theta;
        sum_az += v * a_z * a_z;
        sum_atot += v * a.squaredNorm();
        total_vol += v;
    }

    OphelieBiotVectorFieldDecomposition decomposition;
    decomposition.a_r_rms = std::sqrt(sum_ar / (total_vol + TinyReal));
    decomposition.a_theta_rms = std::sqrt(sum_ath / (total_vol + TinyReal));
    decomposition.a_z_rms = std::sqrt(sum_az / (total_vol + TinyReal));
    decomposition.a_total_rms = std::sqrt(sum_atot / (total_vol + TinyReal));
    decomposition.a_theta_fraction = decomposition.a_theta_rms / (decomposition.a_total_rms + TinyReal);
    return decomposition;
}

struct OphelieBiotSigmaADivergenceDiagnostics
{
    OphelieScalarFieldNormMetrics uncorrected;
    OphelieScalarFieldNormMetrics corrected;
    Real unc_vs_comp_l2 = 0.0;
    Real unc_vs_comp_cosine = 0.0;
    Real normalized_div_unc = 0.0;
    Real normalized_div_comp = 0.0;
    OphelieBiotVectorFieldDecomposition a_decomposition;
};

template <class ExecutionPolicy>
inline OphelieBiotSigmaADivergenceDiagnostics evaluateOphelieBiotSigmaADivergenceDiagnostics(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const Vecd &center, Real radius,
    Real half_height, Real dp)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real boundary_shell_width = Real(2) * dp;
    const Real length_scale = std::max(radius, half_height);

    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdFieldBySigmaCK> scale_sigma_a(
        glass_body, names, names.a_src_real, names.j_imag);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    scale_sigma_a.exec();

    execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, names.j_imag, names.div_j_imag);
    StdVec<Real> div_unc(n, Real(0));
    hostReadScalarField(particles, names.div_j_imag, div_unc.data(), n);

    execOphelieVecdDivergenceCorrected<ExecutionPolicy>(glass_body, inner, names, names.j_imag, names.div_j_imag);
    StdVec<Real> div_comp(n, Real(0));
    hostReadScalarField(particles, names.div_j_imag, div_comp.data(), n);

    OphelieBiotSigmaADivergenceDiagnostics diagnostics;
    diagnostics.uncorrected =
        computeHostScalarFieldNormMetrics(particles, div_unc.data(), n, center, radius, half_height, boundary_shell_width);
    diagnostics.corrected = computeHostScalarFieldNormMetrics(particles, div_comp.data(), n, center, radius,
                                                              half_height, boundary_shell_width);
    diagnostics.unc_vs_comp_l2 = hostVolWeightedRelativeDifference(particles, div_unc.data(), div_comp.data(), n);
    diagnostics.unc_vs_comp_cosine = hostVolWeightedCosine(particles, div_unc.data(), div_comp.data(), n);
    diagnostics.a_decomposition = computeHostCylindricalAFieldDecomposition(particles, names, n, center);

    syncVariableToHost<Vecd>(particles, names.j_imag);
    const Vecd *sigma_a = particles.getVariableDataByName<Vecd>(names.j_imag);
    Real sigma_a_norm = 0.0;
    Real total_vol = 0.0;
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < n; ++i)
    {
        sigma_a_norm += vol[i] * sigma_a[i].squaredNorm();
        total_vol += vol[i];
    }
    const Real sigma_a_rms = std::sqrt(sigma_a_norm / (total_vol + TinyReal));
    const Real denom = sigma_a_rms / (length_scale + TinyReal);
    diagnostics.normalized_div_unc = diagnostics.uncorrected.rms / (denom + TinyReal);
    diagnostics.normalized_div_comp = diagnostics.corrected.rms / (denom + TinyReal);
    return diagnostics;
}

inline void logOphelieBiotSigmaADivergenceDiagnostics(const OphelieBiotSigmaADivergenceDiagnostics &diagnostics)
{
    std::cout << "[ophelie] biot_sigma_a_div: D_unc_l2_all=" << diagnostics.uncorrected.l2_all
              << " D_unc_l2_interior=" << diagnostics.uncorrected.l2_interior
              << " D_unc_l2_boundary=" << diagnostics.uncorrected.l2_boundary
              << " D_unc_rms=" << diagnostics.uncorrected.rms << " D_c_l2_all=" << diagnostics.corrected.l2_all
              << " D_c_l2_interior=" << diagnostics.corrected.l2_interior
              << " D_c_l2_boundary=" << diagnostics.corrected.l2_boundary << " D_c_rms=" << diagnostics.corrected.rms
              << " D_unc_vs_D_c_l2=" << diagnostics.unc_vs_comp_l2
              << " D_unc_vs_D_c_cosine=" << diagnostics.unc_vs_comp_cosine
              << " normalized_div_unc=" << diagnostics.normalized_div_unc
              << " normalized_div_comp=" << diagnostics.normalized_div_comp
              << " A_theta_fraction=" << diagnostics.a_decomposition.a_theta_fraction
              << " A_r_rms=" << diagnostics.a_decomposition.a_r_rms
              << " A_theta_rms=" << diagnostics.a_decomposition.a_theta_rms
              << " A_z_rms=" << diagnostics.a_decomposition.a_z_rms << std::endl;
}

inline void appendOphelieBiotSigmaADivergenceCsv(const std::string &path, Real dp,
                                                 const OphelieBiotSigmaADivergenceDiagnostics &diagnostics)
{
    namespace fs = std::filesystem;
    const bool write_header = path.empty() || !fs::exists(path) || fs::file_size(path) == 0;
    std::ofstream csv(path, std::ios::app);
    if (!csv)
    {
        return;
    }
    if (write_header)
    {
        csv << "dp,operator_kind,l2_all,l2_interior,l2_boundary,linf_all,mean,rms,volume_integral,"
               "normalized_div,unc_vs_comp_l2,unc_vs_comp_cosine,a_theta_fraction,a_r_rms,a_theta_rms,a_z_rms\n";
    }
    auto write_row = [&](const char *operator_kind, const OphelieScalarFieldNormMetrics &metrics, Real normalized_div)
    {
        csv << dp << "," << operator_kind << "," << metrics.l2_all << "," << metrics.l2_interior << ","
            << metrics.l2_boundary << "," << metrics.linf_all << "," << metrics.mean << "," << metrics.rms << ","
            << metrics.volume_integral << "," << normalized_div << "," << diagnostics.unc_vs_comp_l2 << ","
            << diagnostics.unc_vs_comp_cosine << "," << diagnostics.a_decomposition.a_theta_fraction << ","
            << diagnostics.a_decomposition.a_r_rms << "," << diagnostics.a_decomposition.a_theta_rms << ","
            << diagnostics.a_decomposition.a_z_rms << "\n";
    };
    write_row("D_unc", diagnostics.uncorrected, diagnostics.normalized_div_unc);
    write_row("D_c", diagnostics.corrected, diagnostics.normalized_div_comp);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_VECTOR_DIVERGENCE_DIAGNOSTICS_H
