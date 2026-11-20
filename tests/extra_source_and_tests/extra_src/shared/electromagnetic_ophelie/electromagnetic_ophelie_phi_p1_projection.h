/**
 * @file electromagnetic_ophelie_phi_p1_projection.h
 * @brief P1 energy projection for φ (historically Doc216 E5 / Doc217 E-R): solve φ in P1 = span{1,x,y,z} by energy projection
 *        min ||Gφ + a_sign ω A||_σ² using the same SPH G as Doc214 E2.
 *
 * Avoids full-space least-squares onto truncated range(G), which incorrectly
 * projects rotational A near the boundary (Doc216 E4/E5).
 *
 * Doc217 residual semantics:
 * - remaining_field_ratio = ||Gφ+F||_σ / ||F||_σ  (NOT a linear-solver residual;
 *   may be ~1 for a correct rotational solution with φ≈0).
 * - stationarity_rel = ||M a − rhs|| / (||rhs||+eps)  (true reduced-space stationarity).
 * - Linear full-space solver residual is NOT_MEASURED on this path.
 *
 * Volume weight: current assembly uses Σ_i σ_i (·) only (no Vi). Valid as an
 * equal-volume special case on the analytic cylinder lattice; not a general
 * variable-volume quadrature.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_PHI_P1_PROJECTION_H
#define ELECTROMAGNETIC_OPHELIE_PHI_P1_PROJECTION_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_phi_variational.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "update_body_relation.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc216E5P1Report
{
    Real a[4]{0, 0, 0, 0};
    Vecd center = Vecd::Zero();
    /** Legacy name kept for log grepping; equals remaining_field_ratio. */
    Real residual_rel = 0;
    Real remaining_field_ratio = 0;
    Real stationarity_rel = 0;
    Real stationarity_abs = 0;
    bool equal_volume_sigma_only = true;
};

template <class ExecutionPolicy>
inline void execOphelieApplyE2GradientToPhiField(SolidBody &glass_body, Inner<> &inner,
                                                 const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                                 const std::string &phi_field, const std::string &grad_field)
{
    const bool use_corr = ophelieDoc214E2UseCorrectedG(params);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieE2ScalarGradientCK<Inner<>>> compute_grad(
        inner, phi_field, grad_field, names, use_corr);
    compute_grad.exec();
}

/**
 * Solve one component: min_{φ∈P1} ||Gφ + a_sign ω A||_σ², write φ into phi_field.
 * Uses grad_scratch for Gφ of basis and of the solution.
 */
template <class ExecutionPolicy>
inline OphelieDoc216E5P1Report execOphelieDoc216E5P1ProjectionOneComponent(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const std::string &phi_field, const std::string &grad_scratch, const std::string &a_field, Real a_sign)
{
    OphelieDoc216E5P1Report report;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real omega = params.omega();

    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (ophelieDoc214E2UseCorrectedG(params))
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params, Real(1.0e-6), false);
    }

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Vecd>(particles, a_field);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Vecd *A = particles.getVariableDataByName<Vecd>(a_field);
    Real *phi = particles.getVariableDataByName<Real>(phi_field);

    Vecd center = Vecd::Zero();
    for (size_t i = 0; i < n; ++i)
    {
        center += pos[i];
    }
    if (n > 0)
    {
        center /= Real(n);
    }
    report.center = center;

    // Basis columns in host buffer; also stage into phi and apply G.
    std::vector<Real> basis[4];
    for (int k = 0; k < 4; ++k)
    {
        basis[k].assign(n, Real(0));
    }
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd d = pos[i] - center;
        basis[0][i] = Real(1);
        basis[1][i] = d[0];
        basis[2][i] = d[1];
        basis[3][i] = d[2];
    }

    // GB[k][i] = (G basis_k)_i
    std::vector<Vecd> GB[4];
    for (int k = 0; k < 4; ++k)
    {
        GB[k].assign(n, Vecd::Zero());
        for (size_t i = 0; i < n; ++i)
        {
            phi[i] = basis[k][i];
        }
        syncVariableToDevice<Real>(particles, phi_field);
        execOphelieApplyE2GradientToPhiField<ExecutionPolicy>(glass_body, inner, names, params, phi_field,
                                                              grad_scratch);
        syncVariableToHost<Vecd>(particles, grad_scratch);
        Vecd *grad = particles.getVariableDataByName<Vecd>(grad_scratch);
        for (size_t i = 0; i < n; ++i)
        {
            GB[k][i] = grad[i];
        }
    }

    // F_i = a_sign * omega * A_i
    std::vector<Vecd> F(n, Vecd::Zero());
    for (size_t i = 0; i < n; ++i)
    {
        F[i] = a_sign * omega * A[i];
    }

    Eigen::Matrix<Real, 4, 4> M = Eigen::Matrix<Real, 4, 4>::Zero();
    Eigen::Matrix<Real, 4, 1> rhs = Eigen::Matrix<Real, 4, 1>::Zero();
    for (int k = 0; k < 4; ++k)
    {
        for (int l = 0; l < 4; ++l)
        {
            Real acc = 0;
            for (size_t i = 0; i < n; ++i)
            {
                acc += sigma[i] * GB[k][i].dot(GB[l][i]);
            }
            M(k, l) = acc;
        }
        Real rk = 0;
        for (size_t i = 0; i < n; ++i)
        {
            rk -= sigma[i] * GB[k][i].dot(F[i]);
        }
        rhs(k) = rk;
    }
    M(0, 0) += Real(1.0e-8); // weak pin on constant mode (legacy; Doc217 prefers exact gauge later)
    const Eigen::Matrix<Real, 4, 1> a = M.ldlt().solve(rhs);
    for (int k = 0; k < 4; ++k)
    {
        report.a[k] = a(k);
    }

    const Eigen::Matrix<Real, 4, 1> Ma_m_rhs = M * a - rhs;
    report.stationarity_abs = Ma_m_rhs.norm();
    report.stationarity_rel = report.stationarity_abs / (rhs.norm() + TinyReal);

    Real f2 = 0;
    Real r2 = 0;
    for (size_t i = 0; i < n; ++i)
    {
        phi[i] = a(0) * basis[0][i] + a(1) * basis[1][i] + a(2) * basis[2][i] + a(3) * basis[3][i];
        Vecd Gphi = a(0) * GB[0][i] + a(1) * GB[1][i] + a(2) * GB[2][i] + a(3) * GB[3][i];
        const Vecd r = Gphi + F[i];
        f2 += sigma[i] * F[i].squaredNorm();
        r2 += sigma[i] * r.squaredNorm();
    }
    report.remaining_field_ratio = std::sqrt(r2 / (f2 + TinyReal));
    report.residual_rel = report.remaining_field_ratio; // legacy alias
    syncVariableToDevice<Real>(particles, phi_field);
    return report;
}

struct OphelieDoc216E5P1ComplexSolveBundle
{
    OphelieDoc216E5P1Report imag;
    OphelieDoc216E5P1Report realp;
    bool solved_real = false;
};

template <class ExecutionPolicy>
inline OphelieDoc216E5P1ComplexSolveBundle execOphelieDoc216E5P1ProjectionSolveComplex(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    bool skip_real = false)
{
    OphelieDoc216E5P1ComplexSolveBundle out;
    std::cout << "[ophelie][doc216][e5] P1 energy projection ON "
                 "(φ∈span{1,x,y,z}; avoids full-space truncated range(G) projection)"
              << " corrected_G=" << (ophelieDoc214E2UseCorrectedG(params) ? 1 : 0)
              << " weight=sigma_only_equal_volume" << std::endl;

    out.imag = execOphelieDoc216E5P1ProjectionOneComponent<ExecutionPolicy>(
        glass_body, inner, names, params, names.phi_imag, names.grad_phi_imag, names.a_src_real, Real(1));
    std::cout << "[ophelie][doc217][e-r] imag a=[" << out.imag.a[0] << "," << out.imag.a[1] << "," << out.imag.a[2]
              << "," << out.imag.a[3] << "] remaining_field_ratio=" << out.imag.remaining_field_ratio
              << " stationarity_rel=" << out.imag.stationarity_rel
              << " linear_solver_residual=NOT_MEASURED"
              << " (legacy residual_rel=" << out.imag.residual_rel << ")" << std::endl;

    if (!skip_real && params.edge_flux_complex_)
    {
        out.realp = execOphelieDoc216E5P1ProjectionOneComponent<ExecutionPolicy>(
            glass_body, inner, names, params, names.phi_real, names.grad_phi_real, names.a_src_imag, Real(-1));
        out.solved_real = true;
        std::cout << "[ophelie][doc217][e-r] real a=[" << out.realp.a[0] << "," << out.realp.a[1] << ","
                  << out.realp.a[2] << "," << out.realp.a[3]
                  << "] remaining_field_ratio=" << out.realp.remaining_field_ratio
                  << " stationarity_rel=" << out.realp.stationarity_rel
                  << " linear_solver_residual=NOT_MEASURED"
                  << " (legacy residual_rel=" << out.realp.residual_rel << ")" << std::endl;
    }
    return out;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_P1_PROJECTION_H
