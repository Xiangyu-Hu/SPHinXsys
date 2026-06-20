#ifndef ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_H
#define ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OpheliePhiBoundaryGeometryContext
{
    OpheliePhiBoundaryNormalSource normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    Vecd box_center = Vecd(0.0, 0.0, 0.0);
    Vecd box_halfsize = Vecd(0.0, 0.0, 0.0);
    OphelieFrenchReducedCaseParams french;
};

inline Real analyticBoxBoundaryDistance(const Vecd &pos, const Vecd &center, const Vecd &halfsize)
{
    const Vecd r = pos - center;
    return std::min({halfsize[0] - std::abs(r[0]), halfsize[1] - std::abs(r[1]), halfsize[2] - std::abs(r[2])});
}

inline Vecd analyticBoxOutwardNormal(const Vecd &pos, const Vecd &center, const Vecd &halfsize)
{
    const Vecd r = pos - center;
    const Real dist_x = halfsize[0] - std::abs(r[0]);
    const Real dist_y = halfsize[1] - std::abs(r[1]);
    const Real dist_z = halfsize[2] - std::abs(r[2]);
    if (dist_x <= dist_y && dist_x <= dist_z)
    {
        return Vecd(r[0] >= 0.0 ? 1.0 : -1.0, 0.0, 0.0);
    }
    if (dist_y <= dist_z)
    {
        return Vecd(0.0, r[1] >= 0.0 ? 1.0 : -1.0, 0.0);
    }
    return Vecd(0.0, 0.0, r[2] >= 0.0 ? 1.0 : -1.0);
}

inline Real frenchReducedCylinderBoundaryDistance(const Vecd &pos, const OphelieFrenchReducedCaseParams &french)
{
    const Vecd &center = french.glass_center;
    const Real z_lo = center[2] - french.glass_half_height;
    const Real z_hi = center[2] + french.glass_half_height;
    const Real dx = pos[0] - center[0];
    const Real dy = pos[1] - center[1];
    const Real r = std::sqrt(dx * dx + dy * dy);
    const Real dist_side = french.glass_radius - r;
    const Real dist_bottom = pos[2] - z_lo;
    const Real dist_top = z_hi - pos[2];
    return std::min(dist_side, std::min(dist_bottom, dist_top));
}

inline Vecd frenchReducedCylinderOutwardNormal(const Vecd &pos, const OphelieFrenchReducedCaseParams &french)
{
    const Vecd &center = french.glass_center;
    const Real z_lo = center[2] - french.glass_half_height;
    const Real z_hi = center[2] + french.glass_half_height;
    const Real dx = pos[0] - center[0];
    const Real dy = pos[1] - center[1];
    const Real r = std::sqrt(dx * dx + dy * dy);
    const Real dist_side = french.glass_radius - r;
    const Real dist_bottom = pos[2] - z_lo;
    const Real dist_top = z_hi - pos[2];
    if (dist_side <= dist_bottom && dist_side <= dist_top)
    {
        if (r > TinyReal)
        {
            return Vecd(dx / r, dy / r, 0.0);
        }
        return Vecd(0.0, 0.0, 1.0);
    }
    if (dist_bottom <= dist_top)
    {
        return Vecd(0.0, 0.0, -1.0);
    }
    return Vecd(0.0, 0.0, 1.0);
}

inline Real boundaryDistanceFromContext(const Vecd &pos, const OpheliePhiBoundaryGeometryContext &geom)
{
    switch (geom.normal_source)
    {
    case OpheliePhiBoundaryNormalSource::AnalyticCylinder:
        return frenchReducedCylinderBoundaryDistance(pos, geom.french);
    default:
        return analyticBoxBoundaryDistance(pos, geom.box_center, geom.box_halfsize);
    }
}

inline Vecd boundaryOutwardNormalFromContext(const Vecd &pos, const OpheliePhiBoundaryGeometryContext &geom)
{
    switch (geom.normal_source)
    {
    case OpheliePhiBoundaryNormalSource::AnalyticCylinder:
        return frenchReducedCylinderOutwardNormal(pos, geom.french);
    default:
        return analyticBoxOutwardNormal(pos, geom.box_center, geom.box_halfsize);
    }
}

struct OpheliePhiNeumannRhsCorrectionStats
{
    size_t n_boundary = 0;
    Real correction_l2 = 0.0;
    Real correction_max = 0.0;
};

inline bool opheliePhiBoundaryGradNeumannInLhs(const OphelieParameters &params)
{
    return params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann &&
           params.phi_boundary_lhs_grad_neumann_;
}

inline bool opheliePhiBoundaryGradNeumannPostprocess(const OphelieParameters &params)
{
    return params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann &&
           (params.phi_boundary_lhs_grad_neumann_ || params.phi_boundary_grad_neumann_projection_);
}

/** Precompute boundary mask / outward normal / g_n=-ω n·A on host (fixed during phi solve). */
inline size_t setupOpheliePhiBoundaryParticleFields(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                    const OphelieParameters &params,
                                                    const OpheliePhiBoundaryGeometryContext &geom, Real dp)
{
    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
    const Real omega = params.omega();

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);

    Real *mask = particles.getVariableDataByName<Real>(names.phi_boundary_mask);
    Vecd *normal = particles.getVariableDataByName<Vecd>(names.phi_boundary_normal);
    Real *g_n = particles.getVariableDataByName<Real>(names.phi_boundary_gn);

    size_t n_boundary = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dist = boundaryDistanceFromContext(pos[i], geom);
        if (dist <= boundary_width)
        {
            mask[i] = 1.0;
            normal[i] = boundaryOutwardNormalFromContext(pos[i], geom);
            g_n[i] = -omega * normal[i].dot(a_src[i]);
            ++n_boundary;
        }
        else
        {
            mask[i] = 0.0;
            normal[i] = Vecd::Zero();
            g_n[i] = 0.0;
        }
    }

    syncVariableToDevice<Real>(particles, names.phi_boundary_mask);
    syncVariableToDevice<Vecd>(particles, names.phi_boundary_normal);
    syncVariableToDevice<Real>(particles, names.phi_boundary_gn);
#if SPHINXSYS_USE_SYCL
    // Ensure device buffers exist and match host (finalizeLoadIn is a no-op until delegated).
    (void)particles.getVariableByName<Real>(names.phi_boundary_mask)->DelegatedData(execution::par_device);
    (void)particles.getVariableByName<Vecd>(names.phi_boundary_normal)->DelegatedData(execution::par_device);
    (void)particles.getVariableByName<Real>(names.phi_boundary_gn)->DelegatedData(execution::par_device);
    syncVariableToDevice<Real>(particles, names.phi_boundary_mask);
    syncVariableToDevice<Vecd>(particles, names.phi_boundary_normal);
    syncVariableToDevice<Real>(particles, names.phi_boundary_gn);
#endif
    return n_boundary;
}

/** Device: grad_phi_i ← grad_phi_i − n (n·grad_phi_i − g_n) on boundary shell particles. */
class ApplyOpheliePhiBoundaryGradNeumannProjectionCK : public LocalDynamics
{
  public:
    ApplyOpheliePhiBoundaryGradNeumannProjectionCK(SolidBody &sph_body, const OphelieGlassFieldNames &names)
        : LocalDynamics(sph_body),
          dv_boundary_mask_(particles_->template getVariableByName<Real>(names.phi_boundary_mask)),
          dv_boundary_normal_(particles_->template getVariableByName<Vecd>(names.phi_boundary_normal)),
          dv_boundary_gn_(particles_->template getVariableByName<Real>(names.phi_boundary_gn)),
          dv_grad_phi_(particles_->template getVariableByName<Vecd>(names.grad_phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : boundary_mask_(encloser.dv_boundary_mask_->DelegatedData(ex_policy)),
              boundary_normal_(encloser.dv_boundary_normal_->DelegatedData(ex_policy)),
              boundary_gn_(encloser.dv_boundary_gn_->DelegatedData(ex_policy)),
              grad_phi_(encloser.dv_grad_phi_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (boundary_mask_[index_i] < Real(0.5))
            {
                return;
            }
            const Vecd n_out = boundary_normal_[index_i];
            const Real target_gn = boundary_gn_[index_i];
            const Real grad_n = grad_phi_[index_i].dot(n_out);
            grad_phi_[index_i] -= n_out * (grad_n - target_gn);
        }

      protected:
        Real *boundary_mask_;
        Vecd *boundary_normal_;
        Real *boundary_gn_;
        Vecd *grad_phi_;
    };

  protected:
    DiscreteVariable<Real> *dv_boundary_mask_;
    DiscreteVariable<Vecd> *dv_boundary_normal_;
    DiscreteVariable<Real> *dv_boundary_gn_;
    DiscreteVariable<Vecd> *dv_grad_phi_;
};

/**
 * After grad_phi SPH estimate: enforce n·grad_phi = g_n = −ω n·A on boundary shell particles.
 * Caller must ensure Neumann mode / activation flags (see ProjectionDynamics wrapper).
 */
inline void applyOpheliePhiBoundaryGradNeumannProjection(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const OpheliePhiBoundaryGeometryContext &geom, Real dp)
{
    if (params.phi_boundary_mode_ != OpheliePhiBoundaryMode::OneSidedNeumann)
    {
        return;
    }

    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
    const Real omega = params.omega();

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);

    for (size_t i = 0; i < n; ++i)
    {
        if (boundaryDistanceFromContext(pos[i], geom) > boundary_width)
        {
            continue;
        }
        const Vecd n_out = boundaryOutwardNormalFromContext(pos[i], geom);
        const Real g_n = -omega * n_out.dot(a_src[i]);
        const Real grad_n = grad_phi[i].dot(n_out);
        grad_phi[i] -= n_out * (grad_n - g_n);
    }
    syncVariableToDevice<Vecd>(particles, names.grad_phi_imag);
}

/** Host projection using precomputed PhiBoundaryMask/Normal/Gn (device-safe after Interaction CK). */
inline void applyOpheliePhiBoundaryGradNeumannProjectionFromFields(BaseParticles &particles,
                                                                   const OphelieGlassFieldNames &names)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    // Mask/normal/gn are authored on host in setupOpheliePhiBoundaryParticleFields; do not
    // prepareForOutput them here (device copies may be stale and would overwrite host).
    const Real *mask = particles.getVariableDataByName<Real>(names.phi_boundary_mask);
    const Vecd *normal = particles.getVariableDataByName<Vecd>(names.phi_boundary_normal);
    const Real *g_n = particles.getVariableDataByName<Real>(names.phi_boundary_gn);
    Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);

    for (size_t i = 0; i < n; ++i)
    {
        if (mask[i] < Real(0.5))
        {
            continue;
        }
        const Vecd n_out = normal[i];
        const Real grad_n = grad_phi[i].dot(n_out);
        grad_phi[i] -= n_out * (grad_n - g_n[i]);
    }
    syncVariableToDevice<Vecd>(particles, names.grad_phi_imag);
}

template <class ExecutionPolicy>
inline void applyOpheliePhiBoundaryGradNeumannProjectionDynamics(
    SolidBody &glass_body, const OphelieGlassFieldNames &names, const OphelieParameters &params, bool for_lhs,
    const OpheliePhiBoundaryGeometryContext *geom = nullptr, Real dp = 0.0)
{
    const bool active =
        for_lhs ? opheliePhiBoundaryGradNeumannInLhs(params) : opheliePhiBoundaryGradNeumannPostprocess(params);
    if (!active)
    {
        return;
    }
    BaseParticles &particles = glass_body.getBaseParticles();
    if (geom != nullptr)
    {
        applyOpheliePhiBoundaryGradNeumannProjection(particles, names, params, *geom, dp);
        return;
    }
    applyOpheliePhiBoundaryGradNeumannProjectionFromFields(particles, names);
}

/**
 * One-sided Neumann RHS flux closure (first version, host-only):
 *   n·σ∇φ = −ω n·σA  →  g_n = −ω n·A
 *   rhs_i += (V_i / d_b) * σ_i * g_n
 * with d_b = max(dist_to_surface, 0.25 dp).
 */
inline OpheliePhiNeumannRhsCorrectionStats applyOpheliePhiOneSidedNeumannRhsCorrection(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const OpheliePhiBoundaryGeometryContext &geom, Real dp)
{
    OpheliePhiNeumannRhsCorrectionStats stats;
    if (params.phi_boundary_mode_ != OpheliePhiBoundaryMode::OneSidedNeumann ||
        params.phi_boundary_lhs_grad_neumann_)
    {
        return stats;
    }

    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
    const Real omega = params.omega();
    const Real d_b_min = Real(0.25) * dp;

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    Real correction_l2_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dist = boundaryDistanceFromContext(pos[i], geom);
        if (dist > boundary_width)
        {
            continue;
        }
        ++stats.n_boundary;
        const Vecd n_out = boundaryOutwardNormalFromContext(pos[i], geom);
        const Real d_b = std::max(dist, d_b_min);
        const Real g_n = -omega * n_out.dot(a_src[i]);
        const Real correction = (vol[i] / d_b) * sigma[i] * g_n;
        rhs[i] += correction;
        correction_l2_sq += vol[i] * correction * correction;
        stats.correction_max = std::max(stats.correction_max, std::abs(correction));
    }

    stats.correction_l2 = std::sqrt(correction_l2_sq);
    syncVariableToDevice<Real>(particles, names.phi_rhs_imag);
    return stats;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_BOUNDARY_H
