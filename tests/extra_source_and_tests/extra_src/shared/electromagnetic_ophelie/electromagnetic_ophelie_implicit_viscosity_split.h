/**
 * @file electromagnetic_ophelie_implicit_viscosity_split.h
 * @brief Doc175/177: InteractionSplit + DampingPairwise implicit viscosity with local μ.
 *
 * Classical ParallelPolicy path (CPU split). Not a SYCL/CK kernel.
 *
 * Doc177: DampingRate must match ViscousForceCK semi-discrete RHS as dt→0.
 *
 * Per InteractionSplit::exec(dt):
 *  - particle_for_split visits every particle twice (forward + backward color sweeps),
 *    each call with dt_arg = dt/2;
 *  - DampingPairwise::interaction itself does forward+backward neighborhood sweeps.
 * Net from one particle's visits alone ≈ one outer dt of pairwise updates; both fluid
 * endpoints of an Inner pair are visited → another ×2. ViscousForceCK Inner
 * (NoKernelCorrection) contributes geometric ×2 via vec_r·((C_i+C_j)e)=2r.
 * Empirically (Doc177 sine Inner-only operator-diag, factor=1 → |a_imp|/|a_exp|≈2):
 *   Inner DampingRate = 0.5 * μ_harm.
 *
 * Wall contact updates only the fluid particle (no endpoint double-count). CK Wall has
 * geometric+mirror ×4. With InteractionSplit's dual color visit ≈×2 on fluid:
 *   Wall DampingRate = 2 * μ_i.
 *
 * Doc175 primary defect was empty classic cell_data_lists_ (Damping no-op).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_IMPLICIT_VISCOSITY_SPLIT_H
#define ELECTROMAGNETIC_OPHELIE_IMPLICIT_VISCOSITY_SPLIT_H

#include "base_particles.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "particle_dynamics_dissipation.h"

#include <string>

namespace SPH
{
using electromagnetics::ophelie::syncVariableToDevice;
using electromagnetics::ophelie::syncVariableToHost;

/** Match ViscousForceCK Inner after InteractionSplit dual-visit + pairwise endpoints. */
inline constexpr Real kOphelieCkMatchInnerDampingFactor = Real(0.5);
/** Match ViscousForceCK Contact<Wall> after InteractionSplit dual-visit. */
inline constexpr Real kOphelieCkMatchWallDampingFactor = Real(2);

/**
 * Inner fluid–fluid: harmonic mean μ_ij × CK-match factor.
 * Capacity = mass (c=1), matching FixedDampingRate default.
 */
class ParticleViscosityHarmonicDampingRate
{
  public:
    explicit ParticleViscosityHarmonicDampingRate(
        BaseParticles *particles,
        const std::string &mu_name = electromagnetics::ophelie::kOphelieDynamicViscosityField)
        : mu_(particles->getVariableDataByName<Real>(mu_name)),
          mass_(particles->getVariableDataByName<Real>("Mass")) {};

    Real DampingRate(size_t i, size_t j)
    {
        const Real mi = mu_[i];
        const Real mj = mu_[j];
        const Real mu_harm = Real(2) * mi * mj / (mi + mj + TinyReal);
        return kOphelieCkMatchInnerDampingFactor * mu_harm;
    }
    Real DampingRate(size_t i) { return kOphelieCkMatchInnerDampingFactor * mu_[i]; }
    Real Capacity(size_t i) { return mass_[i]; }

  protected:
    Real *mu_;
    Real *mass_;
};

/**
 * Fluid–wall: one-sided μ_i × CK-match wall factor.
 * Never index mu_[j] — j is a wall-body local index.
 */
class ParticleViscosityOneSideDampingRate
{
  public:
    explicit ParticleViscosityOneSideDampingRate(
        BaseParticles *particles,
        const std::string &mu_name = electromagnetics::ophelie::kOphelieDynamicViscosityField)
        : mu_(particles->getVariableDataByName<Real>(mu_name)),
          mass_(particles->getVariableDataByName<Real>("Mass")) {};

    Real DampingRate(size_t i, size_t /*j*/) { return kOphelieCkMatchWallDampingFactor * mu_[i]; }
    Real DampingRate(size_t i) { return kOphelieCkMatchWallDampingFactor * mu_[i]; }
    Real Capacity(size_t i) { return mass_[i]; }

  protected:
    Real *mu_;
    Real *mass_;
};

/** Remove ViscousForceCK contribution from ForcePrior and zero ViscousForce / PreviousViscousForce. */
inline void hostClearViscousForcePriorContribution(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "ForcePrior");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    syncVariableToHost<Vecd>(particles, "PreviousViscousForce");
    Vecd *fp = particles.getVariableDataByName<Vecd>("ForcePrior");
    Vecd *vf = particles.getVariableDataByName<Vecd>("ViscousForce");
    Vecd *pv = particles.getVariableDataByName<Vecd>("PreviousViscousForce");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        // ForcePrior += current - previous  ⇒ undo last viscous by subtracting previous (and any current).
        fp[i] -= pv[i];
        vf[i] = Vecd::Zero();
        pv[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, "ForcePrior");
    syncVariableToDevice<Vecd>(particles, "ViscousForce");
    syncVariableToDevice<Vecd>(particles, "PreviousViscousForce");
}

} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_IMPLICIT_VISCOSITY_SPLIT_H
