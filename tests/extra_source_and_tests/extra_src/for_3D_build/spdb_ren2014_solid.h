#ifndef SPDB_REN2014_SOLID_H
#define SPDB_REN2014_SOLID_H

#include "elastic_solid.h"
#include "spdb_arrhenius.h"

namespace SPH
{
/**
 * Neo-Hookean carrier plus Ren 2014 J2 rate.
 * The inelastic step is ConstituteKernel::ElasticLeftCauchy, the same hook
 * StructureIntegration1stHalf uses for HardeningPlasticSolid.
 * Youngs modulus is a caller-supplied demo_assumption, not a Ren 2014 value.
 * |F-I| or |J-1| beyond the small-strain gate skips the inelastic update.
 */
class Ren2014X20Cr13Solid : public NeoHookeanSolid
{
  public:
    Ren2014X20Cr13Solid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real temperature_K,
                        bool inelastic_enabled = true);
    virtual ~Ren2014X20Cr13Solid() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0);

    class ConstituteKernel : public NeoHookeanSolid::ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, Ren2014X20Cr13Solid &encloser);
        Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0);

      protected:
        Real temperature_K_;
        bool inelastic_enabled_;
        Matd *inverse_plastic_strain_;
        Real *inelastic_eq_;
        Real *range_flag_;
    };

    DiscreteVariable<Matd> *dv_inverse_plastic_strain_;
    DiscreteVariable<Real> *dv_inelastic_eq_;
    DiscreteVariable<Real> *dv_range_flag_;

  protected:
    Real temperature_K_;
    bool inelastic_enabled_;
    Matd *inverse_plastic_strain_;
    Real *inelastic_eq_;
    Real *range_flag_;
};

struct Ren2014PhysicalCommit
{
    size_t particle_count = 0;
    size_t updated = 0;
    size_t rejected = 0;
    size_t large_strain = 0;
    Real mean_sigma_trial_Pa = 0.0;
    Real mean_sigma_relaxed_Pa = 0.0;
    Real mean_inelastic = 0.0;
    Real max_inelastic = 0.0;
    Real mean_rate_per_s = 0.0;
    int accepted = 0;
};

inline Matd ren2014_cauchy_stress(const Matd &F, const Matd &inverse_plastic, Real G0, Real K0, Real *sigma_eq, bool *small_strain)
{
    Matd be = F * inverse_plastic * F.transpose();
    Matd normalized_be = be * math::pow(be.determinant(), -Real(OneOverDimensions));
    const Real J = F.determinant();
    const Real deformation_norm = (F - Matd::Identity()).norm();
    if (small_strain != nullptr)
        *small_strain = deformation_norm <= 0.05 && ABS(J - Real(1.0)) <= Real(0.02) && J > 0.0;
    const Real iso = normalized_be.trace() * OneOverDimensions;
    Matd dev_tau = G0 * (normalized_be - iso * Matd::Identity());
    const Real vol = 0.5 * K0 * (J * J - 1.0);
    Matd tau = dev_tau + vol * Matd::Identity();
    if (sigma_eq != nullptr)
        *sigma_eq = J > 0.0 ? math::sqrt(1.5) * dev_tau.norm() / J : 0.0;
    if (!(J > 0.0))
        return Matd::Zero();
    tau /= J;
    return tau;
}

/** One physical increment from the accepted particle history.
 *  Acoustic relaxation must keep Ren2014X20Cr13Solid inelastic_enabled false,
 *  so this Δt is not added again on later pseudo-time steps.
 *  A rejected particle restores every particle to the arrays from entry.
 */
Ren2014PhysicalCommit commit_ren2014_physical_increment(BaseParticles &particles, Real shear_modulus,
                                                        Real temperature_K, Real dt_physical,
                                                        Real rate_scale = 1.0);

inline Matd relax_ren2014_elastic_left_cauchy(const Matd &F, size_t index_i, Real dt, Real G0, Real temperature_K,
                                              bool inelastic_enabled, Matd *inverse_plastic_strain, Real *inelastic_eq,
                                              Real *range_flag, Real rate_scale = 1.0)
{
    Matd be = F * inverse_plastic_strain[index_i] * F.transpose();
    Real be_det = be.determinant();
    Matd normalized_be = be * math::pow(be_det, -Real(OneOverDimensions));
    const Real J = F.determinant();
    const Real deformation_norm = (F - Matd::Identity()).norm();
    const bool small_strain = deformation_norm <= 0.05 && ABS(J - Real(1.0)) <= Real(0.02);
    const Real dt_eff = dt * rate_scale;
    if (!inelastic_enabled || !small_strain || !(dt_eff > 0.0))
    {
        if (!small_strain)
            range_flag[index_i] = 1.0;
        return normalized_be;
    }

    const Real iso = normalized_be.trace() * OneOverDimensions;
    Matd dev_tau = G0 * (normalized_be - iso * Matd::Identity());
    const Real norm = dev_tau.norm();
    if (!(norm > 0.0) || !(J > 0.0))
        return normalized_be;

    const Real sigma_trial = math::sqrt(1.5) * norm / J;
    // rate_scale folds the demo 850–900 C gate into an effective physical Δt.
    const spdb::EquivalentReturn relaxed = spdb::equivalent_stress_backward_euler(
        sigma_trial, G0, dt_eff, temperature_K, inelastic_eq[index_i], spdb::StartupPolicy::EndpointExtend);
    if (relaxed.accepted != 1)
    {
        range_flag[index_i] = 1.0;
        return normalized_be;
    }

    const Real norm_new = relaxed.sigma_eq_Pa * J / math::sqrt(1.5);
    dev_tau *= norm_new / norm;
    Matd relaxed_be = dev_tau / G0 + iso * Matd::Identity();
    normalized_be = relaxed_be * math::pow(relaxed_be.determinant(), -Real(OneOverDimensions));
    Matd inverse_F = F.inverse();
    inverse_plastic_strain[index_i] = inverse_F * normalized_be * inverse_F.transpose();
    inelastic_eq[index_i] = relaxed.inelastic_eq;
    return normalized_be;
}

} // namespace SPH

#include "spdb_ren2014_solid.hpp"

#endif
