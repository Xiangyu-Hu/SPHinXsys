#ifndef ELECTROMAGNETIC_OPHELIE_SIGMA_T_COUPLING_H
#define ELECTROMAGNETIC_OPHELIE_SIGMA_T_COUPLING_H

#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_source.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

class ZeroOphelieNamedScalarFieldCK : public LocalDynamics
{
  public:
    explicit ZeroOphelieNamedScalarFieldCK(SPHBody &sph_body, const std::string &field_name)
        : LocalDynamics(sph_body), dv_field_(particles_->template getVariableByName<Real>(field_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : field_(encloser.dv_field_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            field_[index_i] = Real(0);
        }

      protected:
        Real *field_;
    };

  protected:
    DiscreteVariable<Real> *dv_field_;
};

struct OphelieSigmaTCouplingOptions
{
    size_t max_iterations = 8;
    Real sigma_relaxation = 0.3;
    Real relative_tol = 0.05;
    Real thermal_dt = 0.1;
    size_t thermal_steps_per_em = 1;
    /** Radial T seed so sigma is spatially non-uniform before heating (machinery gate). */
    Real temperature_seed_delta_k = 100.0;
};

struct OphelieSigmaTCouplingResult
{
    OphelieFrenchEmJouleHeatOneWayResult last_handoff;
    size_t coupling_iterations = 0;
    bool converged = false;
    Real final_sigma_rel = 0.0;
    Real final_power_rel = 0.0;
    Real sigma_min = 0.0;
    Real sigma_max = 0.0;
    Real sigma_mean = 0.0;
    Real joule_power_w = 0.0;
    bool paper_digitized_sigma_law = false;
};

/**
 * Minimal Stage 3.2 loop:
 *   seed/register T → σ(T) → EM(+optional A_glass) → Q → T
 * with under-relaxed sigma updates.
 */
template <class ExecutionPolicy>
inline OphelieSigmaTCouplingResult runFrenchReducedSigmaTEmThermalCoupling(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french, const OphelieJouleHeatOneWayMaterialProps &material,
    const OphelieTemperatureLaw &sigma_law, const OphelieSigmaTCouplingOptions &options,
    bool paper_digitized_sigma_law = false)
{
    OphelieSigmaTCouplingResult result;
    result.paper_digitized_sigma_law = paper_digitized_sigma_law;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    registerOphelieJouleHeatTemperatureField(particles, material.t_initial);
    registerOphelieFrenchTemperatureMaterialFields(particles, params.sigma_glass_);

    if (options.temperature_seed_delta_k > TinyReal)
    {
        StateDynamics<ExecutionPolicy, SeedOphelieRadialTemperatureCK> seed_temperature(
            glass_body, french.glass_center, french.glass_radius, material.t_initial, options.temperature_seed_delta_k);
        seed_temperature.exec();
    }

    StateDynamics<ExecutionPolicy, UpdateOphelieGlassSigmaFromTemperatureCK> update_sigma(
        glass_body, names, sigma_law, Real(1.0));
    update_sigma.exec();

    OphelieProgressLogger progress("sigma_t_coupling");
    StdVec<Real> previous_sigma;
    hostStoreScalarField(particles, names.sigma, previous_sigma, n);
    Real previous_power = 0.0;

    StateDynamics<ExecutionPolicy, UpdateOphelieGlassSigmaFromTemperatureCK> update_sigma_relaxed(
        glass_body, names, sigma_law, options.sigma_relaxation);
    StateDynamics<ExecutionPolicy, ZeroOphelieNamedScalarFieldCK> zero_heating_delta(glass_body,
                                                                                      kOphelieThermalDeltaTField);

    for (size_t iter = 0; iter < options.max_iterations; ++iter)
    {
        result.coupling_iterations = iter + 1;
        if (iter > 0)
        {
            update_sigma_relaxed.exec();
        }

        result.last_handoff = OphelieFrenchEmJouleHeatOneWayResult{};
        runFrenchReducedEmOrSelfInductionForThermalHandoff<ExecutionPolicy>(glass_body, glass_inner, names, params,
                                                                            french, result.last_handoff);

        // Keep Temperature history; reset heating-only DeltaT so interval closure stays valid.
        zero_heating_delta.exec();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
        result.last_handoff.thermal = applyOphelieJouleHeatOneWayTemperatureSteps<ExecutionPolicy>(
            glass_body, particles, q_field, kOphelieTemperatureField, options.thermal_dt, material, n,
            options.thermal_steps_per_em, 1.0e-6, true);

        result.joule_power_w = result.last_handoff.joule_power_w;
        hostOphelieSigmaFieldStats(particles, names.sigma, n, result.sigma_min, result.sigma_max, result.sigma_mean);
        result.final_sigma_rel = hostOphelieRelativeScalarFieldChange(particles, names.sigma, previous_sigma, n);
        result.final_power_rel =
            std::abs(result.joule_power_w - previous_power) / (std::abs(result.joule_power_w) + TinyReal);

        progress.log("outer iter " + std::to_string(result.coupling_iterations) + " P_joule_W=" +
                     std::to_string(result.joule_power_w) + " sigma_mean=" + std::to_string(result.sigma_mean) +
                     " sigma_min=" + std::to_string(result.sigma_min) + " sigma_max=" + std::to_string(result.sigma_max) +
                     " sigma_rel=" + std::to_string(result.final_sigma_rel) +
                     " power_rel=" + std::to_string(result.final_power_rel) +
                     " phi_eq_res_vol=" + std::to_string(result.last_handoff.phi_eq_res_vol));

        hostStoreScalarField(particles, names.sigma, previous_sigma, n);
        previous_power = result.joule_power_w;

        if (iter > 0 && result.final_sigma_rel < options.relative_tol && result.final_power_rel < options.relative_tol)
        {
            result.converged = true;
            break;
        }
    }

    progress.finish("converged=" + std::to_string(result.converged ? 1 : 0) +
                    " iters=" + std::to_string(result.coupling_iterations) +
                    " sigma_rel=" + std::to_string(result.final_sigma_rel) +
                    " power_rel=" + std::to_string(result.final_power_rel));
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_SIGMA_T_COUPLING_H
