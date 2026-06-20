#ifndef ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H
#define ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_thermal_vtp.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Stage 4.1: isotropic conduction + optional cold-wall Dirichlet (no EM feedback). */
struct OphelieThermalDiffusionOneWayOptions
{
    bool enable_diffusion = false;
    bool enable_cold_wall_dirichlet = false;
    /** Shell width = factor * dp (match phi boundary convention). */
    Real boundary_width_factor = 1.5;
    Real pair_weight_regularization = Real(0.01);
};

struct OphelieThermalDiffusionOneWayStepResult : public OphelieJouleHeatOneWayStepResult
{
};

inline void registerOphelieThermalDiffusionAuxFields(BaseParticles &particles, Real thermal_conductivity)
{
    particles.registerStateVariable<Real>(kOphelieThermalLaplaceTField, Real(0));
    particles.registerStateVariable<Real>(kOphelieThermalConductivityField, thermal_conductivity);
    particles.registerStateVariable<Real>(kOphelieThermalBoundaryMaskField, Real(0));
    const size_t n = particles.TotalRealParticles();
    Real *conductivity = particles.getVariableDataByName<Real>(kOphelieThermalConductivityField);
    for (size_t i = 0; i < n; ++i)
    {
        conductivity[i] = thermal_conductivity;
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalLaplaceTField);
    syncVariableToDevice<Real>(particles, kOphelieThermalConductivityField);
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
}

/** Host: mark boundary shell particles for Dirichlet T = T_wall (French cylinder or analytic box). */
inline size_t setupOphelieThermalDirichletBoundaryMask(BaseParticles &particles,
                                                       const OphelieThermalDiffusionOneWayOptions &options,
                                                       const OpheliePhiBoundaryGeometryContext &geom, Real dp)
{
    if (!options.enable_cold_wall_dirichlet)
    {
        return 0;
    }
    const size_t n = particles.TotalRealParticles();
    const Real boundary_width = options.boundary_width_factor * dp;
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *mask = particles.getVariableDataByName<Real>(kOphelieThermalBoundaryMaskField);
    size_t n_boundary = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dist = boundaryDistanceFromContext(pos[i], geom);
        if (dist <= boundary_width)
        {
            mask[i] = 1.0;
            ++n_boundary;
        }
        else
        {
            mask[i] = 0.0;
        }
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalBoundaryMaskField);
    return n_boundary;
}

class ApplyOphelieJouleHeatDiffusionCombinedStepCK : public LocalDynamics
{
  public:
    ApplyOphelieJouleHeatDiffusionCombinedStepCK(SPHBody &sph_body, const std::string &q_field,
                                                 const std::string &delta_t_field,
                                                 const std::string &temperature_field,
                                                 const std::string &laplace_t_field, Real t_initial, Real dt,
                                                 Real rho, Real cp, Real k, bool apply_diffusion)
        : LocalDynamics(sph_body), t_initial_(t_initial), joule_factor_(dt / (rho * cp + TinyReal)),
          diffusion_factor_(apply_diffusion ? -dt * k / (rho * cp + TinyReal) : Real(0)),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_laplace_t_(particles_->template getVariableByName<Real>(laplace_t_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), joule_factor_(encloser.joule_factor_),
              diffusion_factor_(encloser.diffusion_factor_), q_(encloser.dv_q_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              laplace_t_(encloser.dv_laplace_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            delta_t_[index_i] += joule_factor_ * q_[index_i] + diffusion_factor_ * laplace_t_[index_i];
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real joule_factor_;
        Real diffusion_factor_;
        Real *q_;
        Real *delta_t_;
        Real *temperature_;
        Real *laplace_t_;
    };

  protected:
    Real t_initial_;
    Real joule_factor_;
    Real diffusion_factor_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_laplace_t_;
};

class ApplyOphelieThermalDirichletWallCK : public LocalDynamics
{
  public:
    ApplyOphelieThermalDirichletWallCK(SPHBody &sph_body, const std::string &delta_t_field,
                                       const std::string &temperature_field, Real t_wall)
        : LocalDynamics(sph_body), t_wall_(t_wall),
          dv_boundary_mask_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryMaskField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_wall_(encloser.t_wall_), boundary_mask_(encloser.dv_boundary_mask_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (boundary_mask_[index_i] < Real(0.5))
            {
                return;
            }
            delta_t_[index_i] = Real(0);
            temperature_[index_i] = t_wall_;
        }

      protected:
        Real t_wall_;
        Real *boundary_mask_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_wall_;
    DiscreteVariable<Real> *dv_boundary_mask_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieThermalBoundaryComplianceReduceCK
    : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    OphelieThermalBoundaryComplianceReduceCK(SPHBody &sph_body, const std::string &temperature_field, Real t_wall,
                                             Real tolerance)
        : BaseDynamicsType(sph_body), t_wall_(t_wall), tolerance_(tolerance),
          dv_boundary_mask_(particles_->template getVariableByName<Real>(kOphelieThermalBoundaryMaskField)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
        quantity_name_ = "OphelieThermalBoundaryCompliance";
    }

    class FinishDynamics
    {
      public:
        using OutputType = Real;
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser)
        {
            (void)encloser;
        }
        Real Result(const ReduceReturnType &reduced_value)
        {
            return reduced_value.first / (reduced_value.second + TinyReal);
        }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_wall_(encloser.t_wall_), tolerance_(encloser.tolerance_),
              boundary_mask_(encloser.dv_boundary_mask_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (boundary_mask_[index_i] < Real(0.5))
            {
                return ReduceReturnType(0.0, 0.0);
            }
            const Real compliant = std::abs(temperature_[index_i] - t_wall_) <= tolerance_ ? 1.0 : 0.0;
            return ReduceReturnType(compliant, 1.0);
        }

      protected:
        Real t_wall_;
        Real tolerance_;
        Real *boundary_mask_;
        Real *temperature_;
    };

  protected:
    Real t_wall_;
    Real tolerance_;
    DiscreteVariable<Real> *dv_boundary_mask_;
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieThermalMaxTemperatureReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    OphelieThermalMaxTemperatureReduceCK(SPHBody &sph_body, const std::string &temperature_field)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
        quantity_name_ = "OphelieThermalMaxTemperature";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return temperature_[index_i];
        }

      protected:
        Real *temperature_;
    };

  protected:
    DiscreteVariable<Real> *dv_temperature_;
};

template <class ExecutionPolicy>
inline void execOphelieJouleHeatDiffusionOneWayStep(SolidBody &glass_body, Inner<> &glass_inner,
                                                    const std::string &q_field, const std::string &temperature_field,
                                                    Real t_initial, Real dt, const OphelieJouleHeatOneWayMaterialProps &material,
                                                    const OphelieThermalDiffusionOneWayOptions &options)
{
    if (options.enable_diffusion)
    {
        UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
        UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> laplace_temperature(
            glass_inner, temperature_field, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
            options.pair_weight_regularization);
        laplace_temperature.exec();
    }
    else
    {
        StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_laplace(glass_body, kOphelieThermalLaplaceTField);
        zero_laplace.exec();
    }

    StateDynamics<ExecutionPolicy, ApplyOphelieJouleHeatDiffusionCombinedStepCK> combined_step(
        glass_body, q_field, kOphelieThermalDeltaTField, temperature_field, kOphelieThermalLaplaceTField, t_initial, dt,
        material.rho, material.cp, material.k, options.enable_diffusion);
    combined_step.exec();

    if (options.enable_cold_wall_dirichlet)
    {
        StateDynamics<ExecutionPolicy, ApplyOphelieThermalDirichletWallCK> dirichlet_wall(
            glass_body, kOphelieThermalDeltaTField, temperature_field, t_initial);
        dirichlet_wall.exec();
    }
}

template <class ExecutionPolicy>
inline OphelieThermalDiffusionOneWayStepResult execOphelieThermalDiffusionDiagnostics(
    SolidBody &glass_body, const std::string &q_field, Real effective_dt, Real rho, Real cp, size_t n,
    Real t_wall, Real q_threshold, const OphelieThermalDiffusionOneWayOptions &options)
{
    OphelieThermalDiffusionOneWayStepResult result;
    const OphelieJouleHeatOneWayStepResult base =
        execOphelieJouleHeatOneWayClosureDiagnostics<ExecutionPolicy>(
            glass_body, q_field, kOphelieThermalDeltaTField, effective_dt, rho, cp, n, q_threshold);
    result.n_steps = base.n_steps;
    result.max_per_particle_rel_err = base.max_per_particle_rel_err;
    result.mean_delta_t = base.mean_delta_t;
    result.max_delta_t = base.max_delta_t;
    result.total_joule_energy_j = base.total_joule_energy_j;
    result.total_thermal_energy_j = base.total_thermal_energy_j;
    result.energy_balance_rel_err = base.energy_balance_rel_err;
    result.vol_weighted_delta_t = base.vol_weighted_delta_t;
    result.vol_weighted_expected_delta_t = base.vol_weighted_expected_delta_t;
    result.closure_mismatch_vol_fraction = base.closure_mismatch_vol_fraction;
    result.closure_inline_energy_rel_err = base.closure_inline_energy_rel_err;
    if (options.enable_cold_wall_dirichlet)
    {
        ReduceDynamicsCK<ExecutionPolicy, OphelieThermalBoundaryComplianceReduceCK> boundary_compliance(
            glass_body, kOphelieTemperatureField, t_wall, std::max(Real(1.0e-6), Real(1.0e-4) * std::abs(t_wall)));
        result.boundary_dirichlet_compliance = boundary_compliance.exec();
    }
    else
    {
        result.boundary_dirichlet_compliance = 1.0;
    }
    ReduceDynamicsCK<ExecutionPolicy, OphelieThermalMaxTemperatureReduceCK> max_temperature(
        glass_body, kOphelieTemperatureField);
    result.max_temperature = max_temperature.exec();
    return result;
}

inline void copyOphelieJouleHeatOneWayStepResult(OphelieJouleHeatOneWayStepResult &dst,
                                                 const OphelieJouleHeatOneWayStepResult &src)
{
    dst.n_steps = src.n_steps;
    dst.max_per_particle_rel_err = src.max_per_particle_rel_err;
    dst.mean_delta_t = src.mean_delta_t;
    dst.max_delta_t = src.max_delta_t;
    dst.total_joule_energy_j = src.total_joule_energy_j;
    dst.total_thermal_energy_j = src.total_thermal_energy_j;
    dst.energy_balance_rel_err = src.energy_balance_rel_err;
    dst.vol_weighted_delta_t = src.vol_weighted_delta_t;
    dst.vol_weighted_expected_delta_t = src.vol_weighted_expected_delta_t;
    dst.closure_mismatch_vol_fraction = src.closure_mismatch_vol_fraction;
    dst.closure_inline_energy_rel_err = src.closure_inline_energy_rel_err;
    dst.boundary_dirichlet_compliance = src.boundary_dirichlet_compliance;
    dst.max_temperature = src.max_temperature;
}

template <class ExecutionPolicy>
inline OphelieThermalDiffusionOneWayStepResult applyOphelieJouleHeatDiffusionOneWaySteps(
    SolidBody &glass_body, Inner<> &glass_inner, BaseParticles &particles, const std::string &q_field,
    const std::string &temperature_field, Real dt, const OphelieJouleHeatOneWayMaterialProps &material, size_t n,
    size_t n_steps, const OphelieThermalDiffusionOneWayOptions &options, Real q_threshold = 1.0e-6,
    const OphelieThermalVtpRecordingOptions *recording = nullptr)
{
    OphelieThermalDiffusionOneWayStepResult result;
    result.n_steps = n_steps;
    const Real effective_dt = dt * static_cast<Real>(n_steps);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, q_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalDeltaTField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, temperature_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalLaplaceTField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalConductivityField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalBoundaryMaskField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, "VolumetricMeasure");
    syncOphelieJouleHeatOneWayThermalFieldsToDevice(particles, kOphelieThermalDeltaTField, temperature_field);

    for (size_t step = 0; step < n_steps; ++step)
    {
        execOphelieJouleHeatDiffusionOneWayStep<ExecutionPolicy>(glass_body, glass_inner, q_field, temperature_field,
                                                                 material.t_initial, dt, material, options);
        writeOphelieThermalVtpIfDue(recording, step, n_steps);
    }

    result = execOphelieThermalDiffusionDiagnostics<ExecutionPolicy>(
        glass_body, q_field, effective_dt, material.rho, material.cp, n, material.t_initial, q_threshold, options);
    return result;
}

/** French EM + Joule heat with optional isotropic diffusion and cold-crucible Dirichlet shell. */
template <class ExecutionPolicy>
inline OphelieFrenchEmJouleHeatOneWayResult runFrenchReducedEmThenJouleHeatDiffusionOneWay(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french, Real thermal_dt, size_t thermal_steps,
    const OphelieJouleHeatOneWayMaterialProps &material, const OphelieThermalDiffusionOneWayOptions &thermal_options,
    const OphelieThermalVtpRecordingOptions *recording = nullptr)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    OphelieFrenchEmJouleHeatOneWayResult result;
    result.em = runFrenchReducedEmPipeline<ExecutionPolicy>(glass_body, glass_inner, names, params, french);
    result.joule_power_w = result.em.joule_power_raw;
    result.phi_eq_res_vol = result.em.phi_eq_res_vol;

    if (ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        syncOphelieJouleHeatPrimaryForThermalOneWay<ExecutionPolicy>(glass_body, names, params);
    }

    const size_t n = particles.TotalRealParticles();
    registerOphelieJouleHeatTemperatureField(particles, material.t_initial);
    registerOphelieThermalDiffusionAuxFields(particles, material.k);

    if (thermal_options.enable_cold_wall_dirichlet)
    {
        OpheliePhiBoundaryGeometryContext geom;
        geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
        geom.french = french;
        (void)setupOphelieThermalDirichletBoundaryMask(particles, thermal_options, geom, french.dp);
    }

    const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
    const OphelieThermalDiffusionOneWayStepResult thermal_diffusion =
        applyOphelieJouleHeatDiffusionOneWaySteps<ExecutionPolicy>(
            glass_body, glass_inner, particles, q_field, kOphelieTemperatureField, thermal_dt, material, n,
            thermal_steps, thermal_options, 1.0e-6, recording);
    copyOphelieJouleHeatOneWayStepResult(result.thermal, thermal_diffusion);
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_THERMAL_DIFFUSION_ONE_WAY_H
