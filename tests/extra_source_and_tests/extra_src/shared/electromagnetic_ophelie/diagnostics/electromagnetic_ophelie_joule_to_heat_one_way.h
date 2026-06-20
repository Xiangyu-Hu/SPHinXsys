#ifndef ELECTROMAGNETIC_OPHELIE_JOULE_TO_HEAT_ONE_WAY_H
#define ELECTROMAGNETIC_OPHELIE_JOULE_TO_HEAT_ONE_WAY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "base_local_dynamics.h"
#include "simple_algorithms_ck.h"

#include <cmath>
#include <utility>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kOphelieTemperatureField = "Temperature";
/** Increment stored from zero so tiny ΔT is not lost when T ≈ 300 K (fp limit). */
inline constexpr const char *kOphelieThermalDeltaTField = "OphelieThermalDeltaT";
inline constexpr const char *kOphelieThermalLaplaceTField = "OphelieThermalLaplaceT";
inline constexpr const char *kOphelieThermalConductivityField = "OphelieThermalConductivity";
inline constexpr const char *kOphelieThermalBoundaryMaskField = "OphelieThermalBoundaryMask";

/** Molten-glass order-of-magnitude defaults for one-way thermal step (no feedback). */
struct OphelieJouleHeatOneWayMaterialProps
{
    Real rho = 2500.0;
    Real cp = 1200.0;
    /** Isotropic thermal conductivity [W/(m·K)]; Jacoutot molten-glass order-of-magnitude. */
    Real k = 1.0;
    Real t_initial = 300.0;
};

struct OphelieJouleHeatOneWayStepResult
{
    size_t n_steps = 0;
    Real max_per_particle_rel_err = 0.0;
    Real mean_delta_t = 0.0;
    Real max_delta_t = 0.0;
    Real total_joule_energy_j = 0.0;
    Real total_thermal_energy_j = 0.0;
    Real energy_balance_rel_err = 0.0;
    Real vol_weighted_delta_t = 0.0;
    Real vol_weighted_expected_delta_t = 0.0;
    /** Volume fraction where |ΔT - Q·dt/(ρ·cp)| exceeds tolerance (closure diagnostic). */
    Real closure_mismatch_vol_fraction = 0.0;
    Real closure_inline_energy_rel_err = 0.0;
    /** Stage 4.1: cold-wall Dirichlet compliance (1 = all boundary particles at T_wall). */
    Real boundary_dirichlet_compliance = 1.0;
    Real max_temperature = 0.0;
};

/** Particle Q field used as thermal source (complex edge-flux → recon complex Q on device). */
inline std::string ophelieJouleHeatSourceFieldForThermal(const OphelieGlassFieldNames &names,
                                                         const OphelieParameters &params)
{
    if (params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        return names.joule_heat_edge_recon_complex;
    }
    return names.joule_heat;
}

inline Real ophelieJouleHeatOneWayDeltaTExpected(Real q, Real dt, Real rho, Real cp)
{
    return q * dt / (rho * cp + TinyReal);
}

inline void hostOphelieJouleHeatOneWayVolWeightedDeltaT(BaseParticles &particles, const std::string &q_field,
                                                       const std::string &temperature_field, Real t_initial, Real dt,
                                                       Real rho, Real cp, size_t n, Real &delta_t_vol_weighted,
                                                       Real &expected_delta_t_vol_weighted,
                                                       bool host_temperature_authoritative = false)
{
    if (!host_temperature_authoritative)
    {
        syncVariableToHost<Real>(particles, temperature_field);
    }
    syncVariableToHost<Real>(particles, q_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real vol_sum = 0.0;
    delta_t_vol_weighted = 0.0;
    expected_delta_t_vol_weighted = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        vol_sum += vol[i];
        delta_t_vol_weighted += (temperature[i] - t_initial) * vol[i];
        expected_delta_t_vol_weighted += ophelieJouleHeatOneWayDeltaTExpected(q[i], dt, rho, cp) * vol[i];
    }
    delta_t_vol_weighted /= vol_sum + TinyReal;
    expected_delta_t_vol_weighted /= vol_sum + TinyReal;
}

/** One explicit Euler step: delta_T += Q*dt/(rho*cp); T = T0 + delta_T. */
inline void hostApplyOphelieJouleHeatOneWayTemperatureStep(BaseParticles &particles, const std::string &q_field,
                                                            const std::string &delta_t_field,
                                                            const std::string &temperature_field, Real t_initial,
                                                            Real dt, Real rho, Real cp, size_t n,
                                                            bool sync_q_from_device = true)
{
    if (sync_q_from_device)
    {
        syncVariableToHost<Real>(particles, q_field);
    }
    Real *delta_t = particles.getVariableDataByName<Real>(delta_t_field);
    Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real factor = dt / (rho * cp + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        delta_t[i] += q[i] * factor;
        temperature[i] = t_initial + delta_t[i];
    }
}

/** Copy edge-recon E/J/Q into primary fields so VTP / legacy thermal paths see JouleHeat. */
template <class ExecutionPolicy>
inline void syncOphelieJouleHeatPrimaryForThermalOneWay(SolidBody &glass_body, const OphelieGlassFieldNames &names,
                                                        const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, CopyOphelieEdgeReconToPrimaryEJQCK> copy_primary(
        glass_body, names, params.edge_flux_complex_ && ophelieUseEdgeFluxElectromotiveRhs(params));
    copy_primary.exec();
}

class ApplyOphelieJouleHeatOneWayTemperatureStepCK : public LocalDynamics
{
  public:
    ApplyOphelieJouleHeatOneWayTemperatureStepCK(SPHBody &sph_body, const std::string &q_field,
                                                 const std::string &delta_t_field,
                                                 const std::string &temperature_field, Real t_initial, Real dt,
                                                 Real rho, Real cp)
        : LocalDynamics(sph_body), t_initial_(t_initial), factor_(dt / (rho * cp + TinyReal)),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), factor_(encloser.factor_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            delta_t_[index_i] += q_[index_i] * factor_;
            temperature_[index_i] = t_initial_ + delta_t_[index_i];
        }

      protected:
        Real t_initial_;
        Real factor_;
        Real *q_;
        Real *delta_t_;
        Real *temperature_;
    };

  protected:
    Real t_initial_;
    Real dt_;
    Real factor_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
};

template <class ExecutionPolicy>
inline void execOphelieJouleHeatOneWayTemperatureSteps(SolidBody &glass_body, const std::string &q_field,
                                                       const std::string &delta_t_field,
                                                       const std::string &temperature_field, Real t_initial, Real dt,
                                                       Real rho, Real cp, size_t n_steps)
{
    StateDynamics<ExecutionPolicy, ApplyOphelieJouleHeatOneWayTemperatureStepCK> thermal_step(
        glass_body, q_field, delta_t_field, temperature_field, t_initial, dt, rho, cp);
    for (size_t step = 0; step < n_steps; ++step)
    {
        thermal_step.exec();
    }
}

/** Ensure device buffer exists without host→device upload (preserves EM-computed Q on device). */
template <class ExecutionPolicy, typename DataType>
inline void ensureOphelieVariableDelegatedOnDevice(BaseParticles &particles, const std::string &field_name)
{
#if SPHINXSYS_USE_SYCL
    (void)particles.template getVariableByName<DataType>(field_name)->DelegatedData(ExecutionPolicy{});
#else
    (void)particles;
    (void)field_name;
#endif
}

/** Upload host-initialized delta_T / Temperature; Q stays device-authoritative after EM. */
inline void syncOphelieJouleHeatOneWayThermalFieldsToDevice(BaseParticles &particles,
                                                            const std::string &delta_t_field,
                                                            const std::string &temperature_field)
{
    syncVariableToDevice<Real>(particles, delta_t_field);
    syncVariableToDevice<Real>(particles, temperature_field);
}

class OphelieJouleHeatOneWayJouleEnergyReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieJouleHeatOneWayJouleEnergyReduceCK(SPHBody &sph_body, const std::string &q_field, Real effective_dt)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), effective_dt_(effective_dt),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayJouleEnergy";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : effective_dt_(encloser.effective_dt_), q_(encloser.dv_q_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return q_[index_i] * vol_[index_i] * effective_dt_;
        }

      protected:
        Real effective_dt_;
        Real *q_;
        Real *vol_;
    };

  protected:
    Real effective_dt_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_vol_;
};

class OphelieJouleHeatOneWayThermalEnergyReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieJouleHeatOneWayThermalEnergyReduceCK(SPHBody &sph_body, const std::string &delta_t_field, Real rho, Real cp)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), rho_cp_(rho * cp),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayThermalEnergy";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : rho_cp_(encloser.rho_cp_), delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return rho_cp_ * delta_t_[index_i] * vol_[index_i];
        }

      protected:
        Real rho_cp_;
        Real *delta_t_;
        Real *vol_;
    };

  protected:
    Real rho_cp_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_vol_;
};

class OphelieJouleHeatOneWayMaxDeltaTReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    OphelieJouleHeatOneWayMaxDeltaTReduceCK(SPHBody &sph_body, const std::string &delta_t_field)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field))
    {
        quantity_name_ = "OphelieJouleHeatOneWayMaxDeltaT";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return delta_t_[index_i];
        }

      protected:
        Real *delta_t_;
    };

  protected:
    DiscreteVariable<Real> *dv_delta_t_;
};

class OphelieJouleHeatOneWaySumDeltaTReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieJouleHeatOneWaySumDeltaTReduceCK(SPHBody &sph_body, const std::string &delta_t_field)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field))
    {
        quantity_name_ = "OphelieJouleHeatOneWaySumDeltaT";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return delta_t_[index_i];
        }

      protected:
        Real *delta_t_;
    };

  protected:
    DiscreteVariable<Real> *dv_delta_t_;
};

class OphelieJouleHeatOneWayMaxRelErrReduceCK : public BaseLocalDynamicsReduce<ReduceParticleMax, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, UnsignedInt>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceParticleMax, SPHBody>;

    OphelieJouleHeatOneWayMaxRelErrReduceCK(SPHBody &sph_body, const std::string &q_field,
                                            const std::string &delta_t_field, Real effective_dt, Real rho, Real cp,
                                            Real q_threshold)
        : BaseDynamicsType(sph_body), effective_dt_(effective_dt),
          inv_rho_cp_(Real(1.0) / (rho * cp + TinyReal)), q_threshold_(q_threshold),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field))
    {
        quantity_name_ = "OphelieJouleHeatOneWayMaxRelErr";
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
        Real Result(const ReduceReturnType &reduced_value) { return reduced_value.first; }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : effective_dt_(encloser.effective_dt_), inv_rho_cp_(encloser.inv_rho_cp_),
              q_threshold_(encloser.q_threshold_), q_(encloser.dv_q_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (q_[index_i] <= q_threshold_)
            {
                return ReduceReturnType(0.0, index_i);
            }
            const Real expected_delta_t = q_[index_i] * effective_dt_ * inv_rho_cp_;
            if (expected_delta_t <= TinyReal)
            {
                return ReduceReturnType(0.0, index_i);
            }
            return ReduceReturnType(std::abs(delta_t_[index_i] - expected_delta_t) / expected_delta_t, index_i);
        }

      protected:
        Real effective_dt_;
        Real inv_rho_cp_;
        Real q_threshold_;
        Real *q_;
        Real *delta_t_;
    };

  protected:
    Real effective_dt_;
    Real inv_rho_cp_;
    Real q_threshold_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_delta_t_;
};

class OphelieJouleHeatOneWayVolWeightedDeltaReduceCK
    : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    OphelieJouleHeatOneWayVolWeightedDeltaReduceCK(SPHBody &sph_body, const std::string &delta_t_field)
        : BaseDynamicsType(sph_body),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayVolWeightedDeltaT";
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
            : delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return ReduceReturnType(delta_t_[index_i] * vol_[index_i], vol_[index_i]);
        }

      protected:
        Real *delta_t_;
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_vol_;
};

class OphelieJouleHeatOneWayVolWeightedExpectedDeltaReduceCK
    : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    OphelieJouleHeatOneWayVolWeightedExpectedDeltaReduceCK(SPHBody &sph_body, const std::string &q_field,
                                                           Real effective_dt, Real rho, Real cp)
        : BaseDynamicsType(sph_body), effective_dt_(effective_dt), inv_rho_cp_(Real(1.0) / (rho * cp + TinyReal)),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayVolWeightedExpectedDeltaT";
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
            : effective_dt_(encloser.effective_dt_), inv_rho_cp_(encloser.inv_rho_cp_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)), vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real expected_delta_t = q_[index_i] * effective_dt_ * inv_rho_cp_;
            return ReduceReturnType(expected_delta_t * vol_[index_i], vol_[index_i]);
        }

      protected:
        Real effective_dt_;
        Real inv_rho_cp_;
        Real *q_;
        Real *vol_;
    };

  protected:
    Real effective_dt_;
    Real inv_rho_cp_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_vol_;
};

class OphelieJouleHeatOneWayMismatchVolReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieJouleHeatOneWayMismatchVolReduceCK(SPHBody &sph_body, const std::string &q_field,
                                              const std::string &delta_t_field, Real effective_dt, Real rho, Real cp)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), effective_dt_(effective_dt),
          inv_rho_cp_(Real(1.0) / (rho * cp + TinyReal)),
          dv_q_(particles_->template getVariableByName<Real>(q_field)),
          dv_delta_t_(particles_->template getVariableByName<Real>(delta_t_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayMismatchVol";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : effective_dt_(encloser.effective_dt_), inv_rho_cp_(encloser.inv_rho_cp_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real expected_delta_t = q_[index_i] * effective_dt_ * inv_rho_cp_;
            if (std::abs(delta_t_[index_i] - expected_delta_t) <=
                std::max(Real(1.0e-12), Real(1.0e-6) * std::abs(expected_delta_t)))
            {
                return 0.0;
            }
            return vol_[index_i];
        }

      protected:
        Real effective_dt_;
        Real inv_rho_cp_;
        Real *q_;
        Real *delta_t_;
        Real *vol_;
    };

  protected:
    Real effective_dt_;
    Real inv_rho_cp_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_delta_t_;
    DiscreteVariable<Real> *dv_vol_;
};

class OphelieJouleHeatOneWayTotalVolReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OphelieJouleHeatOneWayTotalVolReduceCK(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieJouleHeatOneWayTotalVol";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return vol_[index_i];
        }

      protected:
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_vol_;
};

template <class ExecutionPolicy>
inline OphelieJouleHeatOneWayStepResult execOphelieJouleHeatOneWayClosureDiagnostics(
    SolidBody &glass_body, const std::string &q_field, const std::string &delta_t_field, Real effective_dt, Real rho,
    Real cp, size_t n, Real q_threshold)
{
    OphelieJouleHeatOneWayStepResult result;
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayJouleEnergyReduceCK> joule_energy(
        glass_body, q_field, effective_dt);
    result.total_joule_energy_j = joule_energy.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayThermalEnergyReduceCK> thermal_energy(
        glass_body, delta_t_field, rho, cp);
    result.total_thermal_energy_j = thermal_energy.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayMaxDeltaTReduceCK> max_delta_t(glass_body, delta_t_field);
    result.max_delta_t = max_delta_t.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWaySumDeltaTReduceCK> sum_delta_t(glass_body, delta_t_field);
    result.mean_delta_t = sum_delta_t.exec() / static_cast<Real>(n + TinyReal);
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayMaxRelErrReduceCK> max_rel_err(
        glass_body, q_field, delta_t_field, effective_dt, rho, cp, q_threshold);
    result.max_per_particle_rel_err = max_rel_err.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayVolWeightedDeltaReduceCK> vol_weighted_delta(
        glass_body, delta_t_field);
    result.vol_weighted_delta_t = vol_weighted_delta.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayVolWeightedExpectedDeltaReduceCK> vol_weighted_expected(
        glass_body, q_field, effective_dt, rho, cp);
    result.vol_weighted_expected_delta_t = vol_weighted_expected.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayMismatchVolReduceCK> mismatch_vol(
        glass_body, q_field, delta_t_field, effective_dt, rho, cp);
    const Real mismatch_vol_sum = mismatch_vol.exec();
    ReduceDynamicsCK<ExecutionPolicy, OphelieJouleHeatOneWayTotalVolReduceCK> total_vol(glass_body);
    result.closure_mismatch_vol_fraction = mismatch_vol_sum / (total_vol.exec() + TinyReal);
    result.closure_inline_energy_rel_err =
        std::abs(result.total_thermal_energy_j - result.total_joule_energy_j) /
        (std::abs(result.total_joule_energy_j) + TinyReal);
    result.energy_balance_rel_err = result.closure_inline_energy_rel_err;
    return result;
}

inline Real hostMaxJouleHeatOneWayTemperatureRelativeError(BaseParticles &particles, const std::string &q_field,
                                                           const std::string &temperature_field, Real t_initial,
                                                           Real dt, Real rho, Real cp, size_t n,
                                                           Real q_threshold = 0.0, bool host_temperature_authoritative = true)
{
    if (!host_temperature_authoritative)
    {
        syncVariableToHost<Real>(particles, temperature_field);
    }
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    Real max_rel_err = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        if (q[i] <= q_threshold)
        {
            continue;
        }
        const Real delta_t_expected = ophelieJouleHeatOneWayDeltaTExpected(q[i], dt, rho, cp);
        if (delta_t_expected <= TinyReal)
        {
            continue;
        }
        const Real delta_t = temperature[i] - t_initial;
        max_rel_err = std::max(max_rel_err, std::abs(delta_t - delta_t_expected) / delta_t_expected);
    }
    return max_rel_err;
}

inline void registerOphelieJouleHeatTemperatureField(BaseParticles &particles, Real t_initial)
{
    particles.registerStateVariable<Real>(kOphelieThermalDeltaTField, Real(0));
    particles.registerStateVariable<Real>(kOphelieTemperatureField, t_initial);
    const size_t n = particles.TotalRealParticles();
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    Real *temperature = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    for (size_t i = 0; i < n; ++i)
    {
        delta_t[i] = Real(0);
        temperature[i] = t_initial;
    }
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
}

inline void hostRegisterOphelieTemperatureField(BaseParticles &particles, Real t_initial)
{
    registerOphelieJouleHeatTemperatureField(particles, t_initial);
}

inline Real hostOphelieJouleHeatOneWayTotalThermalEnergy(BaseParticles &particles, const std::string &temperature_field,
                                                         Real t_initial, Real rho, Real cp, size_t n,
                                                         bool host_temperature_authoritative = true)
{
    if (!host_temperature_authoritative)
    {
        syncVariableToHost<Real>(particles, temperature_field);
    }
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real energy = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        energy += rho * cp * (temperature[i] - t_initial) * vol[i];
    }
    return energy;
}

inline Real hostOphelieJouleHeatOneWayTotalJouleEnergy(BaseParticles &particles, const std::string &q_field, Real dt,
                                                       size_t n, bool host_q_authoritative = false)
{
    if (!host_q_authoritative)
    {
        syncVariableToHost<Real>(particles, q_field);
    }
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real energy = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        energy += q[i] * vol[i] * dt;
    }
    return energy;
}

/** Volume-weighted energy closure; requires host-authoritative Q/T and synced VolumetricMeasure. */
inline void hostOphelieJouleHeatOneWayEnergyClosureFromSyncedFields(
    const Real *q, const Real *temperature, const Real *vol, Real t_initial, Real effective_dt, Real rho, Real cp,
    size_t n, OphelieJouleHeatOneWayStepResult &result)
{
    const Real inv_rho_cp = Real(1.0) / (rho * cp + TinyReal);
    Real vol_sum = 0.0;
    Real vol_weighted_delta_t = 0.0;
    Real vol_weighted_expected_delta_t = 0.0;
    Real mismatch_vol = 0.0;
    result.total_joule_energy_j = 0.0;
    result.total_thermal_energy_j = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        vol_sum += vol[i];
        const Real delta_t = temperature[i] - t_initial;
        const Real expected_delta_t = q[i] * effective_dt * inv_rho_cp;
        vol_weighted_delta_t += delta_t * vol[i];
        vol_weighted_expected_delta_t += expected_delta_t * vol[i];
        result.total_joule_energy_j += q[i] * vol[i] * effective_dt;
        result.total_thermal_energy_j += rho * cp * delta_t * vol[i];
        if (std::abs(delta_t - expected_delta_t) > std::max(Real(1.0e-12), Real(1.0e-6) * std::abs(expected_delta_t)))
        {
            mismatch_vol += vol[i];
        }
    }
    result.vol_weighted_delta_t = vol_weighted_delta_t / (vol_sum + TinyReal);
    result.vol_weighted_expected_delta_t = vol_weighted_expected_delta_t / (vol_sum + TinyReal);
    result.closure_mismatch_vol_fraction = mismatch_vol / (vol_sum + TinyReal);
    result.closure_inline_energy_rel_err =
        std::abs(result.total_thermal_energy_j - result.total_joule_energy_j) /
        (std::abs(result.total_joule_energy_j) + TinyReal);
    result.energy_balance_rel_err = result.closure_inline_energy_rel_err;
}

/** Frozen-source explicit Euler on device: T += Q*dt/(rho*cp); closure via device ReduceDynamicsCK. */
template <class ExecutionPolicy>
inline OphelieJouleHeatOneWayStepResult applyOphelieJouleHeatOneWayTemperatureSteps(
    SolidBody &glass_body, BaseParticles &particles, const std::string &q_field, const std::string &temperature_field,
    Real dt, const OphelieJouleHeatOneWayMaterialProps &material, size_t n, size_t n_steps,
    Real q_threshold = 1.0e-6)
{
    OphelieJouleHeatOneWayStepResult result;
    result.n_steps = n_steps;
    const Real effective_dt = dt * static_cast<Real>(n_steps);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, q_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, kOphelieThermalDeltaTField);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, temperature_field);
    ensureOphelieVariableDelegatedOnDevice<ExecutionPolicy, Real>(particles, "VolumetricMeasure");
    syncOphelieJouleHeatOneWayThermalFieldsToDevice(particles, kOphelieThermalDeltaTField, temperature_field);
    execOphelieJouleHeatOneWayTemperatureSteps<ExecutionPolicy>(
        glass_body, q_field, kOphelieThermalDeltaTField, temperature_field, material.t_initial, dt, material.rho,
        material.cp, n_steps);
    OphelieJouleHeatOneWayStepResult diagnostics = execOphelieJouleHeatOneWayClosureDiagnostics<ExecutionPolicy>(
        glass_body, q_field, kOphelieThermalDeltaTField, effective_dt, material.rho, material.cp, n, q_threshold);
    result.max_per_particle_rel_err = diagnostics.max_per_particle_rel_err;
    result.mean_delta_t = diagnostics.mean_delta_t;
    result.max_delta_t = diagnostics.max_delta_t;
    result.total_joule_energy_j = diagnostics.total_joule_energy_j;
    result.total_thermal_energy_j = diagnostics.total_thermal_energy_j;
    result.energy_balance_rel_err = diagnostics.energy_balance_rel_err;
    result.vol_weighted_delta_t = diagnostics.vol_weighted_delta_t;
    result.vol_weighted_expected_delta_t = diagnostics.vol_weighted_expected_delta_t;
    result.closure_mismatch_vol_fraction = diagnostics.closure_mismatch_vol_fraction;
    result.closure_inline_energy_rel_err = diagnostics.closure_inline_energy_rel_err;
    return result;
}

struct OphelieFrenchEmJouleHeatOneWayResult
{
    OphelieFrenchEmSolveResult em;
    OphelieJouleHeatOneWayStepResult thermal;
    Real joule_power_w = 0.0;
    Real phi_eq_res_vol = 0.0;
};

/** French reduced EM solve then frozen-Q thermal one-way (no EM feedback). */
template <class ExecutionPolicy>
inline OphelieFrenchEmJouleHeatOneWayResult runFrenchReducedEmThenJouleHeatOneWay(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french, Real thermal_dt, size_t thermal_steps,
    const OphelieJouleHeatOneWayMaterialProps &material = OphelieJouleHeatOneWayMaterialProps())
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

    const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
    result.thermal = applyOphelieJouleHeatOneWayTemperatureSteps<ExecutionPolicy>(
        glass_body, particles, q_field, kOphelieTemperatureField, thermal_dt, material, n, thermal_steps);
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_JOULE_TO_HEAT_ONE_WAY_H
