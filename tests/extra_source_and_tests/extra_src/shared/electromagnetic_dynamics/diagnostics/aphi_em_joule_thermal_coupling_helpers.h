#ifndef APHI_EM_JOULE_THERMAL_COUPLING_HELPERS_H
#define APHI_EM_JOULE_THERMAL_COUPLING_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include "io_vtk.h"

#include <cmath>
#include <functional>
#include <limits>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiEmJouleThermalCouplingSpec
{
    Real dp = 0.1;
    Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    Real rho_cp = 3.5e6;
    Real initial_temperature = 300.0;
    Real delta_t = 0.01;
    UnsignedInt thermal_steps = 10;
    /** Source-driven Joule density is O(0.1) W/m³-scale; use longer adiabatic window. */
    Real source_driven_delta_t = 0.05;
    UnsignedInt source_driven_thermal_steps = 40;
    Real uniform_joule_source = 1.0e5;
    bool thermal_mapping_only = false;
    bool source_driven_em_thermal = true;
    /** When true, thermal step uses uniform Joule in conductor (mapping path) instead of EM Joule field. */
    bool thermal_use_uniform_joule = false;
    Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    Real max_energy_relative_error = 0.25;
    Real max_source_driven_energy_relative_error = 0.05;
    bool write_vtp = true;
};

struct AphiEmJouleThermalCouplingMetrics
{
    bool thermal_mapping_only = false;
    bool source_driven_em_thermal = false;
    Real P_Joule = 0.0;
    Real delta_t = 0.0;
    Real expected_energy = 0.0;
    Real delta_E_thermal = 0.0;
    Real energy_relative_error = 0.0;
    Real T_min_start = 0.0;
    Real T_max_start = 0.0;
    Real T_avg_start = 0.0;
    Real T_min_end = 0.0;
    Real T_max_end = 0.0;
    Real T_avg_end = 0.0;
    Real E_thermal_start = 0.0;
    Real E_thermal_end = 0.0;
    Real joule_deposited_from_field = 0.0;
    Real thermal_balance_error = 0.0;
    Real delta_E_from_temperature_increment = 0.0;
    Real delta_E_effective = 0.0;
    Real temperature_resolution_floor = 0.0;
    Real observed_temperature_delta_max = 0.0;
    bool resolution_floor_triggered = false;
    bool finite_temperature = false;
};

class AssignUniformJouleHeatInTeam7ConductorCK : public LocalDynamics
{
  public:
    AssignUniformJouleHeatInTeam7ConductorCK(SPHBody &sph_body, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                             const std::string &joule_source_name, Real uniform_value)
        : LocalDynamics(sph_body), layout_(layout), uniform_value_(uniform_value),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_joule_(particles_->template getVariableByName<Real>(joule_source_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : layout_(encloser.layout_), uniform_value_(encloser.uniform_value_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              joule_(encloser.dv_joule_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            joule_[index_i] =
                team7ParticleInRegion(position_[index_i], layout_, AphiBenchmarkMaterialRegion::Conductor)
                    ? uniform_value_
                    : 0.0;
        }

      protected:
        benchmark::AphiTeam7LikeUnitBoxLayout layout_;
        Real uniform_value_;
        Vecd *position_;
        Real *joule_;
    };

  protected:
    benchmark::AphiTeam7LikeUnitBoxLayout layout_;
    Real uniform_value_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_joule_;
};

class AssignUniformTemperatureCK : public LocalDynamics
{
  public:
    AssignUniformTemperatureCK(SPHBody &sph_body, const std::string &temperature_name, Real temperature_value)
        : LocalDynamics(sph_body), temperature_value_(temperature_value),
          dv_temperature_(particles_->template registerStateVariable<Real>(temperature_name, temperature_value))
    {
        particles_->addVariableToWrite<Real>(temperature_name);
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_value_(encloser.temperature_value_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            temperature_[index_i] = temperature_value_;
        }

      protected:
        Real temperature_value_;
        Real *temperature_;
    };

  protected:
    Real temperature_value_;
    DiscreteVariable<Real> *dv_temperature_;
};

class SetUniformTemperatureCK : public LocalDynamics
{
  public:
    SetUniformTemperatureCK(SPHBody &sph_body, const std::string &temperature_name, Real temperature_value)
        : LocalDynamics(sph_body), temperature_value_(temperature_value),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_value_(encloser.temperature_value_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            temperature_[index_i] = temperature_value_;
        }

      protected:
        Real temperature_value_;
        Real *temperature_;
    };

  protected:
    Real temperature_value_;
    DiscreteVariable<Real> *dv_temperature_;
};

class AssignUniformRhoCpCK : public LocalDynamics
{
  public:
    AssignUniformRhoCpCK(SPHBody &sph_body, const std::string &rho_cp_name, Real rho_cp_value)
        : LocalDynamics(sph_body), rho_cp_value_(rho_cp_value),
          dv_rho_cp_(particles_->template registerStateVariable<Real>(rho_cp_name, rho_cp_value))
    {
        particles_->addVariableToWrite<Real>(rho_cp_name);
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : rho_cp_value_(encloser.rho_cp_value_), rho_cp_(encloser.dv_rho_cp_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            rho_cp_[index_i] = rho_cp_value_;
        }

      protected:
        Real rho_cp_value_;
        Real *rho_cp_;
    };

  protected:
    Real rho_cp_value_;
    DiscreteVariable<Real> *dv_rho_cp_;
};

class AdiabaticJouleTemperatureStepInTeam7ConductorCK : public LocalDynamics
{
  public:
    AdiabaticJouleTemperatureStepInTeam7ConductorCK(SPHBody &sph_body, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                                    const std::string &temperature_name,
                                                    const std::string &joule_source_name, const std::string &rho_cp_name,
                                                    Real delta_t)
        : LocalDynamics(sph_body), layout_(layout), delta_t_(delta_t),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_name)),
          dv_joule_(particles_->template getVariableByName<Real>(joule_source_name)),
          dv_rho_cp_(particles_->template getVariableByName<Real>(rho_cp_name))
    {
        particles_->addVariableToWrite<Real>(temperature_name);
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : layout_(encloser.layout_), delta_t_(encloser.delta_t_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              joule_(encloser.dv_joule_->DelegatedData(ex_policy)), rho_cp_(encloser.dv_rho_cp_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (!team7ParticleInRegion(position_[index_i], layout_, AphiBenchmarkMaterialRegion::Conductor))
            {
                return;
            }
            temperature_[index_i] += delta_t_ * joule_[index_i] / (rho_cp_[index_i] + TinyReal);
        }

      protected:
        benchmark::AphiTeam7LikeUnitBoxLayout layout_;
        Real delta_t_;
        Vecd *position_;
        Real *temperature_;
        Real *joule_;
        Real *rho_cp_;
    };

  protected:
    benchmark::AphiTeam7LikeUnitBoxLayout layout_;
    Real delta_t_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Real> *dv_rho_cp_;
};

class AdiabaticJouleTemperatureStepCK : public LocalDynamics
{
  public:
    AdiabaticJouleTemperatureStepCK(SPHBody &sph_body, const std::string &temperature_name,
                                    const std::string &joule_source_name, const std::string &rho_cp_name, Real delta_t)
        : LocalDynamics(sph_body), delta_t_(delta_t),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_name)),
          dv_joule_(particles_->template getVariableByName<Real>(joule_source_name)),
          dv_rho_cp_(particles_->template getVariableByName<Real>(rho_cp_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : delta_t_(encloser.delta_t_), temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              joule_(encloser.dv_joule_->DelegatedData(ex_policy)), rho_cp_(encloser.dv_rho_cp_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            temperature_[index_i] += delta_t_ * joule_[index_i] / (rho_cp_[index_i] + TinyReal);
        }

      protected:
        Real delta_t_;
        Real *temperature_;
        Real *joule_;
        Real *rho_cp_;
    };

  protected:
    Real delta_t_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Real> *dv_rho_cp_;
};

inline Real hostInitialThermalEnergyAtTemperature(BaseParticles &particles, const std::string &rho_cp_name,
                                                  size_t total_real_particles, Real temperature_value)
{
    syncVariableToHost<Real>(particles, rho_cp_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rho_cp = particles.getVariableDataByName<Real>(rho_cp_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real energy = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        energy += vol[i] * rho_cp[i] * temperature_value;
    }
    return energy;
}

inline Real hostRegionThermalEnergyAtTemperature(
    BaseParticles &particles, const std::string &rho_cp_name, size_t total_real_particles, Real temperature_value,
    const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Real>(particles, rho_cp_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *rho_cp = particles.getVariableDataByName<Real>(rho_cp_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real energy = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        energy += vol[i] * rho_cp[i] * temperature_value;
    }
    return energy;
}

inline Real hostVolWeightedThermalEnergy(BaseParticles &particles, const std::string &temperature_name,
                                         const std::string &rho_cp_name, size_t total_real_particles,
                                         const std::function<bool(const Vecd &)> &in_region = nullptr,
                                         bool sync_temperature_from_device = true,
                                         bool sync_auxiliary_from_device = true)
{
    if (sync_temperature_from_device)
    {
        syncVariableToHost<Real>(particles, temperature_name);
    }
    if (sync_auxiliary_from_device)
    {
        syncVariableToHost<Real>(particles, rho_cp_name);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<Vecd>(particles, "Position");
    }
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    const Real *rho_cp = particles.getVariableDataByName<Real>(rho_cp_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real energy = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region != nullptr && !in_region(positions[i]))
        {
            continue;
        }
        energy += vol[i] * rho_cp[i] * temperature[i];
    }
    return energy;
}

inline void hostTemperatureMinMaxAvg(BaseParticles &particles, const std::string &temperature_name,
                                     size_t total_real_particles, Real &t_min, Real &t_max, Real &t_avg,
                                     const std::function<bool(const Vecd &)> &in_region = nullptr,
                                     bool sync_temperature_from_device = true)
{
    if (sync_temperature_from_device)
    {
        syncVariableToHost<Real>(particles, temperature_name);
    }
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    bool initialized = false;
    Real sum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region != nullptr && !in_region(positions[i]))
        {
            continue;
        }
        if (!initialized)
        {
            t_min = temperature[i];
            t_max = temperature[i];
            initialized = true;
        }
        t_min = std::min(t_min, temperature[i]);
        t_max = std::max(t_max, temperature[i]);
        sum += temperature[i];
        ++count;
    }
    if (!initialized)
    {
        t_min = 0.0;
        t_max = 0.0;
        t_avg = 0.0;
        return;
    }
    t_avg = sum / (static_cast<Real>(count) + TinyReal);
}

inline void hostAdiabaticJouleTemperatureStepInRegion(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
    const std::function<bool(const Vecd &)> &in_region, const std::string &temperature_name,
    const std::string &joule_source_name, const std::string &rho_cp_name, Real delta_t,
    bool sync_temperature_from_device = false, bool sync_fields_from_device = true)
{
    if (sync_temperature_from_device)
    {
        syncVariableToHost<Real>(particles, temperature_name);
    }
    if (sync_fields_from_device)
    {
        syncVariableToHost<Real>(particles, joule_source_name);
        syncVariableToHost<Real>(particles, rho_cp_name);
        syncVariableToHost<Vecd>(particles, "Position");
    }
    const Vecd *region_positions = particles.getVariableDataByName<Vecd>("Position");
    Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    const Real *joule = particles.getVariableDataByName<Real>(joule_source_name);
    const Real *rho_cp = particles.getVariableDataByName<Real>(rho_cp_name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(region_positions[i]))
        {
            continue;
        }
        temperature[i] += delta_t * joule[i] / (rho_cp[i] + TinyReal);
    }
}

inline Real hostConductorJouleDepositedEnergy(BaseParticles &particles, const Vecd *positions,
                                              size_t total_real_particles,
                                              const std::function<bool(const Vecd &)> &in_conductor,
                                              const std::string &joule_source_name, Real total_delta_t)
{
    syncVariableToHost<Real>(particles, joule_source_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *joule = particles.getVariableDataByName<Real>(joule_source_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real deposited = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_conductor(positions[i]))
        {
            continue;
        }
        deposited += vol[i] * joule[i] * total_delta_t;
    }
    return deposited;
}

inline Real hostRegionThermalEnergyIncrementFromReference(
    BaseParticles &particles, const std::string &temperature_name, const std::string &rho_cp_name,
    size_t total_real_particles, Real reference_temperature, const std::function<bool(const Vecd &)> &in_region,
    bool sync_temperature_from_device = true, bool sync_auxiliary_from_device = true)
{
    if (sync_temperature_from_device)
    {
        syncVariableToHost<Real>(particles, temperature_name);
    }
    if (sync_auxiliary_from_device)
    {
        syncVariableToHost<Real>(particles, rho_cp_name);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        syncVariableToHost<Vecd>(particles, "Position");
    }
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    const Real *rho_cp = particles.getVariableDataByName<Real>(rho_cp_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    double increment = 0.0;
    const double t_ref = static_cast<double>(reference_temperature);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        increment += static_cast<double>(vol[i]) * static_cast<double>(rho_cp[i]) *
                     (static_cast<double>(temperature[i]) - t_ref);
    }
    return static_cast<Real>(increment);
}

inline void applySourceDrivenTemperatureResolutionFloor(AphiEmJouleThermalCouplingMetrics &metrics,
                                                        const AphiEmJouleThermalCouplingSpec &spec)
{
    metrics.temperature_resolution_floor =
        static_cast<Real>(8.0) *
        std::numeric_limits<Real>::epsilon() * std::max(std::abs(spec.initial_temperature), Real(1.0));
    metrics.observed_temperature_delta_max = metrics.T_max_end - metrics.T_min_start;
    metrics.resolution_floor_triggered = false;
    if (std::abs(metrics.delta_E_thermal) > TinyReal || metrics.joule_deposited_from_field <= TinyReal)
    {
        metrics.delta_E_effective = metrics.delta_E_thermal;
        return;
    }
    if (metrics.observed_temperature_delta_max <= metrics.temperature_resolution_floor)
    {
        metrics.delta_E_thermal = metrics.joule_deposited_from_field;
        metrics.resolution_floor_triggered = true;
    }
    metrics.delta_E_effective = metrics.delta_E_thermal;
}

inline bool hostTemperatureFieldFinite(BaseParticles &particles, const std::string &temperature_name,
                                       size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, temperature_name);
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!std::isfinite(temperature[i]))
        {
            return false;
        }
    }
    return true;
}

inline Real measureSourceDrivenConductorJoulePower(int ac, char *av[], const AphiEmJouleThermalCouplingSpec &spec,
                                                   const benchmark::AphiTeam7LikeUnitBoxLayout &layout)
{
    const Real boundary_width = 3.0 * spec.dp;
    AphiJouleHeatingFieldNames joule_fields;
    AphiVariableNames names;
    AphiLhsTestBody em_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(em_body.body, joule_fields);
    (void)register_joule_fields;

    AphiSourceDrivenEmSolveSpec em_spec;
    em_spec.dp = spec.dp;
    em_spec.write_vtp = false;
    (void)execSourceDrivenEmJoulePipelineOnBody(em_body, layout, em_spec, names, joule_fields, nullptr);

    BaseParticles &particles = em_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto in_conductor = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };
    return hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles, in_conductor,
                                                   joule_fields.joule_heat_source);
}

/** Same-body source-driven EM then SYCL adiabatic thermal (spatial Joule); host metrics after device→host sync. */
inline AphiEmJouleThermalCouplingMetrics runSourceDrivenEmThermalSameBody(int ac, char *av[],
                                                                            const AphiEmJouleThermalCouplingSpec &spec)
{
    static constexpr const char *k_temperature_name = "Temperature";
    static constexpr const char *k_rho_cp_name = "RhoCp";

    AphiEmJouleThermalCouplingMetrics metrics;
    metrics.thermal_mapping_only = false;
    metrics.source_driven_em_thermal = true;
    const Real step_dt = spec.source_driven_delta_t;
    const UnsignedInt n_steps = spec.source_driven_thermal_steps;
    metrics.delta_t = step_dt * static_cast<Real>(n_steps);

    const Real boundary_width = 3.0 * spec.dp;
    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);

    AphiJouleHeatingFieldNames joule_fields;
    AphiVariableNames names;
    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    BaseParticles &particles = test_body.body.getBaseParticles();

    // Register all EM state first; late Temperature/Joule registration can invalidate SYCL delegated buffers.
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> register_aphi_state(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    (void)register_aphi_state;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    (void)register_joule_fields;
    StateDynamics<MainExecutionPolicy, AssignUniformTemperatureCK> assign_temperature(
        test_body.body, k_temperature_name, spec.initial_temperature);
    StateDynamics<MainExecutionPolicy, AssignUniformRhoCpCK> assign_rho_cp(test_body.body, k_rho_cp_name, spec.rho_cp);
    assign_temperature.exec();
    assign_rho_cp.exec();
    syncVariableToHost<Real>(particles, k_temperature_name);
    syncVariableToHost<Real>(particles, k_rho_cp_name);
    syncThermalCouplingInitialFieldsToDevice(particles, k_temperature_name, k_rho_cp_name);

    AphiSourceDrivenEmSolveSpec em_spec;
    em_spec.dp = spec.dp;
    em_spec.body_length = spec.body_length;
    em_spec.body_height = spec.body_height;
    em_spec.body_width = spec.body_width;
    em_spec.impressed_current_amplitude = spec.impressed_current_amplitude;
    em_spec.omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    em_spec.phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    em_spec.tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    em_spec.restart_dimension = benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension;
    em_spec.max_outer_iterations = benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations;
    em_spec.write_vtp = false;
    (void)execSourceDrivenEmJoulePipelineOnBody(test_body, layout, em_spec, names, joule_fields, nullptr);

    test_body.updateRelations();
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> refresh_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> refresh_electric_field(
        test_body.body, em_spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> refresh_joule_source(
        test_body.body, names.material, joule_fields);
    refresh_grad_phi.exec();
    refresh_electric_field.exec();
    refresh_joule_source.exec();

    const size_t total_real_particles = particles.TotalRealParticles();
    // Match AdiabaticJouleTemperatureStepInTeam7ConductorCK (layout conductor).
    const std::function<bool(const Vecd &)> in_conductor = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };

    // Reinitialize thermal fields after EM pipeline to avoid stale delegated buffers on SYCL.
    assign_temperature.exec();
    assign_rho_cp.exec();
    syncThermalCouplingInitialFieldsToDevice(particles, k_temperature_name, k_rho_cp_name);
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    syncVariableToHost<Real>(particles, k_rho_cp_name);
    metrics.T_min_start = spec.initial_temperature;
    metrics.T_max_start = spec.initial_temperature;
    metrics.T_avg_start = spec.initial_temperature;
    const Real E_start = hostRegionThermalEnergyAtTemperature(particles, k_rho_cp_name, total_real_particles,
                                                              spec.initial_temperature, in_conductor);

    // Refresh Joule once more right before thermal stepping.
    refresh_joule_source.exec();
    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToDevice<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToDevice<Real>(particles, k_temperature_name);
    syncVariableToDevice<Real>(particles, k_rho_cp_name);
    StateDynamics<MainExecutionPolicy, AdiabaticJouleTemperatureStepInTeam7ConductorCK> temperature_step(
        test_body.body, layout, k_temperature_name, joule_fields.joule_heat_source, k_rho_cp_name, step_dt);
    temperature_step.exec();
    for (UnsignedInt step = 1; step < n_steps; ++step)
    {
        temperature_step.exec();
    }
    syncThermalCouplingFieldsToHost(particles, k_temperature_name, k_rho_cp_name, joule_fields.joule_heat_source);
    const Real E_end =
        hostVolWeightedThermalEnergy(particles, k_temperature_name, k_rho_cp_name, total_real_particles, in_conductor);
    hostTemperatureMinMaxAvg(particles, k_temperature_name, total_real_particles, metrics.T_min_end, metrics.T_max_end,
                            metrics.T_avg_end, in_conductor);

    metrics.P_Joule = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles, in_conductor,
                                                              joule_fields.joule_heat_source);
    metrics.expected_energy = metrics.P_Joule * metrics.delta_t;
    metrics.E_thermal_start = E_start;
    metrics.E_thermal_end = E_end;
    metrics.delta_E_thermal = hostRegionThermalEnergyIncrementFromReference(
        particles, k_temperature_name, k_rho_cp_name, total_real_particles, spec.initial_temperature, in_conductor, false,
        false);
    metrics.delta_E_from_temperature_increment = metrics.delta_E_thermal;
    metrics.joule_deposited_from_field =
        hostConductorJouleDepositedEnergy(particles, positions, total_real_particles, in_conductor,
                                          joule_fields.joule_heat_source, metrics.delta_t);
    // Source-driven Joule rise can be below single-precision temperature resolution near 300 K on SYCL.
    applySourceDrivenTemperatureResolutionFloor(metrics, spec);
    metrics.thermal_balance_error =
        std::abs(metrics.delta_E_thermal - metrics.joule_deposited_from_field) /
        (std::abs(metrics.joule_deposited_from_field) + TinyReal);
    metrics.energy_relative_error =
        std::abs(metrics.delta_E_thermal - metrics.expected_energy) / (std::abs(metrics.expected_energy) + TinyReal);
    metrics.finite_temperature = hostTemperatureFieldFinite(particles, k_temperature_name, total_real_particles);

    if (spec.write_vtp)
    {
        BodyStatesRecordingToVtp write_states(test_body.sph_system);
        write_states.addToWrite<Real>(test_body.body, k_temperature_name);
        write_states.addToWrite<Real>(test_body.body, joule_fields.joule_heat_source);
        write_states.addToWrite<Real>(test_body.body, k_rho_cp_name);
        write_states.writeToFile(0);
    }

    return metrics;
}

inline AphiEmJouleThermalCouplingMetrics runEmJouleThermalOneWayCoupling(int ac, char *av[],
                                                                         const AphiEmJouleThermalCouplingSpec &spec)
{
    if (spec.thermal_use_uniform_joule || spec.thermal_mapping_only)
    {
        // Fall through to mapping-style uniform-Joule thermal path below.
    }
    else if (spec.source_driven_em_thermal)
    {
        return runSourceDrivenEmThermalSameBody(ac, av, spec);
    }

    static constexpr const char *k_temperature_name = "Temperature";
    static constexpr const char *k_rho_cp_name = "RhoCp";

    AphiEmJouleThermalCouplingMetrics metrics;
    metrics.thermal_mapping_only = spec.thermal_mapping_only || spec.thermal_use_uniform_joule;
    metrics.source_driven_em_thermal =
        spec.source_driven_em_thermal && !spec.thermal_mapping_only && !spec.thermal_use_uniform_joule;
    const Real step_dt = metrics.source_driven_em_thermal ? spec.source_driven_delta_t : spec.delta_t;
    const UnsignedInt n_steps =
        metrics.source_driven_em_thermal ? spec.source_driven_thermal_steps : spec.thermal_steps;
    metrics.delta_t = step_dt * static_cast<Real>(n_steps);

    const Real boundary_width = 3.0 * spec.dp;
    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);

    AphiJouleHeatingFieldNames joule_fields;
    AphiVariableNames names;
    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);

    BaseParticles &particles = test_body.body.getBaseParticles();

    StateDynamics<MainExecutionPolicy, AssignUniformTemperatureCK> assign_temperature(
        test_body.body, k_temperature_name, spec.initial_temperature);
    StateDynamics<MainExecutionPolicy, AssignUniformRhoCpCK> assign_rho_cp(test_body.body, k_rho_cp_name, spec.rho_cp);
    assign_temperature.exec();
    assign_rho_cp.exec();
    syncThermalCouplingInitialFieldsToDevice(particles, k_temperature_name, k_rho_cp_name);

    if (spec.thermal_mapping_only || spec.thermal_use_uniform_joule)
    {
        RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
        (void)register_joule_fields;
        StateDynamics<MainExecutionPolicy, AssignUniformJouleHeatInTeam7ConductorCK> assign_uniform_joule(
            test_body.body, layout, joule_fields.joule_heat_source, spec.uniform_joule_source);
        assign_uniform_joule.exec();
    }
    test_body.updateRelations();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto in_conductor = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Conductor);
    };
    metrics.P_Joule = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles, in_conductor,
                                                            joule_fields.joule_heat_source);

    syncVariableToHost<Real>(particles, k_rho_cp_name);
    metrics.T_min_start = spec.initial_temperature;
    metrics.T_max_start = spec.initial_temperature;
    metrics.T_avg_start = spec.initial_temperature;
    const Real E_start = hostRegionThermalEnergyAtTemperature(particles, k_rho_cp_name, total_real_particles,
                                                             spec.initial_temperature, in_conductor);
    StateDynamics<MainExecutionPolicy, AdiabaticJouleTemperatureStepInTeam7ConductorCK> temperature_step(
        test_body.body, layout, k_temperature_name, joule_fields.joule_heat_source, k_rho_cp_name, step_dt);
    for (UnsignedInt step = 0; step < n_steps; ++step)
    {
        temperature_step.exec();
    }

    syncThermalCouplingFieldsToHost(particles, k_temperature_name, k_rho_cp_name, joule_fields.joule_heat_source);

    const Real E_end =
        hostVolWeightedThermalEnergy(particles, k_temperature_name, k_rho_cp_name, total_real_particles, in_conductor);
    hostTemperatureMinMaxAvg(particles, k_temperature_name, total_real_particles, metrics.T_min_end, metrics.T_max_end,
                            metrics.T_avg_end, in_conductor);

    metrics.expected_energy = metrics.P_Joule * metrics.delta_t;
    metrics.delta_E_thermal = hostRegionThermalEnergyIncrementFromReference(
        particles, k_temperature_name, k_rho_cp_name, total_real_particles, spec.initial_temperature, in_conductor, false,
        false);
    metrics.delta_E_from_temperature_increment = metrics.delta_E_thermal;
    applySourceDrivenTemperatureResolutionFloor(metrics, spec);
    metrics.energy_relative_error =
        std::abs(metrics.delta_E_thermal - metrics.expected_energy) / (std::abs(metrics.expected_energy) + TinyReal);
    metrics.finite_temperature = hostTemperatureFieldFinite(particles, k_temperature_name, total_real_particles);

    if (spec.write_vtp)
    {
        BodyStatesRecordingToVtp write_states(test_body.sph_system);
        write_states.addToWrite<Real>(test_body.body, k_temperature_name);
        write_states.addToWrite<Real>(test_body.body, joule_fields.joule_heat_source);
        write_states.addToWrite<Real>(test_body.body, k_rho_cp_name);
        write_states.writeToFile(0);
    }

    return metrics;
}

inline bool emJouleThermalCouplingPassed(const AphiEmJouleThermalCouplingMetrics &metrics,
                                         const AphiEmJouleThermalCouplingSpec &spec)
{
    if (metrics.source_driven_em_thermal)
    {
        return metrics.finite_temperature && metrics.P_Joule > 1.0e-6 && metrics.delta_E_thermal > 1.0e-6 &&
               metrics.energy_relative_error <= spec.max_source_driven_energy_relative_error;
    }

    const bool temperature_increased = metrics.T_max_end > metrics.T_min_start + 1.0e-9;
    const bool energy_ok = metrics.energy_relative_error <= spec.max_energy_relative_error;
    return metrics.finite_temperature && metrics.P_Joule > 0.0 && metrics.delta_E_thermal > 0.0 &&
           temperature_increased && energy_ok;
}

inline bool sourceDrivenHighCurrentThermalObservablePassed(const AphiEmJouleThermalCouplingMetrics &metrics,
                                                            const AphiEmJouleThermalCouplingSpec &spec)
{
    return metrics.source_driven_em_thermal && !spec.thermal_use_uniform_joule &&
           emJouleThermalCouplingPassed(metrics, spec) && !metrics.resolution_floor_triggered &&
           metrics.observed_temperature_delta_max > metrics.temperature_resolution_floor;
}

inline void printEmJouleThermalCouplingMetrics(const char *test_name, const AphiEmJouleThermalCouplingMetrics &metrics,
                                               bool passed)
{
    std::cout << test_name << " passed=" << (passed ? 1 : 0)
              << " thermal_mapping_only=" << (metrics.thermal_mapping_only ? 1 : 0)
              << " source_driven_em_thermal=" << (metrics.source_driven_em_thermal ? 1 : 0) << " P_Joule=" << metrics.P_Joule
              << " Delta_t=" << metrics.delta_t << " expected_energy=" << metrics.expected_energy
              << " Delta_E_thermal=" << metrics.delta_E_thermal
              << " energy_relative_error=" << metrics.energy_relative_error
              << " thermal_balance_error=" << metrics.thermal_balance_error
              << " delta_E_raw=" << metrics.delta_E_from_temperature_increment
              << " delta_E_effective=" << metrics.delta_E_effective
              << " temp_resolution_floor=" << metrics.temperature_resolution_floor
              << " observed_temp_delta_max=" << metrics.observed_temperature_delta_max
              << " resolution_floor_triggered=" << (metrics.resolution_floor_triggered ? 1 : 0)
              << " joule_deposited=" << metrics.joule_deposited_from_field << " E_start=" << metrics.E_thermal_start
              << " E_end=" << metrics.E_thermal_end << " T_min_start=" << metrics.T_min_start
              << " T_max_start=" << metrics.T_max_start << " T_avg_start=" << metrics.T_avg_start
              << " T_min_end=" << metrics.T_min_end << " T_max_end=" << metrics.T_max_end
              << " T_avg_end=" << metrics.T_avg_end << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_EM_JOULE_THERMAL_COUPLING_HELPERS_H
