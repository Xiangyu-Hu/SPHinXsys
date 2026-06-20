/**
 * @file rh200_fake_joule_heat.h
 * @brief Step 1 fake JouleHeat source + energy budget helpers for RH200 glass stirring case.
 */
#ifndef RH200_FAKE_JOULE_HEAT_H
#define RH200_FAKE_JOULE_HEAT_H

#include "base_local_dynamics.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace rh200
{

inline constexpr const char *kJouleHeatField = "JouleHeat";
inline constexpr const char *kEnergyBudgetCsv = "./output/rh200_energy_budget.csv";

enum class FakeJouleMode
{
    Off,
    Uniform,
    AnalyticPosition,
    EmFixed,
    EmGrid,
};

inline bool jouleHeatEnabled(FakeJouleMode mode)
{
    return mode != FakeJouleMode::Off;
}

inline FakeJouleMode parseFakeJouleMode(const std::string &text)
{
    if (text == "uniform")
    {
        return FakeJouleMode::Uniform;
    }
    if (text == "analytic-position" || text == "analytic_position")
    {
        return FakeJouleMode::AnalyticPosition;
    }
    if (text == "em-fixed" || text == "em_fixed" || text == "ophelie-em" || text == "ophelie_em")
    {
        return FakeJouleMode::EmFixed;
    }
    if (text == "em-grid" || text == "em_grid" || text == "ophelie-em-grid" || text == "ophelie_em_grid")
    {
        return FakeJouleMode::EmGrid;
    }
    return FakeJouleMode::Off;
}

inline const char *fakeJouleModeName(FakeJouleMode mode)
{
    switch (mode)
    {
    case FakeJouleMode::Uniform:
        return "uniform";
    case FakeJouleMode::AnalyticPosition:
        return "analytic-position";
    case FakeJouleMode::EmFixed:
        return "em-fixed";
    case FakeJouleMode::EmGrid:
        return "em-grid";
    default:
        return "off";
    }
}

struct FakeJouleHeatGeometry
{
    Vecd center = Vecd::Zero();
    Real sigma_xy = Real(0.1);
    Real length_z = Real(0.1);
};

/** Q_i = shape_i with sum(Q*V) = target_power after normalization. */
class AssignAnalyticFakeJouleHeatShapeCK : public LocalDynamics
{
  public:
    AssignAnalyticFakeJouleHeatShapeCK(SPHBody &sph_body, const FakeJouleHeatGeometry &geom)
        : LocalDynamics(sph_body), center_(geom.center), inv_two_sigma2_(Real(1.0) / (2.0 * geom.sigma_xy * geom.sigma_xy + TinyReal)),
          half_length_z_(Real(0.5) * std::max(geom.length_z, Real(1.0e-6))),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField)),
          dv_pos_(particles_->template getVariableByName<Vecd>("Position"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : center_(encloser.center_), inv_two_sigma2_(encloser.inv_two_sigma2_),
              half_length_z_(encloser.half_length_z_), joule_(encloser.dv_joule_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd rel = pos_[index_i] - center_;
            const Real r2_xy = rel[0] * rel[0] + rel[1] * rel[1];
            const Real radial = std::exp(-r2_xy * inv_two_sigma2_);
            const Real z_weight = Real(0.35) + Real(0.65) * std::exp(-std::abs(rel[2]) / (half_length_z_ + TinyReal));
            joule_[index_i] = radial * z_weight;
        }

      protected:
        Vecd center_;
        Real inv_two_sigma2_;
        Real half_length_z_;
        Real *joule_;
        Vecd *pos_;
    };

  protected:
    Vecd center_;
    Real inv_two_sigma2_;
    Real half_length_z_;
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Vecd> *dv_pos_;
};

class ScaleFakeJouleHeatCK : public LocalDynamics
{
  public:
    explicit ScaleFakeJouleHeatCK(SPHBody &sph_body, Real scale)
        : LocalDynamics(sph_body), scale_(scale),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField))
    {
    }

    void setScale(Real scale) { scale_ = scale; }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : scale_(encloser.scale_), joule_(encloser.dv_joule_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            joule_[index_i] *= scale_;
        }

      protected:
        Real scale_;
        Real *joule_;
    };

  protected:
    Real scale_;
    DiscreteVariable<Real> *dv_joule_;
};

/** Explicit source: T += Q * dt / (rho0 * cp). Uses reference density for consistent energy bookkeeping. */
class ApplyFakeJouleHeatToTemperatureCK : public LocalDynamics
{
  public:
    ApplyFakeJouleHeatToTemperatureCK(SPHBody &sph_body, Real cp, Real rho0)
        : LocalDynamics(sph_body), cp_(cp), rho0_(rho0),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField)),
          dv_temperature_(particles_->template getVariableByName<Real>("Temperature"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : cp_(encloser.cp_), rho0_(encloser.rho0_), joule_(encloser.dv_joule_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            const Real rho_cp = rho0_ * cp_ + TinyReal;
            temperature_[index_i] += joule_[index_i] * dt / rho_cp;
        }

      protected:
        Real cp_;
        Real rho0_;
        Real *joule_;
        Real *temperature_;
    };

  protected:
    Real cp_;
    Real rho0_;
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Real> *dv_temperature_;
};

class FakeJoulePowerReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    FakeJoulePowerReduceCK(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "Rh200FakeJoulePower";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : joule_(encloser.dv_joule_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return joule_[index_i] * vol_[index_i];
        }

      protected:
        Real *joule_;
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Real> *dv_vol_;
};

class FakeJouleRawPowerReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    FakeJouleRawPowerReduceCK(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "Rh200FakeJouleRawPower";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : joule_(encloser.dv_joule_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return joule_[index_i] * vol_[index_i];
        }

      protected:
        Real *joule_;
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Real> *dv_vol_;
};

class ThermalEnergyAboveInitialReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    ThermalEnergyAboveInitialReduceCK(SPHBody &sph_body, Real t_initial, Real cp, Real rho0)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), t_initial_(t_initial), cp_(cp), rho0_(rho0),
          dv_temperature_(particles_->template getVariableByName<Real>("Temperature")),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "Rh200ThermalEnergyAboveInitial";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : t_initial_(encloser.t_initial_), cp_(encloser.cp_), rho0_(encloser.rho0_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return cp_ * rho0_ * (temperature_[index_i] - t_initial_) * vol_[index_i];
        }

      protected:
        Real t_initial_;
        Real cp_;
        Real rho0_;
        Real *temperature_;
        Real *vol_;
    };

  protected:
    Real t_initial_;
    Real cp_;
    Real rho0_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_vol_;
};

class TemperatureMaxReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    explicit TemperatureMaxReduceCK(SPHBody &sph_body)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>("Temperature"))
    {
        quantity_name_ = "Rh200TemperatureMax";
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

class TemperatureVolWeightedMeanReduceCK : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    explicit TemperatureVolWeightedMeanReduceCK(SPHBody &sph_body)
        : BaseDynamicsType(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>("Temperature")),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "Rh200TemperatureVolWeightedMean";
    }

    class FinishDynamics
    {
      public:
        using OutputType = Real;
        template <class EncloserType>
        FinishDynamics(EncloserType &)
        {
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
            : temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return ReduceReturnType(temperature_[index_i] * vol_[index_i], vol_[index_i]);
        }

      protected:
        Real *temperature_;
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_vol_;
};

inline FakeJouleHeatGeometry fakeJouleHeatGeometryFromBounds(const BoundingBoxd &bounds)
{
    FakeJouleHeatGeometry geom;
    const Vecd bmin(bounds.lower_[0], bounds.lower_[1], bounds.lower_[2]);
    const Vecd bmax(bounds.upper_[0], bounds.upper_[1], bounds.upper_[2]);
    const Vecd size = bmax - bmin;
    geom.center = Real(0.5) * (bmin + bmax);
    geom.sigma_xy = Real(0.22) * std::max(size[0], size[1]);
    geom.length_z = size[2];
    return geom;
}

inline void writeEnergyBudgetCsvHeader(const std::string &path)
{
    std::ofstream out(path, std::ios::out);
    out << "time,P_joule_W,thermal_energy_J,E_joule_integrated_J,T_mean_K,T_max_K,P_wall_loss_W\n";
}

inline void appendEnergyBudgetCsv(const std::string &path, Real time, Real p_joule, Real thermal_energy,
                                  Real e_joule_integrated, Real t_mean, Real t_max)
{
    std::ofstream out(path, std::ios::app);
    out << time << "," << p_joule << "," << thermal_energy << "," << e_joule_integrated << "," << t_mean << ","
        << t_max << ",0\n";
}

template <class AssignDynamics, class ScaleDynamics, class RawPowerReduceDynamics>
inline Real normalizeAnalyticFakeJouleHeat(SPHBody &glass, Real target_power_w, AssignDynamics &assign_shape,
                                           ScaleDynamics &scale_joule, RawPowerReduceDynamics &raw_power_reduce)
{
    (void)glass;
    assign_shape.exec();
    const Real raw_power = raw_power_reduce.exec();
    const Real scale = target_power_w / (raw_power + TinyReal);
    scale_joule.setScale(scale);
    scale_joule.exec();
    return scale;
}

template <class ExecutionPolicy>
inline void registerUniformFakeJouleHeat(FluidBody &glass, Real target_power_w)
{
    BaseParticles &particles = glass.getBaseParticles();
    particles.registerStateVariable<Real>(kJouleHeatField, Real(0));
    const size_t n = particles.TotalRealParticles();
    Real total_vol = Real(0);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < n; ++i)
    {
        total_vol += vol[i];
    }
    const Real q_uniform = target_power_w / (total_vol + TinyReal);
    Real *joule = particles.getVariableDataByName<Real>(kJouleHeatField);
    for (size_t i = 0; i < n; ++i)
    {
        joule[i] = q_uniform;
    }
#if SPHINXSYS_USE_SYCL
    (void)particles.getVariableByName<Real>(kJouleHeatField)->DelegatedData(ExecutionPolicy{});
#endif
    std::cout << "[rh200] fake joule uniform: target_power=" << target_power_w << " W, total_vol=" << total_vol
              << " m^3, q_uniform=" << q_uniform << " W/m^3" << std::endl;
}

} // namespace rh200
} // namespace SPH

#endif // RH200_FAKE_JOULE_HEAT_H
