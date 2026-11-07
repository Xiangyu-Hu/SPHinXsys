/**
 * @file electromagnetic_ophelie_boussinesq.h
 * @brief Stage 4.1 Boussinesq buoyancy force: F = -m β (T − T0) g.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H
#define ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H

#include "force_prior_ck.h"
#include "force_prior_ck.hpp"
#include "external_force.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "sphinxsys_variable.h"

#include <cmath>
#include <iostream>
#include <memory>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Boussinesq buoyancy as ForcePrior companion to GravityForceCK.
 * F_b = −mass · β · (T − T_ref) · g_vector
 * with g_vector = Gravity::InducedAcceleration (usually (0,0,−g)).
 */
class BoussinesqBuoyancyForceCK : public LocalDynamics, public ForcePriorCK
{
  public:
    BoussinesqBuoyancyForceCK(SPHBody &sph_body, const Gravity &gravity, Real beta, Real t_ref,
                              const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), ForcePriorCK(this->particles_, "BoussinesqBuoyancyForceCK"), gravity_(gravity),
          beta_(beta), t_ref_(t_ref),
          sv_force_tref_(std::make_unique<SingleVariable<Real>>("OphelieBoussinesqForceTref", t_ref)),
          sv_physical_time_(&sph_system_->svPhysicalTime()),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_mass_(particles_->getVariableByName<Real>("Mass")),
          dv_temperature_(particles_->getVariableByName<Real>(temperature_field))
    {
    }

    /**
     * Doc189: force reference is a device-shared scalar. The cached UpdateKernel stores the pointer
     * (not a copied Real), so later setValue is visible on the next exec without rebuilding ForcePrior.
     * EOS/material Tref is a different quantity and is not changed here.
     */
    void setForceReferenceTemperature(Real t_ref)
    {
        t_ref_ = t_ref;
        sv_force_tref_->setValue(t_ref);
    }
    Real getForceReferenceTemperature() const { return sv_force_tref_->getValue(); }

    class UpdateKernel : public ForcePriorCK::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : ForcePriorCK::UpdateKernel(ex_policy, encloser), gravity_(encloser.gravity_), beta_(encloser.beta_),
              t_ref_(encloser.sv_force_tref_->DelegatedData(ex_policy)),
              physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)), mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd g_acc = gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
            const Real t_ref = *t_ref_;
            this->current_force_[index_i] = -mass_[index_i] * beta_ * (temperature_[index_i] - t_ref) * g_acc;
            ForcePriorCK::UpdateKernel::update(index_i, dt);
        }

      protected:
        Gravity gravity_;
        Real beta_;
        Real *t_ref_;
        Real *physical_time_;
        Vecd *pos_;
        Real *mass_;
        Real *temperature_;
    };

  protected:
    Gravity gravity_;
    Real beta_;
    Real t_ref_;
    std::unique_ptr<SingleVariable<Real>> sv_force_tref_;
    SingleVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieVelocityMaxReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    explicit OphelieVelocityMaxReduceCK(SPHBody &sph_body)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_vel_(particles_->template getVariableByName<Vecd>("Velocity"))
    {
        quantity_name_ = "OphelieVelocityMax";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : vel_(encloser.dv_vel_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return vel_[index_i].norm();
        }

      protected:
        Vecd *vel_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_vel_;
};

class OphelieTemperatureMeanReduceCK : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    explicit OphelieTemperatureMeanReduceCK(SPHBody &sph_body,
                                            const std::string &temperature_field = "Temperature")
        : BaseDynamicsType(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieTemperatureMean";
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

/**
 * P1-03: hot-bottom linear profile T = T0 + dT * (1 − ξ), ξ=(z−z_bottom)/H ∈ [0,1].
 * Syncs ThermalDeltaT = T − T0 when that field is present.
 */
inline void hostPrescribeHotBottomLinearTemperature(BaseParticles &particles, const Vecd &glass_center,
                                                    Real glass_half_height, Real t0, Real dT,
                                                    const std::string &temperature_field = "Temperature",
                                                    const std::string &delta_t_field = "ThermalDeltaT")
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, temperature_field);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    Real *delta_t = nullptr;
    try
    {
        syncVariableToHost<Real>(particles, delta_t_field);
        delta_t = particles.getVariableDataByName<Real>(delta_t_field);
    }
    catch (const std::exception &)
    {
        delta_t = nullptr;
    }
    const Real z_bottom = glass_center[2] - glass_half_height;
    const Real H = Real(2) * glass_half_height;
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real xi = std::min(Real(1), std::max(Real(0), (pos[i][2] - z_bottom) / (H + TinyReal)));
        const Real Ti = t0 + dT * (Real(1) - xi);
        temp[i] = Ti;
        if (delta_t != nullptr)
        {
            delta_t[i] = Ti - t0;
        }
    }
    syncVariableToDevice<Real>(particles, temperature_field);
    if (delta_t != nullptr)
    {
        syncVariableToDevice<Real>(particles, delta_t_field);
    }
}

/**
 * Doc77 §3: smooth radial ΔT (non-hydrostatic drive).
 * T = T0 + dT * (1 − r/R) — hot axis, cold rim (EM-like). Optional mild z taper:
 * multiply by 0.5*(1+cos(π(ξ−0.5))) so mid-height is strongest (optional via z_taper).
 */
inline void hostPrescribeRadialHotCoreTemperature(BaseParticles &particles, const Vecd &glass_center,
                                                  Real glass_radius, Real glass_half_height, Real t0, Real dT,
                                                  bool z_taper = false,
                                                  const std::string &temperature_field = "Temperature",
                                                  const std::string &delta_t_field = "ThermalDeltaT")
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, temperature_field);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    Real *delta_t = nullptr;
    try
    {
        syncVariableToHost<Real>(particles, delta_t_field);
        delta_t = particles.getVariableDataByName<Real>(delta_t_field);
    }
    catch (const std::exception &)
    {
        delta_t = nullptr;
    }
    const Real z_bottom = glass_center[2] - glass_half_height;
    const Real H = Real(2) * glass_half_height;
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - glass_center[0];
        const Real dy = pos[i][1] - glass_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real rn = std::min(Real(1), std::max(Real(0), r / (glass_radius + TinyReal)));
        Real amp = Real(1) - rn;
        if (z_taper)
        {
            const Real xi = std::min(Real(1), std::max(Real(0), (pos[i][2] - z_bottom) / (H + TinyReal)));
            amp *= Real(0.5) * (Real(1) + std::cos(Pi * (xi - Real(0.5))));
        }
        const Real Ti = t0 + dT * amp;
        temp[i] = Ti;
        if (delta_t != nullptr)
        {
            delta_t[i] = Ti - t0;
        }
    }
    syncVariableToDevice<Real>(particles, temperature_field);
    if (delta_t != nullptr)
    {
        syncVariableToDevice<Real>(particles, delta_t_field);
    }
}

struct OphelieBoussinesqForceAudit
{
    Real f_z_sum = 0.0;
    Real f_z_analytic = 0.0;
    Real rel_err = 0.0;
    Real mean_dT = 0.0;
    Real mass_sum = 0.0;
    size_t n = 0;
};

/** After BoussinesqBuoyancyForceCK::exec, audit Σ F_z vs −Σ m β (T−T0) g_z. */
inline OphelieBoussinesqForceAudit hostAuditBoussinesqForcePrior(
    BaseParticles &particles, Real beta, Real t_ref, Real g_signed_z,
    const std::string &temperature_field = "Temperature",
    const std::string &force_field = "BoussinesqBuoyancyForceCK")
{
    OphelieBoussinesqForceAudit a;
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, force_field);
    const Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *f = particles.getVariableDataByName<Vecd>(force_field);
    a.n = particles.TotalRealParticles();
    double fz = 0.0, fz_an = 0.0, mdT = 0.0, msum = 0.0;
    for (size_t i = 0; i < a.n; ++i)
    {
        const Real dT = temp[i] - t_ref;
        fz += f[i][2];
        fz_an += -mass[i] * beta * dT * g_signed_z;
        mdT += mass[i] * dT;
        msum += mass[i];
    }
    a.f_z_sum = static_cast<Real>(fz);
    a.f_z_analytic = static_cast<Real>(fz_an);
    a.mass_sum = static_cast<Real>(msum);
    a.mean_dT = msum > TinyReal ? static_cast<Real>(mdT / msum) : Real(0);
    a.rel_err = std::abs(a.f_z_sum - a.f_z_analytic) / (std::abs(a.f_z_analytic) + TinyReal);
    return a;
}

inline void logOphelieBoussinesqForceAudit(const OphelieBoussinesqForceAudit &a, const char *tag)
{
    std::cout << "[ophelie][boussinesq-audit] " << tag << " n=" << a.n << " mean_dT=" << a.mean_dT
              << " Fz_sum=" << a.f_z_sum << " Fz_analytic=" << a.f_z_analytic << " rel_err=" << a.rel_err
              << " mass=" << a.mass_sum << std::endl;
}

/**
 * Doc189 particle-wise buoyancy gate. Reads the existing dedicated force field.
 * Does not execute dynamics and does not write ForcePrior.
 * F_ref = -m β (T-Tf) g_acc, g_acc=(0,0,-g) with g>0.
 */
struct OphelieBuoyancyFieldGate
{
    double E_L2 = 0.0;
    double E_net = 0.0;
    double max_abs_N = 0.0;
    double net_abs_N = 0.0;
    double sum_fref_norm_N = 0.0;
    double sum_abs_fref_N = 0.0;
    double ffloor_N = 1.0e-8;
    double tref_requested_K = 0.0;
    double tref_inferred_K = 0.0;
    double tmean_K = 0.0;
    bool near_zero_force = false;
    bool algebraic_ok = false;
    bool dynamic_reference_ok = false;
};

inline OphelieBuoyancyFieldGate hostAuditBuoyancyFieldGate(BaseParticles &particles, Real beta, Real t_force,
                                                          Real g_positive,
                                                          const std::string &temperature_field = "Temperature")
{
    OphelieBuoyancyFieldGate g;
    g.tref_requested_K = static_cast<double>(t_force);
    g.ffloor_N = 1.0e-8;
    syncVariableToHost<Real>(particles, temperature_field);
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Vecd>(particles, "BoussinesqBuoyancyForceCK");
    const Real *temp = particles.getVariableDataByName<Real>(temperature_field);
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *f = particles.getVariableDataByName<Vecd>("BoussinesqBuoyancyForceCK");
    const size_t n = particles.TotalRealParticles();
    double sum_e2 = 0.0, sum_f2 = 0.0, sum_abs = 0.0;
    double ex = 0.0, ey = 0.0, ez = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;
    double wT = 0.0, w = 0.0, sum_fz = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double m = static_cast<double>(mass[i]);
        const double Ti = static_cast<double>(temp[i]);
        const double dT = Ti - static_cast<double>(t_force);
        const double fz_ref = m * static_cast<double>(beta) * static_cast<double>(g_positive) * dT;
        const double ei_x = static_cast<double>(f[i][0]);
        const double ei_y = static_cast<double>(f[i][1]);
        const double ei_z = static_cast<double>(f[i][2]) - fz_ref;
        const double e2 = ei_x * ei_x + ei_y * ei_y + ei_z * ei_z;
        sum_e2 += e2;
        sum_f2 += fz_ref * fz_ref;
        sum_abs += std::abs(fz_ref);
        g.max_abs_N = std::max(g.max_abs_N, std::sqrt(e2));
        ex += ei_x;
        ey += ei_y;
        ez += ei_z;
        fx += static_cast<double>(f[i][0]);
        fy += static_cast<double>(f[i][1]);
        fz += fz_ref;
        sum_fz += static_cast<double>(f[i][2]);
        wT += m * Ti;
        w += m;
    }
    g.sum_fref_norm_N = std::sqrt(sum_f2);
    g.sum_abs_fref_N = sum_abs;
    g.net_abs_N = std::sqrt(ex * ex + ey * ey + ez * ez);
    g.tmean_K = w > 0.0 ? wT / w : 0.0;
    const double denom_m = w * static_cast<double>(beta) * static_cast<double>(g_positive);
    g.tref_inferred_K = denom_m != 0.0 ? g.tmean_K - sum_fz / denom_m : g.tmean_K;
    const double sum_fref_vec = std::sqrt(fx * fx + fy * fy + fz * fz);
    g.near_zero_force = g.sum_fref_norm_N < g.ffloor_N;
    // Mean-removed buoyancy has a tiny net force. Do not divide that residual by a 1e-8 N floor.
    // Absolute net screen 1e-5 N is far below the old 0.007 N stale-reference residual.
    const double net_abs_tol_N = 1.0e-5;
    const bool net_reference_cancelled = sum_fref_vec < std::max(g.ffloor_N, 1.0e-4 * g.sum_fref_norm_N);
    if (g.near_zero_force)
    {
        g.E_L2 = std::sqrt(sum_e2);
        g.E_net = g.net_abs_N;
        g.algebraic_ok = g.max_abs_N <= 1.0e-6 && g.net_abs_N <= net_abs_tol_N;
    }
    else
    {
        g.E_L2 = std::sqrt(sum_e2 / std::max(sum_f2, g.ffloor_N * g.ffloor_N));
        if (net_reference_cancelled)
        {
            g.E_net = g.net_abs_N;
            g.algebraic_ok = g.E_L2 <= 1.0e-5 && g.net_abs_N <= net_abs_tol_N;
        }
        else
        {
            g.E_net = g.net_abs_N / std::max(sum_fref_vec, g.ffloor_N);
            // Relative net ratio is unstable when ||sum Fref|| is a tiny residual of cancelling particles.
            g.algebraic_ok = g.E_L2 <= 1.0e-5 && (g.E_net <= 1.0e-5 || g.net_abs_N <= net_abs_tol_N);
        }
    }
    g.dynamic_reference_ok = std::abs(g.tref_inferred_K - g.tref_requested_K) <= 1.0e-2;
    return g;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H
