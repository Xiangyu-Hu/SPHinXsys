#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Fixed-capacity, device-copyable temperature law.
 *
 * The digitized Jacoutot curves must be installed explicitly by the case
 * configuration. A one-point law is intentionally constant and is useful for
 * controlled regressions only; it is not a substitute for the paper's Fig. 2.
 */
class OphelieTemperatureLaw
{
  public:
    static constexpr size_t kMaxPoints = 16;

    std::array<Real, kMaxPoints> temperature_k_{};
    std::array<Real, kMaxPoints> value_{};
    size_t point_count_ = 0;
    bool logarithmic_value_ = false;

    void clear()
    {
        point_count_ = 0;
        logarithmic_value_ = false;
    }

    void setPoint(size_t index, Real temperature_k, Real value)
    {
        if (index >= kMaxPoints)
        {
            return;
        }
        temperature_k_[index] = temperature_k;
        value_[index] = value;
        point_count_ = std::max(point_count_, index + 1);
    }

    Real evaluate(Real temperature_k) const
    {
        if (point_count_ == 0)
        {
            return Real(0);
        }
        if (point_count_ == 1 || temperature_k <= temperature_k_[0])
        {
            return value_[0];
        }
        const size_t last = point_count_ - 1;
        if (temperature_k >= temperature_k_[last])
        {
            return value_[last];
        }
        for (size_t i = 0; i < last; ++i)
        {
            if (temperature_k <= temperature_k_[i + 1])
            {
                const Real t0 = temperature_k_[i];
                const Real t1 = temperature_k_[i + 1];
                const Real weight = (temperature_k - t0) / (t1 - t0 + TinyReal);
                if (logarithmic_value_ && value_[i] > TinyReal && value_[i + 1] > TinyReal)
                {
                    return std::exp((Real(1) - weight) * std::log(value_[i]) + weight * std::log(value_[i + 1]));
                }
                return (Real(1) - weight) * value_[i] + weight * value_[i + 1];
            }
        }
        return value_[last];
    }
};

struct OphelieFrenchGlassMaterialLaws
{
    OphelieTemperatureLaw rho;
    OphelieTemperatureLaw beta;
    OphelieTemperatureLaw mu;
    OphelieTemperatureLaw cp;
    OphelieTemperatureLaw k;
    OphelieTemperatureLaw sigma;
    bool paper_digitized_ = false;
};

/**
 * Exact Table-1 point from Jacoutot et al. (2008), T=1473 K.
 *
 * This is deliberately a constant regression law. Production French cases
 * should prefer thesis σ(T) laws (Fig. III.12 / IV.10) via
 * makeJacoutotThesisSigmaTemperatureLaw* and set paper_digitized_=true.
 *
 * Note: CEP paper Figure 2 is the coupling flowchart, not σ(T).
 */
inline OphelieFrenchGlassMaterialLaws makeJacoutot1473KReferenceMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws;
    laws.rho.setPoint(0, Real(1473), Real(2750));
    laws.beta.setPoint(0, Real(1473), Real(1.0e-5));
    laws.mu.setPoint(0, Real(1473), Real(4));
    laws.cp.setPoint(0, Real(1473), Real(1150));
    laws.k.setPoint(0, Real(1473), Real(4));
    laws.sigma.setPoint(0, Real(1473), Real(16));
    laws.mu.logarithmic_value_ = true;
    laws.sigma.logarithmic_value_ = true;
    return laws;
}

/**
 * Provisional Arrhenius-like sigma(T) anchored at Jacoutot Table-1:
 *   sigma(1473 K)=16 S/m, Ea/R = 15000 K.
 *
 * This is NOT thesis Fig. III.12 / IV.10. It exists so Stage 3.2 can exercise
 * T→σ→EM→Q→T machinery before switching to the thesis closed-form laws.
 */
inline OphelieTemperatureLaw makeProvisionalJacoutotLikeSigmaTemperatureLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real t_ref = Real(1473);
    constexpr Real sigma_ref = Real(16);
    constexpr Real ea_over_r = Real(15000);
    const Real temperatures[5] = {Real(1273), Real(1373), Real(1473), Real(1573), Real(1673)};
    for (size_t i = 0; i < 5; ++i)
    {
        const Real t = temperatures[i];
        const Real value = sigma_ref * std::exp(ea_over_r * (Real(1) / t_ref - Real(1) / t));
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * Jacoutot thesis Fig. III.12 / law (III.6) [C_2], natural / static bath track:
 *   log10(sigma) = 3.7921 - 3179.8 / T
 * Sampled onto OphelieTemperatureLaw with log-space interpolation.
 * At 1473 K this gives ~43 S/m (Table-1's 16 is a representative point, not this law).
 *
 * GPT 2026-08-11: raw III.12 is sensitivity-only for French Natural CEP2008 reproduction.
 */
inline OphelieTemperatureLaw makeJacoutotThesisSigmaTemperatureLawIII12()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real a = Real(3.7921);
    constexpr Real b = Real(3179.8);
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1229), Real(1300),
                                   Real(1400), Real(1473), Real(1500), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        const Real value = std::pow(Real(10), a - b / t);
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * French Natural (CEP 2008) production σ(T):
 *   journal Table-1 anchor σ(1473K)=16, thesis III.12 temperature shape, clipped to [1e-6, 30].
 *   sigma_nat(T) = 16 * 10^(3179.8 * (1/1473 - 1/T))
 */
inline OphelieTemperatureLaw makeFrenchNaturalJournalAnchoredSigmaTemperatureLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real t_ref = Real(1473);
    constexpr Real sigma_ref = Real(16);
    constexpr Real b = Real(3179.8);
    constexpr Real sigma_min = Real(1.0e-6);
    constexpr Real sigma_max = Real(30);
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1550), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        const Real raw = sigma_ref * std::pow(Real(10), b * (Real(1) / t_ref - Real(1) / t));
        const Real value = std::min(sigma_max, std::max(sigma_min, raw));
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

/**
 * Jacoutot thesis Fig. IV.10, stirred Uox2 + RuO2 glass [C_8]:
 *   log10(sigma) = 4.05726 - 3923.73 / T
 */
inline OphelieTemperatureLaw makeJacoutotThesisSigmaTemperatureLawIV10()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    constexpr Real a = Real(4.05726);
    constexpr Real b = Real(3923.73);
    const Real temperatures[11] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 11; ++i)
    {
        const Real t = temperatures[i];
        const Real value = std::pow(Real(10), a - b / t);
        sigma.setPoint(i, t, value);
    }
    return sigma;
}

inline OphelieFrenchGlassMaterialLaws makeProvisionalJacoutotLikeMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeProvisionalJacoutotLikeSigmaTemperatureLaw();
    laws.paper_digitized_ = false;
    return laws;
}

inline OphelieFrenchGlassMaterialLaws makeJacoutotThesisIII12MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeJacoutotThesisSigmaTemperatureLawIII12();
    laws.paper_digitized_ = true;
    return laws;
}

/** CEP2008 Natural production material: Table-1 constants + journal-anchored σ(T). */
inline OphelieFrenchGlassMaterialLaws makeFrenchNaturalCep2008MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeFrenchNaturalJournalAnchoredSigmaTemperatureLaw();
    laws.paper_digitized_ = true;
    return laws;
}

inline OphelieFrenchGlassMaterialLaws makeJacoutotThesisIV10MaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeJacoutot1473KReferenceMaterialLaws();
    laws.sigma = makeJacoutotThesisSigmaTemperatureLawIV10();
    laws.paper_digitized_ = true;
    return laws;
}

inline constexpr const char *kOphelieSigmaRawTemperatureField = "OphelieSigmaRawTemperature";

inline void registerOphelieFrenchTemperatureMaterialFields(BaseParticles &particles, Real initial_sigma)
{
    particles.registerStateVariable<Real>(kOphelieSigmaRawTemperatureField, initial_sigma);
}

inline void hostOphelieSigmaFieldStats(BaseParticles &particles, const std::string &sigma_field, size_t n, Real &sigma_min,
                                       Real &sigma_max, Real &sigma_mean)
{
    syncVariableToHost<Real>(particles, sigma_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *sigma = particles.getVariableDataByName<Real>(sigma_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    sigma_min = sigma[0];
    sigma_max = sigma[0];
    Real vol_sum = 0.0;
    Real weighted = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sigma_min = std::min(sigma_min, sigma[i]);
        sigma_max = std::max(sigma_max, sigma[i]);
        vol_sum += vol[i];
        weighted += sigma[i] * vol[i];
    }
    sigma_mean = weighted / (vol_sum + TinyReal);
}

inline Real hostOphelieRelativeScalarFieldChange(BaseParticles &particles, const std::string &field_name,
                                                 const StdVec<Real> &previous, size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    const Real *field = particles.getVariableDataByName<Real>(field_name);
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real delta = field[i] - previous[i];
        numerator += delta * delta;
        denominator += field[i] * field[i];
    }
    return std::sqrt(numerator / (denominator + TinyReal));
}

inline void hostStoreScalarField(BaseParticles &particles, const std::string &field_name, StdVec<Real> &storage,
                                 size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    const Real *field = particles.getVariableDataByName<Real>(field_name);
    storage.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        storage[i] = field[i];
    }
}

/**
 * Applies sigma(T) with pointwise under-relaxation:
 * sigma^(n+1)=(1-alpha)sigma^n+alpha*sigma_raw(T).
 */
class UpdateOphelieGlassSigmaFromTemperatureCK : public LocalDynamics
{
  public:
    UpdateOphelieGlassSigmaFromTemperatureCK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
                                             const OphelieTemperatureLaw &sigma_law, Real relaxation_alpha,
                                             const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), sigma_law_(sigma_law),
          relaxation_alpha_(std::max(Real(0), std::min(relaxation_alpha, Real(1)))),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
          dv_sigma_raw_(particles_->template getVariableByName<Real>(kOphelieSigmaRawTemperatureField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_law_(encloser.sigma_law_), relaxation_alpha_(encloser.relaxation_alpha_),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              sigma_raw_(encloser.dv_sigma_raw_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            const Real raw = sigma_law_.evaluate(temperature_[index_i]);
            sigma_raw_[index_i] = std::max(raw, TinyReal);
            sigma_[index_i] = (Real(1) - relaxation_alpha_) * sigma_[index_i] +
                              relaxation_alpha_ * sigma_raw_[index_i];
        }

      protected:
        OphelieTemperatureLaw sigma_law_;
        Real relaxation_alpha_;
        Real *temperature_;
        Real *sigma_;
        Real *sigma_raw_;
    };

  protected:
    OphelieTemperatureLaw sigma_law_;
    Real relaxation_alpha_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_sigma_raw_;
};

/** Diagnostic seed: only Temperature is set; DeltaT stays 0 so Joule closure remains heating-only. */
class SeedOphelieRadialTemperatureCK : public LocalDynamics
{
  public:
    SeedOphelieRadialTemperatureCK(SPHBody &sph_body, const Vecd &center, Real glass_radius, Real t0, Real delta_t,
                                   const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), center_(center), glass_radius_(std::max(glass_radius, TinyReal)), t0_(t0),
          delta_t_(delta_t), dv_pos_(particles_->template getVariableByName<Vecd>("Position")),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : center_(encloser.center_), glass_radius_(encloser.glass_radius_), t0_(encloser.t0_),
              delta_t_(encloser.delta_t_), pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = Real(0))
        {
            (void)dt;
            const Real dx = pos_[index_i][0] - center_[0];
            const Real dy = pos_[index_i][1] - center_[1];
            const Real r = std::sqrt(dx * dx + dy * dy);
            const Real frac = std::min(r / glass_radius_, Real(1));
            temperature_[index_i] = t0_ + delta_t_ * frac;
        }

      protected:
        Vecd center_;
        Real glass_radius_;
        Real t0_;
        Real delta_t_;
        Vecd *pos_;
        Real *temperature_;
    };

  protected:
    Vecd center_;
    Real glass_radius_;
    Real t0_;
    Real delta_t_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_temperature_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_MATERIAL_LAWS_H
