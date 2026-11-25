#ifndef SPDB_ARRHENIUS_CK_HPP
#define SPDB_ARRHENIUS_CK_HPP

#include "base_data_type.h"

namespace spdb
{
constexpr double kLn2 = 0.693147180559945309417;
constexpr double kExpLimit = 700.0;

inline double log_sinh(double x)
{
    if (x < 1.0e-4)
    {
        const double x2 = x * x;
        return SPH::math::log(x) + SPH::math::log1p(x2 / 6.0 + x2 * x2 / 120.0);
    }
    if (x > 20.0)
        return x - kLn2 + SPH::math::log1p(-SPH::math::exp(-2.0 * x));
    return SPH::math::log(SPH::math::sinh(x));
}

inline double asinh_exp(double z)
{
    if (z > 0.0)
        return z + SPH::math::log(1.0 + SPH::math::sqrt(1.0 + SPH::math::exp(-2.0 * z)));
    if (z < -40.0)
        return SPH::math::exp(z);
    const double y = SPH::math::exp(z);
    return SPH::math::log(y + SPH::math::sqrt(y * y + 1.0));
}

inline void flag_window(double temperature_K, double strain_eval, double inelastic_strain, StartupPolicy policy,
                        ArrheniusStatus *status)
{
    if (status == nullptr)
        return;
    if (temperature_K < Ren2014X20Cr13::fit_temperature_min_K || temperature_K > Ren2014X20Cr13::fit_temperature_max_K)
    {
        status->temperature_outside = true;
        status->ok = false;
        status->diagnostic = "temperature_outside_1173_1423K";
    }
    if (inelastic_strain > Ren2014X20Cr13::fit_strain_max)
    {
        status->strain_above_fit = true;
        status->ok = false;
        status->diagnostic = "strain_above_0.75";
    }
    if (policy == StartupPolicy::EndpointExtend && inelastic_strain < Ren2014X20Cr13::fit_strain_min)
        status->endpoint_extension = true;
    if (policy == StartupPolicy::PolyToZero && strain_eval < Ren2014X20Cr13::fit_strain_min)
        status->poly_below_fit = true;
}

inline double Ren2014X20Cr13::horner(const double coeff[6], double strain)
{
    double value = 0.0;
    for (int k = 5; k >= 0; --k)
        value = value * strain + coeff[k];
    return value;
}

inline ArrheniusCoeffs Ren2014X20Cr13::coefficients(double strain_eval)
{
    ArrheniusCoeffs c;
    c.alpha_per_MPa = horner(alpha_MPa, strain_eval);
    c.n = horner(n_coeff, strain_eval);
    c.Q_J_per_mol = horner(Q_coeff, strain_eval);
    c.lnA = horner(lnA_coeff, strain_eval);
    return c;
}

inline double Ren2014X20Cr13::evaluation_strain(double inelastic_strain, StartupPolicy policy)
{
    if (policy == StartupPolicy::EndpointExtend)
        return inelastic_strain > fit_strain_min ? inelastic_strain : fit_strain_min;
    return inelastic_strain;
}

inline double Ren2014X20Cr13::inelastic_rate(double temperature_K, double equivalent_stress_MPa, double inelastic_strain,
                                             StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status)
{
    ArrheniusStatus local;
    ArrheniusStatus *st = status != nullptr ? status : &local;
    *st = ArrheniusStatus{};
    if (inelastic_strain < 0.0)
    {
        st->ok = false;
        st->diagnostic = "negative_inelastic_strain";
        return 0.0;
    }
    if (equivalent_stress_MPa == 0.0)
    {
        flag_window(temperature_K, evaluation_strain(inelastic_strain, policy), inelastic_strain, policy, st);
        st->diagnostic = st->ok ? "zero_stress" : st->diagnostic;
        return 0.0;
    }
    if (equivalent_stress_MPa < 0.0)
    {
        st->ok = false;
        st->diagnostic = "negative_equivalent_stress";
        return 0.0;
    }
    const double strain_eval = evaluation_strain(inelastic_strain, policy);
    flag_window(temperature_K, strain_eval, inelastic_strain, policy, st);
    if (st->strain_above_fit)
        return 0.0;
    const ArrheniusCoeffs c = coefficients(strain_eval);
    if (!(c.alpha_per_MPa > 0.0) || !(c.n > 0.0))
    {
        st->ok = false;
        st->diagnostic = "nonpositive_alpha_or_n";
        return 0.0;
    }
    const double x = c.alpha_per_MPa * equivalent_stress_MPa;
    const double ln_rate = c.lnA - c.Q_J_per_mol / (gas_constant_J_per_mol_K * temperature_K) + c.n * log_sinh(x);
    if (ln_rate < -kExpLimit)
    {
        st->rate_underflow = true;
        st->diagnostic = "rate_underflow_recorded_as_zero";
        return 0.0;
    }
    const double rate = SPH::math::exp(ln_rate);
    if (rate < fit_strain_rate_min_per_s)
    {
        st->low_rate_extrapolation = true;
        if (!allow_low_rate_extrapolation)
        {
            st->ok = false;
            st->diagnostic = "rate_below_0.01_extrapolation_refused";
        }
        else if (st->ok)
            st->diagnostic = "low_rate_extrapolation";
    }
    return rate;
}

inline double Ren2014X20Cr13::inelastic_rate_Pa(double temperature_K, double equivalent_stress_Pa, double inelastic_strain,
                                                StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status)
{
    return inelastic_rate(temperature_K, equivalent_stress_Pa * 1.0e-6, inelastic_strain, policy,
                          allow_low_rate_extrapolation, status);
}

inline double Ren2014X20Cr13::flow_stress_MPa(double temperature_K, double strain_rate_per_s, double inelastic_strain,
                                              StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status)
{
    ArrheniusStatus local;
    ArrheniusStatus *st = status != nullptr ? status : &local;
    *st = ArrheniusStatus{};
    if (!(strain_rate_per_s > 0.0))
    {
        st->ok = false;
        st->diagnostic = "nonpositive_strain_rate";
        return 0.0;
    }
    if (inelastic_strain < 0.0)
    {
        st->ok = false;
        st->diagnostic = "negative_inelastic_strain";
        return 0.0;
    }
    const double strain_eval = evaluation_strain(inelastic_strain, policy);
    flag_window(temperature_K, strain_eval, inelastic_strain, policy, st);
    if (strain_rate_per_s < fit_strain_rate_min_per_s || strain_rate_per_s > fit_strain_rate_max_per_s)
    {
        st->low_rate_extrapolation = strain_rate_per_s < fit_strain_rate_min_per_s;
        if (strain_rate_per_s < fit_strain_rate_min_per_s && !allow_low_rate_extrapolation)
        {
            st->ok = false;
            st->diagnostic = "rate_below_0.01_extrapolation_refused";
            return 0.0;
        }
        if (strain_rate_per_s > fit_strain_rate_max_per_s)
        {
            st->ok = false;
            st->diagnostic = "rate_above_10";
        }
    }
    if (st->strain_above_fit)
        return 0.0;
    const ArrheniusCoeffs c = coefficients(strain_eval);
    if (!(c.alpha_per_MPa > 0.0) || !(c.n > 0.0))
    {
        st->ok = false;
        st->diagnostic = "nonpositive_alpha_or_n";
        return 0.0;
    }
    const double z = (SPH::math::log(strain_rate_per_s) + c.Q_J_per_mol / (gas_constant_J_per_mol_K * temperature_K) - c.lnA) / c.n;
    return asinh_exp(z) / c.alpha_per_MPa;
}

inline double Ren2014X20Cr13::flow_stress_Pa(double temperature_K, double strain_rate_per_s, double inelastic_strain,
                                             StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status)
{
    return flow_stress_MPa(temperature_K, strain_rate_per_s, inelastic_strain, policy, allow_low_rate_extrapolation, status) *
           1.0e6;
}

inline double Ren2014X20Cr13::von_mises_Pa(double sxx, double syy, double szz, double syz, double szx, double sxy)
{
    const double mean = (sxx + syy + szz) / 3.0;
    const double dxx = sxx - mean;
    const double dyy = syy - mean;
    const double dzz = szz - mean;
    const double second = dxx * dxx + dyy * dyy + dzz * dzz + 2.0 * (syz * syz + szx * szx + sxy * sxy);
    return SPH::math::sqrt(1.5 * second);
}

// Fixed total strain, J2: sigma_new = sigma_trial - 3 G * rate(sigma_new) * dt.
// Bounded bisection so the same function can run inside a SYCL kernel.
struct EquivalentReturn
{
    double sigma_eq_Pa = 0.0;
    double inelastic_eq = 0.0;
    double rate_per_s = 0.0;
    int accepted = 0;
};

inline EquivalentReturn equivalent_stress_backward_euler(double sigma_trial_Pa, double shear_modulus_Pa, double dt_s,
                                                         double temperature_K, double inelastic_eq, StartupPolicy policy)
{
    EquivalentReturn out;
    out.inelastic_eq = inelastic_eq;
    out.sigma_eq_Pa = sigma_trial_Pa > 0.0 ? sigma_trial_Pa : 0.0;
    if (!(sigma_trial_Pa > 0.0) || !(shear_modulus_Pa > 0.0) || !(dt_s > 0.0) || !(inelastic_eq >= 0.0))
    {
        out.accepted = sigma_trial_Pa >= 0.0 && inelastic_eq >= 0.0 ? 1 : 0;
        return out;
    }

    auto residual = [&](double sigma_Pa, double *rate_out) {
        const double p_trial_rate =
            Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_Pa, inelastic_eq, policy, true, nullptr);
        const double p_next = inelastic_eq + dt_s * p_trial_rate;
        const double rate =
            Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_Pa, p_next, policy, true, nullptr);
        if (rate_out != nullptr)
            *rate_out = rate;
        return sigma_Pa - sigma_trial_Pa + 3.0 * shear_modulus_Pa * dt_s * rate;
    };

    double unused_rate = 0.0;
    const double f_hi = residual(sigma_trial_Pa, &unused_rate);
    if (inelastic_eq + dt_s * unused_rate > Ren2014X20Cr13::fit_strain_max)
        return out;
    double lo = 0.0;
    double hi = sigma_trial_Pa;
    if (f_hi < 0.0)
        return out;
    for (int it = 0; it < 60; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        const double f_mid = residual(mid, nullptr);
        if (f_mid > 0.0)
            hi = mid;
        else
            lo = mid;
    }
    const double sigma = 0.5 * (lo + hi);
    double rate = 0.0;
    residual(sigma, &rate);
    const double p_next = inelastic_eq + dt_s * rate;
    if (p_next > Ren2014X20Cr13::fit_strain_max)
        return out;
    out.sigma_eq_Pa = sigma;
    out.inelastic_eq = p_next;
    out.rate_per_s = rate;
    out.accepted = 1;
    return out;
}

} // namespace spdb

#endif
