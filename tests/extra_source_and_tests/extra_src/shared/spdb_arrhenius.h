/**
 * @file spdb_arrhenius.h
 * @brief Ren et al. 2014 X20Cr13 hot-compression Arrhenius model.
 *
 * Literature reference material. Not a calibrated 7200 s creep law.
 * Coefficients: Ren F C, Chen J, Chen F, Trans. Nonferrous Met. Soc. China
 * 24 (2014) 1407–1413. DOI 10.1016/S1003-6326(14)63206-4.
 */
#ifndef SPDB_ARRHENIUS_H
#define SPDB_ARRHENIUS_H

namespace spdb
{

enum class StartupPolicy
{
    EndpointExtend = 0, // evaluate coefficients at max(strain, 0.05); history stays unchanged
    PolyToZero = 1      // evaluate the printed polynomial below 0.05
};

struct ArrheniusStatus
{
    bool ok = true;
    bool low_rate_extrapolation = false;
    bool endpoint_extension = false;
    bool poly_below_fit = false;
    bool strain_above_fit = false;
    bool temperature_outside = false;
    bool rate_underflow = false;
    const char *diagnostic = "ok";
};

struct ArrheniusCoeffs
{
    double alpha_per_MPa = 0.0;
    double n = 0.0;
    double Q_J_per_mol = 0.0;
    double lnA = 0.0;
};

class Ren2014X20Cr13
{
  public:
    static constexpr const char *model_name = "ren2014_x20cr13";
    static constexpr const char *doi = "10.1016/S1003-6326(14)63206-4";
    static constexpr double gas_constant_J_per_mol_K = 8.314;
    static constexpr double fit_temperature_min_K = 1173.0;
    static constexpr double fit_temperature_max_K = 1423.0;
    static constexpr double fit_strain_rate_min_per_s = 0.01;
    static constexpr double fit_strain_rate_max_per_s = 10.0;
    static constexpr double fit_strain_min = 0.05;
    static constexpr double fit_strain_max = 0.75;

    // p(ε) = c0 + c1 ε + ... + c5 ε^5. α is per MPa. Q is J/mol. ln A uses natural log, A in 1/s.
    static constexpr double alpha_MPa[6] = {0.01277, -0.04834, 0.18608, -0.36932, 0.37353, -0.15155};
    static constexpr double n_coeff[6] = {7.19785, -12.73951, 57.97832, -145.77045, 168.05383, -71.91728};
    static constexpr double Q_coeff[6] = {454181.0, -319488.0, -28978.6, 1777650.0, -3257280.0, 1745480.0};
    static constexpr double lnA_coeff[6] = {40.18281, -29.51957, 2.19563, 134.56396, -247.06753, 130.34468};

    static double horner(const double coeff[6], double strain);
    static ArrheniusCoeffs coefficients(double strain_eval);
    static double evaluation_strain(double inelastic_strain, StartupPolicy policy);

    // Equivalent stress is nonnegative. Compression sign must not enter sinh.
    static double flow_stress_MPa(double temperature_K, double strain_rate_per_s, double inelastic_strain,
                                  StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status);
    static double flow_stress_Pa(double temperature_K, double strain_rate_per_s, double inelastic_strain,
                                 StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status);
    static double inelastic_rate(double temperature_K, double equivalent_stress_MPa, double inelastic_strain,
                                 StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status);
    static double inelastic_rate_Pa(double temperature_K, double equivalent_stress_Pa, double inelastic_strain,
                                    StartupPolicy policy, bool allow_low_rate_extrapolation, ArrheniusStatus *status);

    // von Mises equivalent of a Cauchy stress. Pure pressure returns 0.
    static double von_mises_Pa(double sxx, double syy, double szz, double syz, double szx, double sxy);
};

struct HoldIntegration
{
    double final_strain = 0.0;
    double final_rate_per_s = 0.0;
    double time_reached_s = 0.0;
    int accepted_steps = 0;
    int stopped = 0;
    double low_rate_step_fraction = 0.0;
    const char *stop_reason = "completed";
};

HoldIntegration integrate_constant_stress(double temperature_K, double equivalent_stress_MPa, double duration_s,
                                          double dt_s, StartupPolicy policy, bool allow_low_rate_extrapolation);

// Small-strain uniaxial compression branch of J2:
// sigma_zz = E * (eps_zz - eps_in_zz), eps_in_zz = -p, p_dot = inelastic_rate(|sigma_zz|).
// One backward-Euler step. The caller commits the returned history only when accepted == 1.
// youngs_modulus_Pa is a demo_assumption, not a Ren 2014 value.
struct UniaxialState
{
    double eps_zz = 0.0;
    double inelastic_eq = 0.0;
    double eps_in_zz = 0.0;
    double sigma_zz_Pa = 0.0;
    double sigma_eq_Pa = 0.0;
    double inelastic_rate_per_s = 0.0;
    int iterations = 0;
    int accepted = 0;
    const char *diagnostic = "not_run";
};

UniaxialState uniaxial_backward_euler(double temperature_K, double youngs_modulus_Pa, double eps_zz,
                                      double inelastic_eq, double dt_s, StartupPolicy policy);

// Fixed Cauchy stress. E is not used to update p; it only sets the reported elastic strain.
UniaxialState uniaxial_constant_stress_step(double temperature_K, double youngs_modulus_Pa, double sigma_zz_Pa,
                                            double inelastic_eq, double dt_s, StartupPolicy policy);

int run_ren2014_material_point();

} // namespace spdb

// Inline rate evaluation. Included by SYCL ConstituteKernel translation units.
#include "spdb_arrhenius_ck.hpp"

#endif
