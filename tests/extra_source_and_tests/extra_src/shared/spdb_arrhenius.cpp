#include "spdb_arrhenius.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace spdb
{
HoldIntegration integrate_constant_stress(double temperature_K, double equivalent_stress_MPa, double duration_s,
                                          double dt_s, StartupPolicy policy, bool allow_low_rate_extrapolation)
{
    HoldIntegration result;
    if (!(dt_s > 0.0) || !(duration_s > 0.0))
    {
        result.stopped = 1;
        result.stop_reason = "nonpositive_time";
        return result;
    }
    double time = 0.0;
    double strain = 0.0;
    int low_rate_steps = 0;
    while (time < duration_s - 0.5 * dt_s)
    {
        double step = dt_s;
        if (time + step > duration_s)
            step = duration_s - time;
        ArrheniusStatus status;
        const double rate = Ren2014X20Cr13::inelastic_rate(temperature_K, equivalent_stress_MPa, strain, policy,
                                                           allow_low_rate_extrapolation, &status);
        if (!status.ok && status.strain_above_fit)
        {
            result.stopped = 1;
            result.stop_reason = status.diagnostic;
            break;
        }
        if (!status.ok && !status.low_rate_extrapolation && !status.endpoint_extension && !status.poly_below_fit)
        {
            result.stopped = 1;
            result.stop_reason = status.diagnostic;
            break;
        }
        const double trial = strain + rate * step;
        if (trial > Ren2014X20Cr13::fit_strain_max)
        {
            result.stopped = 1;
            result.stop_reason = "strain_would_exceed_0.75";
            result.final_rate_per_s = rate;
            break;
        }
        strain = trial;
        time += step;
        result.final_rate_per_s = rate;
        result.accepted_steps++;
        if (status.low_rate_extrapolation)
            low_rate_steps++;
    }
    result.final_strain = strain;
    result.time_reached_s = time;
    result.low_rate_step_fraction =
        result.accepted_steps > 0 ? static_cast<double>(low_rate_steps) / static_cast<double>(result.accepted_steps) : 0.0;
    return result;
}

namespace
{
double compressive_sigma_eq_Pa(double youngs_modulus_Pa, double eps_zz, double inelastic_eq)
{
    return std::max(0.0, -youngs_modulus_Pa * (eps_zz + inelastic_eq));
}

UniaxialState rejected_uniaxial(double eps_zz, double inelastic_eq, const char *diagnostic)
{
    UniaxialState state;
    state.eps_zz = eps_zz;
    state.inelastic_eq = inelastic_eq;
    state.eps_in_zz = -inelastic_eq;
    state.diagnostic = diagnostic;
    state.accepted = 0;
    return state;
}
} // namespace

UniaxialState uniaxial_backward_euler(double temperature_K, double youngs_modulus_Pa, double eps_zz,
                                      double inelastic_eq, double dt_s, StartupPolicy policy)
{
    if (!(youngs_modulus_Pa > 0.0) || !(dt_s > 0.0) || !(inelastic_eq >= 0.0))
        return rejected_uniaxial(eps_zz, inelastic_eq, "invalid_uniaxial_input");
    if (eps_zz > 0.0)
        return rejected_uniaxial(eps_zz, inelastic_eq, "compression_branch_only");

    const double p_zero_stress = -eps_zz;
    auto residual = [&](double p, ArrheniusStatus *status) {
        const double sigma_eq = compressive_sigma_eq_Pa(youngs_modulus_Pa, eps_zz, p);
        return p - inelastic_eq -
               dt_s * Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_eq, p, policy, true, status);
    };

    ArrheniusStatus status;
    const double f_start = residual(inelastic_eq, &status);
    if (!status.ok && status.strain_above_fit)
        return rejected_uniaxial(eps_zz, inelastic_eq, status.diagnostic);
    if (!(p_zero_stress > inelastic_eq))
    {
        UniaxialState state;
        state.eps_zz = eps_zz;
        state.inelastic_eq = inelastic_eq;
        state.eps_in_zz = -inelastic_eq;
        state.sigma_zz_Pa = youngs_modulus_Pa * (eps_zz + inelastic_eq);
        state.sigma_eq_Pa = compressive_sigma_eq_Pa(youngs_modulus_Pa, eps_zz, inelastic_eq);
        state.inelastic_rate_per_s = 0.0;
        state.accepted = 1;
        state.diagnostic = "zero_or_tensile_trial_no_increment";
        return state;
    }
    if (f_start >= 0.0)
    {
        UniaxialState state;
        state.eps_zz = eps_zz;
        state.inelastic_eq = inelastic_eq;
        state.eps_in_zz = -inelastic_eq;
        state.sigma_zz_Pa = youngs_modulus_Pa * (eps_zz + inelastic_eq);
        state.sigma_eq_Pa = compressive_sigma_eq_Pa(youngs_modulus_Pa, eps_zz, inelastic_eq);
        state.inelastic_rate_per_s =
            Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, state.sigma_eq_Pa, inelastic_eq, policy, true, &status);
        state.accepted = 1;
        state.diagnostic = status.diagnostic;
        return state;
    }

    double lo = inelastic_eq;
    double hi = p_zero_stress;
    ArrheniusStatus hi_status;
    const double f_hi = residual(hi, &hi_status);
    if (!std::isfinite(f_hi) || f_hi < 0.0)
        return rejected_uniaxial(eps_zz, inelastic_eq, "strain_step_bracket_failed");

    int iterations = 0;
    for (; iterations < 80; ++iterations)
    {
        const double mid = 0.5 * (lo + hi);
        ArrheniusStatus mid_status;
        const double f_mid = residual(mid, &mid_status);
        if (!std::isfinite(f_mid))
            return rejected_uniaxial(eps_zz, inelastic_eq, "non_finite_residual");
        if (mid_status.strain_above_fit)
            return rejected_uniaxial(eps_zz, inelastic_eq, mid_status.diagnostic);
        if (f_mid > 0.0)
            hi = mid;
        else
            lo = mid;
    }
    const double p = 0.5 * (lo + hi);
    if (p > Ren2014X20Cr13::fit_strain_max)
        return rejected_uniaxial(eps_zz, inelastic_eq, "strain_would_exceed_0.75");

    UniaxialState state;
    state.eps_zz = eps_zz;
    state.inelastic_eq = p;
    state.eps_in_zz = -p;
    state.sigma_zz_Pa = youngs_modulus_Pa * (eps_zz + p);
    state.sigma_eq_Pa = compressive_sigma_eq_Pa(youngs_modulus_Pa, eps_zz, p);
    state.inelastic_rate_per_s =
        Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, state.sigma_eq_Pa, p, policy, true, &status);
    state.iterations = iterations;
    state.accepted = 1;
    state.diagnostic = status.diagnostic;
    return state;
}

UniaxialState uniaxial_constant_stress_step(double temperature_K, double youngs_modulus_Pa, double sigma_zz_Pa,
                                            double inelastic_eq, double dt_s, StartupPolicy policy)
{
    if (!(youngs_modulus_Pa > 0.0) || !(dt_s > 0.0) || !(inelastic_eq >= 0.0) || !(sigma_zz_Pa <= 0.0))
        return rejected_uniaxial(0.0, inelastic_eq, "invalid_stress_step");

    const double sigma_eq_Pa = -sigma_zz_Pa;
    auto residual = [&](double p, ArrheniusStatus *status) {
        return p - inelastic_eq -
               dt_s * Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_eq_Pa, p, policy, true, status);
    };

    ArrheniusStatus status;
    const double rate_start =
        Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_eq_Pa, inelastic_eq, policy, true, &status);
    if (!status.ok && status.strain_above_fit)
        return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq, status.diagnostic);
    if (!(rate_start > 0.0))
    {
        UniaxialState state;
        state.inelastic_eq = inelastic_eq;
        state.eps_in_zz = -inelastic_eq;
        state.eps_zz = sigma_zz_Pa / youngs_modulus_Pa + state.eps_in_zz;
        state.sigma_zz_Pa = sigma_zz_Pa;
        state.sigma_eq_Pa = sigma_eq_Pa;
        state.accepted = 1;
        state.diagnostic = status.diagnostic;
        return state;
    }

    double lo = inelastic_eq;
    double hi = inelastic_eq + dt_s * rate_start;
    int iterations = 0;
    for (int expand = 0; expand < 24; ++expand)
    {
        if (hi > Ren2014X20Cr13::fit_strain_max)
            return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq,
                                     "strain_would_exceed_0.75");
        ArrheniusStatus hi_status;
        const double f_hi = residual(hi, &hi_status);
        ++iterations;
        if (hi_status.strain_above_fit)
            return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq, hi_status.diagnostic);
        if (f_hi >= 0.0)
            break;
        const double rate_hi =
            Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_eq_Pa, hi, policy, true, &hi_status);
        const double next = inelastic_eq + dt_s * std::max(rate_hi, rate_start);
        if (!(next > hi))
            return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq, "stress_step_bracket_failed");
        hi = next;
        if (expand == 23)
        {
            ArrheniusStatus last_status;
            if (residual(hi, &last_status) < 0.0)
                return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq,
                                         "stress_step_bracket_failed");
        }
    }

    for (int it = 0; it < 80; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        ArrheniusStatus mid_status;
        const double f_mid = residual(mid, &mid_status);
        ++iterations;
        if (mid_status.strain_above_fit)
            return rejected_uniaxial(sigma_zz_Pa / youngs_modulus_Pa - inelastic_eq, inelastic_eq, mid_status.diagnostic);
        if (f_mid > 0.0)
            hi = mid;
        else
            lo = mid;
    }
    const double p = 0.5 * (lo + hi);
    UniaxialState state;
    state.inelastic_eq = p;
    state.eps_in_zz = -p;
    state.eps_zz = sigma_zz_Pa / youngs_modulus_Pa + state.eps_in_zz;
    state.sigma_zz_Pa = sigma_zz_Pa;
    state.sigma_eq_Pa = sigma_eq_Pa;
    state.inelastic_rate_per_s =
        Ren2014X20Cr13::inelastic_rate_Pa(temperature_K, sigma_eq_Pa, p, policy, true, &status);
    state.iterations = iterations;
    state.accepted = 1;
    state.diagnostic = status.diagnostic;
    return state;
}

namespace
{
struct CheckRow
{
    double strain;
    double alpha;
    double n;
    double Q;
    double lnA;
    double stress_MPa;
};

const CheckRow kChecks[] = {
    {0.05, 0.010774322203125, 6.6886268560375, 438336.5472125, 38.72763763065, 79.5777533},
    {0.10, 0.0094633175, 6.3739979602, 423411.7908, 37.3639699538, 89.2865360},
    {0.50, 0.0075646875, 5.357318125, 360364.85, 31.423978125, 106.4643916},
    {0.75, 0.007601826171875, 4.8660883984375, 331797.4671875, 28.805037890625, 101.5959056},
};

double rel_err(double got, double ref)
{
    if (ref == 0.0)
        return std::abs(got);
    return std::abs(got - ref) / std::abs(ref);
}

bool fail_if(bool bad, int &failures, const std::string &message)
{
    if (!bad)
        return false;
    std::cerr << "FAIL " << message << "\n";
    failures++;
    return true;
}
} // namespace

int run_ren2014_material_point()
{
    std::cout << "SPDB ren2014_x20cr13 material point\n"
              << "  DOI " << Ren2014X20Cr13::doi << "\n"
              << "  Literature hot-compression reference. Not a 950 C / 7200 s process prediction.\n"
              << "  Default startup is endpoint extension: coefficients use max(strain, 0.05), history starts at 0.\n";

    int failures = 0;
    std::ofstream csv("material_point_reference.csv");
    csv << std::setprecision(16);
    csv << "record,strain,T_K,strain_rate_per_s,stress_MPa,alpha_per_MPa,n,Q_J_per_mol,lnA,"
           "reference,relative_error,diagnostic\n";

    const double T_check = 1223.15;
    const double rate_check = 0.01;
    for (const CheckRow &row : kChecks)
    {
        const ArrheniusCoeffs c = Ren2014X20Cr13::coefficients(row.strain);
        ArrheniusStatus status;
        const double stress = Ren2014X20Cr13::flow_stress_MPa(T_check, rate_check, row.strain, StartupPolicy::PolyToZero,
                                                              false, &status);
        const double e_alpha = rel_err(c.alpha_per_MPa, row.alpha);
        const double e_n = rel_err(c.n, row.n);
        const double e_Q = rel_err(c.Q_J_per_mol, row.Q);
        const double e_lnA = rel_err(c.lnA, row.lnA);
        const double e_stress = rel_err(stress, row.stress_MPa);
        const double worst = std::max(e_alpha, std::max(e_n, std::max(e_Q, std::max(e_lnA, e_stress))));
        fail_if(worst > 1.0e-6, failures, "checkpoint strain " + std::to_string(row.strain) + " rel " + std::to_string(worst));
        std::cout << "checkpoint eps=" << row.strain << " sigma_MPa=" << std::setprecision(10) << stress
                  << " worst_rel=" << worst << "\n";
        csv << "checkpoint," << row.strain << "," << T_check << "," << rate_check << "," << stress << ","
            << c.alpha_per_MPa << "," << c.n << "," << c.Q_J_per_mol << "," << c.lnA << "," << row.stress_MPa << ","
            << e_stress << "," << status.diagnostic << "\n";
    }

    const double temperatures[] = {1173.0, 1223.0, 1423.0};
    const double strains[] = {0.05, 0.5, 0.75};
    const double rates[] = {0.01, 0.1, 1.0, 10.0};
    double worst_roundtrip = 0.0;
    double worst_pa = 0.0;
    for (double T : temperatures)
        for (double strain : strains)
            for (double rate : rates)
            {
                ArrheniusStatus forward_status;
                const double stress_MPa = Ren2014X20Cr13::flow_stress_MPa(T, rate, strain, StartupPolicy::PolyToZero, false,
                                                                          &forward_status);
                ArrheniusStatus back_status;
                const double rate_back = Ren2014X20Cr13::inelastic_rate(T, stress_MPa, strain, StartupPolicy::PolyToZero, false,
                                                                        &back_status);
                const double roundtrip = rel_err(rate_back, rate);
                worst_roundtrip = std::max(worst_roundtrip, roundtrip);
                const double stress_Pa = Ren2014X20Cr13::flow_stress_Pa(T, rate, strain, StartupPolicy::PolyToZero, false, nullptr);
                const double pa_mismatch = rel_err(stress_Pa, stress_MPa * 1.0e6);
                const double rate_from_pa =
                    Ren2014X20Cr13::inelastic_rate_Pa(T, stress_Pa, strain, StartupPolicy::PolyToZero, false, nullptr);
                const double rate_pa_mismatch = rel_err(rate_from_pa, rate_back);
                worst_pa = std::max(worst_pa, std::max(pa_mismatch, rate_pa_mismatch));
                csv << "roundtrip," << strain << "," << T << "," << rate << "," << stress_MPa << ",,,,," << rate_back << ","
                    << roundtrip << "," << forward_status.diagnostic << "\n";
            }
    fail_if(worst_roundtrip > 1.0e-6, failures, "roundtrip worst rel " + std::to_string(worst_roundtrip));
    fail_if(worst_pa > 1.0e-12, failures, "Pa/MPa mismatch " + std::to_string(worst_pa));
    std::cout << "roundtrip worst_rel=" << worst_roundtrip << " pa_mismatch=" << worst_pa << "\n";

    ArrheniusStatus zero_status;
    const double zero_rate =
        Ren2014X20Cr13::inelastic_rate(T_check, 0.0, 0.0, StartupPolicy::EndpointExtend, true, &zero_status);
    fail_if(zero_rate != 0.0, failures, "zero stress did not return zero rate");

    const double hydrostatic = Ren2014X20Cr13::von_mises_Pa(10.0e6, 10.0e6, 10.0e6, 0.0, 0.0, 0.0);
    const double uniaxial = Ren2014X20Cr13::von_mises_Pa(0.0, 0.0, -10.0e6, 0.0, 0.0, 0.0);
    fail_if(hydrostatic > 1.0e-8, failures, "hydrostatic von Mises is not zero");
    fail_if(rel_err(uniaxial, 10.0e6) > 1.0e-12, failures, "uniaxial von Mises");
    ArrheniusStatus hydro_rate_status;
    const double hydro_rate = Ren2014X20Cr13::inelastic_rate_Pa(T_check, hydrostatic, 0.0, StartupPolicy::EndpointExtend, true,
                                                                &hydro_rate_status);
    fail_if(hydro_rate != 0.0, failures, "hydrostatic equivalent did not give zero rate");

    const double ref_rate = 4.295261254e-9;
    const double ref_strain = 3.092588103e-5;
    ArrheniusStatus hold_status;
    const double rate_10 = Ren2014X20Cr13::inelastic_rate(T_check, 10.0, 0.0, StartupPolicy::EndpointExtend, true, &hold_status);
    fail_if(rel_err(rate_10, ref_rate) > 1.0e-6, failures, "10 MPa endpoint rate");
    std::cout << "hold10_rate=" << std::setprecision(12) << rate_10 << " diagnostic=" << hold_status.diagnostic
              << " low_rate=" << hold_status.low_rate_extrapolation << "\n";

    const double dts[] = {60.0, 30.0, 15.0};
    double previous_strain = 0.0;
    bool have_previous = false;
    for (double dt : dts)
    {
        const HoldIntegration hold =
            integrate_constant_stress(T_check, 10.0, 7200.0, dt, StartupPolicy::EndpointExtend, true);
        const double err = rel_err(hold.final_strain, ref_strain);
        fail_if(err > 1.0e-3, failures, "10 MPa hold dt " + std::to_string(dt));
        if (have_previous)
            fail_if(rel_err(hold.final_strain, previous_strain) > 1.0e-3, failures, "hold step halving");
        previous_strain = hold.final_strain;
        have_previous = true;
        std::cout << "hold10 dt=" << dt << " strain=" << hold.final_strain << " rel_to_constant_rate=" << err
                  << " low_rate_fraction=" << hold.low_rate_step_fraction << "\n";
        csv << "hold10_endpoint," << hold.final_strain << "," << T_check << "," << hold.final_rate_per_s << ",10,,,,"
            << ref_strain << "," << err << ",dt=" << dt << "\n";
    }

    const double extra_stresses[] = {20.0, 40.0};
    for (double stress : extra_stresses)
    {
        double previous = 0.0;
        bool have_step = false;
        for (double dt : dts)
        {
            const HoldIntegration hold =
                integrate_constant_stress(T_check, stress, 7200.0, dt, StartupPolicy::EndpointExtend, true);
            if (have_step && hold.final_strain != 0.0)
                fail_if(rel_err(hold.final_strain, previous) > 0.02, failures,
                        "stress " + std::to_string(stress) + " step halving");
            previous = hold.final_strain;
            have_step = true;
            std::cout << "sensitivity stress_MPa=" << stress << " dt=" << dt << " strain=" << hold.final_strain
                      << " stopped=" << hold.stopped << " reason=" << hold.stop_reason << "\n";
            csv << "hold_sensitivity," << hold.final_strain << "," << T_check << "," << hold.final_rate_per_s << ","
                << stress << ",,,,," << hold.stop_reason << ",dt=" << dt << "\n";
        }
    }

    const HoldIntegration poly =
        integrate_constant_stress(T_check, 10.0, 7200.0, 15.0, StartupPolicy::PolyToZero, true);
    std::cout << "poly_to_zero_10MPa strain=" << poly.final_strain << " rate=" << poly.final_rate_per_s << "\n";
    csv << "hold10_poly_to_zero," << poly.final_strain << "," << T_check << "," << poly.final_rate_per_s
        << ",10,,,,,sensitivity\n";

    // Demo elastic assumption. Not from Ren 2014. Inelastic increment does not depend on E.
    constexpr double E_demo_Pa = 1.0e11;
    constexpr double nu_demo = 0.30;
    constexpr double rho_demo = 7750.0;
    std::cout << "demo_assumption E_Pa=" << E_demo_Pa << " nu=" << nu_demo << " rho=" << rho_demo
              << " source=demo_assumption_not_ren2014\n";
    for (double scale : {0.5, 1.0, 2.0})
    {
        const double E = scale * E_demo_Pa;
        const double elastic_strain = -10.0e6 / E;
        const double inelastic_strain = previous_strain;
        std::cout << "E_scale=" << scale << " elastic_strain=" << elastic_strain
                  << " inelastic_strain=" << inelastic_strain
                  << " total_if_unloaded_residual=" << inelastic_strain << "\n";
        csv << "elastic_sensitivity," << inelastic_strain << "," << T_check << ",," << 10.0 << ",,,,," << elastic_strain
            << ",E_scale=" << scale << "\n";
    }

    // J2 direction in Cauchy space. Uniaxial compression: Dzz = -edot_eq, trace(Dvp) = 0.
    {
        const double szz = -40.0e6;
        const double seq = Ren2014X20Cr13::von_mises_Pa(0.0, 0.0, szz, 0.0, 0.0, 0.0);
        const double mean = szz / 3.0;
        const double edot = 1.0e-3;
        const double Dzz = 1.5 * (edot / seq) * (szz - mean);
        const double Dxx = 1.5 * (edot / seq) * (0.0 - mean);
        fail_if(rel_err(Dzz, -edot) > 1.0e-12, failures, "uniaxial J2 Dzz");
        fail_if(rel_err(Dxx, 0.5 * edot) > 1.0e-12, failures, "uniaxial J2 Dxx");
        fail_if(std::abs(2.0 * Dxx + Dzz) > 1.0e-12 * edot, failures, "J2 inelastic volume");
        std::cout << "j2_direction Dzz=" << Dzz << " Dxx=" << Dxx << "\n";
        csv << "j2_direction," << Dzz << "," << T_check << "," << edot << "," << (szz / 1.0e6) << ",,,,," << -edot << ","
            << rel_err(Dzz, -edot) << ",small_strain_cauchy\n";
    }

    // Device-callable J2 return: sigma_new = sigma_trial - 3 G rate(sigma_new) dt.
    {
        const double G = E_demo_Pa / (2.0 * (1.0 + nu_demo));
        const double sigma_trial = 80.0e6;
        const double dt_return = 0.05;
        const EquivalentReturn returned = equivalent_stress_backward_euler(
            sigma_trial, G, dt_return, T_check, 0.0, StartupPolicy::EndpointExtend);
        const double residual = returned.sigma_eq_Pa - sigma_trial + 3.0 * G * dt_return * returned.rate_per_s;
        fail_if(returned.accepted != 1, failures, "equivalent return not accepted");
        fail_if(std::abs(residual) / sigma_trial > 1.0e-6, failures, "equivalent return residual");
        fail_if(rel_err(returned.inelastic_eq, dt_return * returned.rate_per_s) > 1.0e-8, failures,
                "equivalent return strain");
        const EquivalentReturn zero_trial =
            equivalent_stress_backward_euler(0.0, G, dt_return, T_check, 0.0, StartupPolicy::EndpointExtend);
        fail_if(zero_trial.sigma_eq_Pa != 0.0 || zero_trial.rate_per_s != 0.0, failures, "zero trial return");
        std::cout << "equivalent_return sigma_MPa=" << (returned.sigma_eq_Pa / 1.0e6)
                  << " rate=" << returned.rate_per_s << " residual_rel=" << (residual / sigma_trial) << "\n";
        csv << "equivalent_return," << returned.inelastic_eq << "," << T_check << "," << returned.rate_per_s << ","
            << (returned.sigma_eq_Pa / 1.0e6) << ",,,,," << (sigma_trial / 1.0e6) << "," << (residual / sigma_trial)
            << ",j2_fixed_strain\n";
    }

    // Imposed axial rate -0.01. The first 1e-5 of strain is still elastic; at strain 0.2 the
    // stress should sit on the Arrhenius flow curve. E is demo_assumption.
    const double imposed_rate = -0.01;
    const double tiny_dt = 1.0e-3;
    const UniaxialState tiny = uniaxial_backward_euler(T_check, E_demo_Pa, imposed_rate * tiny_dt, 0.0, tiny_dt,
                                                       StartupPolicy::EndpointExtend);
    fail_if(tiny.accepted != 1, failures, "elastic transient step");
    fail_if(rel_err(tiny.sigma_eq_Pa, E_demo_Pa * std::abs(tiny.eps_zz)) > 1.0e-6, failures, "elastic transient stress");
    const double flow_window_Pa =
        Ren2014X20Cr13::flow_stress_Pa(T_check, 0.01, 0.05, StartupPolicy::EndpointExtend, false, nullptr);
    fail_if(!(tiny.sigma_eq_Pa < 0.2 * flow_window_Pa), failures, "elastic transient already at flow stress");
    std::cout << "elastic_transient sigma_eq_Pa=" << tiny.sigma_eq_Pa << " flow_window_Pa=" << flow_window_Pa << "\n";
    csv << "elastic_transient," << tiny.eps_zz << "," << T_check << "," << imposed_rate << "," << (tiny.sigma_eq_Pa / 1.0e6)
        << ",,,,," << (E_demo_Pa * std::abs(tiny.eps_zz) / 1.0e6) << ","
        << rel_err(tiny.sigma_eq_Pa, E_demo_Pa * std::abs(tiny.eps_zz)) << ",E_demo\n";

    struct Marched
    {
        UniaxialState final_state;
        int steps = 0;
        int stopped = 0;
    };
    auto march_rate = [&](double E, double dt, double duration) {
        Marched out;
        double t = 0.0;
        double eps = 0.0;
        double p = 0.0;
        while (t < duration - 0.5 * dt)
        {
            const double step = std::min(dt, duration - t);
            const UniaxialState trial = uniaxial_backward_euler(T_check, E, eps + imposed_rate * step, p, step,
                                                                StartupPolicy::EndpointExtend);
            if (trial.accepted != 1)
            {
                out.stopped = 1;
                out.final_state = trial;
                return out;
            }
            eps = trial.eps_zz;
            p = trial.inelastic_eq;
            t += step;
            out.steps++;
            out.final_state = trial;
        }
        return out;
    };

    const double rate_duration = 0.2 / 0.01;
    const Marched rate_coarse = march_rate(E_demo_Pa, 0.05, rate_duration);
    const Marched rate_fine = march_rate(E_demo_Pa, 0.025, rate_duration);
    fail_if(rate_coarse.stopped || rate_fine.stopped, failures, "constant strain-rate march stopped");
    const double flow_coarse = Ren2014X20Cr13::flow_stress_Pa(T_check, 0.01, rate_coarse.final_state.inelastic_eq,
                                                              StartupPolicy::EndpointExtend, false, nullptr);
    const double flow_fine = Ren2014X20Cr13::flow_stress_Pa(T_check, 0.01, rate_fine.final_state.inelastic_eq,
                                                            StartupPolicy::EndpointExtend, false, nullptr);
    fail_if(rel_err(rate_coarse.final_state.sigma_eq_Pa, flow_coarse) > 1.0e-3, failures, "strain-rate vs flow dt=0.05");
    fail_if(rel_err(rate_fine.final_state.sigma_eq_Pa, flow_fine) > 1.0e-3, failures, "strain-rate vs flow dt=0.025");
    fail_if(rel_err(rate_coarse.final_state.sigma_eq_Pa, rate_fine.final_state.sigma_eq_Pa) > 5.0e-3, failures,
            "strain-rate step halving");
    std::cout << "strain_rate dt=0.05 sigma_MPa=" << (rate_coarse.final_state.sigma_eq_Pa / 1.0e6)
              << " p=" << rate_coarse.final_state.inelastic_eq << " flow_MPa=" << (flow_coarse / 1.0e6)
              << " steps=" << rate_coarse.steps << "\n";
    std::cout << "strain_rate dt=0.025 sigma_MPa=" << (rate_fine.final_state.sigma_eq_Pa / 1.0e6)
              << " p=" << rate_fine.final_state.inelastic_eq << " flow_MPa=" << (flow_fine / 1.0e6) << "\n";
    csv << "strain_rate_flow," << rate_fine.final_state.inelastic_eq << "," << T_check << ",0.01,"
        << (rate_fine.final_state.sigma_eq_Pa / 1.0e6) << ",,,,," << (flow_fine / 1.0e6) << ","
        << rel_err(rate_fine.final_state.sigma_eq_Pa, flow_fine) << ",dt=0.025\n";

    for (double scale : {0.5, 1.0, 2.0})
    {
        const Marched scaled = march_rate(scale * E_demo_Pa, 0.05, rate_duration);
        fail_if(scaled.stopped, failures, "strain-rate E sensitivity stopped");
        fail_if(rel_err(scaled.final_state.sigma_eq_Pa, rate_coarse.final_state.sigma_eq_Pa) > 1.0e-2, failures,
                "strain-rate E sensitivity");
        csv << "strain_rate_E," << scaled.final_state.inelastic_eq << "," << T_check << ",0.01,"
            << (scaled.final_state.sigma_eq_Pa / 1.0e6) << ",,,,," << scale << ",,"
            << "E_scale=" << scale << "\n";
    }

    auto hold_at_stress = [&](double E, double sigma_zz_Pa, double dt_hold) {
        Marched out;
        double t = 0.0;
        double p = 0.0;
        while (t < 7200.0 - 0.5 * dt_hold)
        {
            const double step = std::min(dt_hold, 7200.0 - t);
            const UniaxialState trial =
                uniaxial_constant_stress_step(T_check, E, sigma_zz_Pa, p, step, StartupPolicy::EndpointExtend);
            if (trial.accepted != 1)
            {
                out.stopped = 1;
                out.final_state = trial;
                return out;
            }
            p = trial.inelastic_eq;
            t += step;
            out.steps++;
            out.final_state = trial;
        }
        return out;
    };

    const Marched loaded = hold_at_stress(E_demo_Pa, -10.0e6, 15.0);
    fail_if(loaded.stopped, failures, "10 MPa implicit hold stopped");
    fail_if(rel_err(loaded.final_state.inelastic_eq, ref_strain) > 1.0e-6, failures, "implicit hold vs constant rate");
    const UniaxialState unloaded =
        uniaxial_backward_euler(T_check, E_demo_Pa, -loaded.final_state.inelastic_eq, loaded.final_state.inelastic_eq, 60.0,
                                StartupPolicy::EndpointExtend);
    fail_if(unloaded.accepted != 1, failures, "unload step");
    fail_if(rel_err(unloaded.inelastic_eq, loaded.final_state.inelastic_eq) > 1.0e-12, failures, "unload changed inelastic strain");
    fail_if(unloaded.sigma_eq_Pa > 1.0e-8, failures, "unload left residual stress");
    fail_if(unloaded.inelastic_rate_per_s > 0.0, failures, "unload left a creep rate");
    const double recovered = loaded.final_state.eps_zz - unloaded.eps_zz;
    fail_if(rel_err(recovered, -10.0e6 / E_demo_Pa) > 1.0e-10, failures, "elastic recovery");
    std::cout << "load_hold_unload p=" << loaded.final_state.inelastic_eq
              << " eps_loaded=" << loaded.final_state.eps_zz << " eps_residual=" << unloaded.eps_zz
              << " recovered=" << recovered << "\n";
    csv << "load_hold_unload," << unloaded.eps_zz << "," << T_check << "," << unloaded.inelastic_rate_per_s << ",-10,,,,"
        << loaded.final_state.inelastic_eq << "," << rel_err(recovered, -10.0e6 / E_demo_Pa) << ",residual_is_inelastic\n";

    for (double scale : {0.5, 2.0})
    {
        const Marched scaled_hold = hold_at_stress(scale * E_demo_Pa, -10.0e6, 15.0);
        fail_if(scaled_hold.stopped, failures, "hold E sensitivity stopped");
        fail_if(rel_err(scaled_hold.final_state.inelastic_eq, loaded.final_state.inelastic_eq) > 1.0e-10, failures,
                "hold inelastic depends on E");
        const double elastic = -10.0e6 / (scale * E_demo_Pa);
        fail_if(rel_err(scaled_hold.final_state.eps_zz - scaled_hold.final_state.eps_in_zz, elastic) > 1.0e-10, failures,
                "hold elastic strain");
    }

    for (double strain = 0.05; strain <= 0.75 + 1.0e-12; strain += 0.05)
    {
        for (double rate : rates)
        {
            ArrheniusStatus status;
            const double stress =
                Ren2014X20Cr13::flow_stress_MPa(1223.0, rate, strain, StartupPolicy::PolyToZero, false, &status);
            const ArrheniusCoeffs c = Ren2014X20Cr13::coefficients(strain);
            csv << "flow_curve," << strain << ",1223," << rate << "," << stress << "," << c.alpha_per_MPa << "," << c.n
                << "," << c.Q_J_per_mol << "," << c.lnA << ",,," << status.diagnostic << "\n";
        }
    }

    std::cout << "material_point_failures=" << failures << "\n";
    csv.close();
    return failures == 0 ? 0 : 1;
}

} // namespace spdb
