/**
 * @file electromagnetic_ophelie_client_material_laws.h
 * @brief Client melt material laws after Doc229 reply (2026-11-20).
 *
 * Viscosity (Pa·s):  mu(T) = exp(-4 + 8500/T), T in K, window 1050–1800 K.
 * Conductivity:      sigma_client(T) = exp(4.5 - 8500/T)  [natural exp; NOT pow(10,...)].
 *
 * French Table-1 constants remain separate material_ids for borrowed/fallback presets.
 * Do not mutate French regression IDs.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_CLIENT_MATERIAL_LAWS_H
#define ELECTROMAGNETIC_OPHELIE_CLIENT_MATERIAL_LAWS_H

#include "electromagnetic_ophelie_french_material_laws.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Stable IDs for configs / effective-case.json (do not rename silently). */
inline constexpr const char *kClientMaterialIdMuSigmaProvisional = "client_ln_sigma_provisional";
inline constexpr const char *kClientMaterialIdFrenchSigma16Borrowed = "french_table1_sigma16_borrowed";
inline constexpr const char *kClientMaterialIdFrenchSigmaTVerified = "french_sigma_T_verified_source";

inline constexpr Real kClientTWindowMinK = Real(1050);
inline constexpr Real kClientTWindowMaxK = Real(1800);
inline constexpr Real kClientNcProbeTargetK = Real(1273.15); // ~1000 °C thermocouple target
inline constexpr Real kClientStirPeakTrendK = Real(1673.15); // ~1400 °C qualitative trend only

inline Real clientDynamicViscosityPaS(Real temperature_k)
{
    return std::exp(Real(-4) + Real(8500) / temperature_k);
}

/** Client conductivity with natural exponential (Doc229 reply figure). */
inline Real clientConductivityLnExpSm(Real temperature_k)
{
    return std::exp(Real(4.5) - Real(8500) / temperature_k);
}

/** Equivalent log10 form of the same law (for logging only). */
inline Real clientConductivityLog10Equivalent(Real temperature_k)
{
    return (Real(4.5) - Real(8500) / temperature_k) / std::log(Real(10));
}

inline bool clientTemperatureInDeclaredWindow(Real temperature_k)
{
    return temperature_k >= kClientTWindowMinK && temperature_k <= kClientTWindowMaxK;
}

struct OphelieClientFormulaCheckRow
{
    Real t_k = 0;
    Real t_c = 0;
    Real mu = 0;
    Real sigma = 0;
    Real mu_ref = 0;
    Real sigma_ref = 0;
    Real mu_rel = 0;
    Real sigma_rel = 0;
    bool ok = false;
};

/**
 * Doc229 §2 / §7.3 recompute table.
 * Targets: 1473 K → mu≈5.87400286, sigma≈0.28068105 (relative tol 1e-4).
 */
inline std::vector<OphelieClientFormulaCheckRow> buildClientFormulaCheckRows()
{
    struct Node
    {
        Real t_k;
        Real mu_ref;
        Real sigma_ref;
    };
    const Node nodes[] = {
        {Real(1050), Real(60.05363594), Real(0.02745415)},
        {kClientNcProbeTargetK, Real(14.53201178), Real(0.11345444)},
        {Real(1473), Real(5.87400286), Real(0.28068105)},
        {kClientStirPeakTrendK, Real(2.94538010), Real(0.55976520)},
        {Real(1800), Real(2.05900369), Real(0.80073740)},
    };

    std::vector<OphelieClientFormulaCheckRow> rows;
    rows.reserve(5);
    constexpr Real tol = Real(1.0e-4);
    for (const Node &n : nodes)
    {
        OphelieClientFormulaCheckRow r;
        r.t_k = n.t_k;
        r.t_c = n.t_k - Real(273.15);
        r.mu = clientDynamicViscosityPaS(n.t_k);
        r.sigma = clientConductivityLnExpSm(n.t_k);
        r.mu_ref = n.mu_ref;
        r.sigma_ref = n.sigma_ref;
        r.mu_rel = std::abs(r.mu - r.mu_ref) / (std::abs(r.mu_ref) + TinyReal);
        r.sigma_rel = std::abs(r.sigma - r.sigma_ref) / (std::abs(r.sigma_ref) + TinyReal);
        r.ok = (r.mu_rel <= tol) && (r.sigma_rel <= tol);
        rows.push_back(r);
    }
    return rows;
}

inline bool writeClientFormulaCheckCsv(const std::string &path, std::string *error = nullptr)
{
    const auto rows = buildClientFormulaCheckRows();
    std::ofstream out(path);
    if (!out)
    {
        if (error)
            *error = "cannot open " + path;
        return false;
    }
    out << "T_K,T_C,mu_Pa_s,sigma_S_per_m,mu_ref,sigma_ref,mu_rel_err,sigma_rel_err,ok\n";
    out << std::setprecision(10);
    bool all_ok = true;
    for (const auto &r : rows)
    {
        all_ok = all_ok && r.ok;
        out << r.t_k << ',' << r.t_c << ',' << r.mu << ',' << r.sigma << ',' << r.mu_ref << ',' << r.sigma_ref
            << ',' << r.mu_rel << ',' << r.sigma_rel << ',' << (r.ok ? 1 : 0) << '\n';
    }
    // Consistency probe: mu*sigma = exp(0.5) for the shared Arrhenius pair.
    const Real t = Real(1473);
    const Real prod = clientDynamicViscosityPaS(t) * clientConductivityLnExpSm(t);
    const Real prod_ref = std::exp(Real(0.5));
    out << "# mu*sigma_at_1473=" << prod << " exp(0.5)=" << prod_ref
        << " rel=" << std::abs(prod - prod_ref) / prod_ref << " (consistency only)\n";
    return all_ok;
}

inline OphelieTemperatureLaw makeClientViscosityTemperatureLaw()
{
    OphelieTemperatureLaw mu;
    mu.logarithmic_value_ = true;
    const Real t[] = {Real(1050), Real(1200), Real(1273.15), Real(1400), Real(1473), Real(1600),
                      Real(1673.15), Real(1800)};
    for (size_t i = 0; i < 8; ++i)
        mu.setPoint(i, t[i], clientDynamicViscosityPaS(t[i]));
    return mu;
}

inline OphelieTemperatureLaw makeClientLnExpConductivityTemperatureLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    const Real t[] = {Real(1050), Real(1200), Real(1273.15), Real(1400), Real(1473), Real(1600),
                      Real(1673.15), Real(1800)};
    for (size_t i = 0; i < 8; ++i)
        sigma.setPoint(i, t[i], clientConductivityLnExpSm(t[i]));
    return sigma;
}

/**
 * Client provisional melt laws:
 * - mu/sigma from Doc229 corrected formulas
 * - rho=2700 provisional (2.7 g/cm³ interpretation; unit pending client confirm)
 * - cp/k/beta borrowed from French Table-1
 */
inline OphelieFrenchGlassMaterialLaws makeClientProvisionalMuSigmaMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws;
    laws.rho.setPoint(0, kClientNcProbeTargetK, Real(2700));
    laws.beta.setPoint(0, kClientNcProbeTargetK, Real(1.0e-5));
    laws.mu = makeClientViscosityTemperatureLaw();
    laws.cp.setPoint(0, kClientNcProbeTargetK, Real(1150));
    laws.k.setPoint(0, kClientNcProbeTargetK, Real(4));
    laws.sigma = makeClientLnExpConductivityTemperatureLaw();
    laws.paper_digitized_ = false;
    laws.fig2_thermal_rheology_digitized_ = false;
    laws.mu_source_ = "client Doc229: mu=exp(-4+8500/T), window 1050-1800 K";
    laws.cp_source_ = "borrowed French Table-1 cp=1150 (authorized reference; not client measurement)";
    laws.k_source_ = "borrowed French Table-1 k=4 (authorized reference; not client measurement)";
    laws.sigma_source_ = "client Doc229 reply figure: sigma=exp(4.5-8500/T) [ln]; provisional / may be wrong";
    return laws;
}

/** Development path: freeze mu at a reference T, keep client or French sigma as installed. */
inline void freezeClientDynamicViscosityAtTemperature(OphelieFrenchGlassMaterialLaws &laws, Real t_ref_k)
{
    const Real mu = clientDynamicViscosityPaS(t_ref_k);
    laws.mu.clear();
    laws.mu.setPoint(0, t_ref_k, mu);
    laws.mu.logarithmic_value_ = true;
    laws.fig2_thermal_rheology_digitized_ = false;
    laws.mu_source_ = "client frozen mu at reference T (NOT full mu(T) rheology)";
}

/**
 * Borrowed French Table-1 sigma=16 path for client geometry development.
 * Density still client provisional 2700; mu still client formula (caller may freeze).
 */
inline OphelieFrenchGlassMaterialLaws makeClientGeometryFrenchSigma16BorrowedMaterialLaws()
{
    OphelieFrenchGlassMaterialLaws laws = makeClientProvisionalMuSigmaMaterialLaws();
    laws.sigma.clear();
    laws.sigma.setPoint(0, Real(1473), Real(16));
    laws.sigma.logarithmic_value_ = true;
    laws.sigma_source_ = "french_table1_sigma16_borrowed on client geometry (not client measurement)";
    return laws;
}

inline void logClientMaterialIdentity(const char *material_id, const OphelieFrenchGlassMaterialLaws &laws,
                                      Real t_probe_k)
{
    std::cout << "[ophelie][client-mat] material_id=" << material_id
              << " T_probe_K=" << t_probe_k
              << " mu=" << laws.mu.evaluate(t_probe_k)
              << " sigma=" << laws.sigma.evaluate(t_probe_k)
              << " rho=" << laws.rho.evaluate(t_probe_k)
              << " cp=" << laws.cp.evaluate(t_probe_k)
              << " k=" << laws.k.evaluate(t_probe_k)
              << " beta=" << laws.beta.evaluate(t_probe_k) << std::endl;
    std::cout << "[ophelie][client-mat] sources"
              << " mu=[" << laws.mu_source_ << "]"
              << " sigma=[" << laws.sigma_source_ << "]"
              << " cp=[" << laws.cp_source_ << "]"
              << " k=[" << laws.k_source_ << "]" << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLIENT_MATERIAL_LAWS_H
