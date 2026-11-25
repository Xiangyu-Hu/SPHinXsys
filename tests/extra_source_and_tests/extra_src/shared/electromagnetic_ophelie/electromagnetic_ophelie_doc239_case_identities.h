/**
 * @file electromagnetic_ophelie_doc239_case_identities.h
 * @brief Doc239/239A case identities and Jacoutot Ch.3 material/sigma adapters.
 *
 * Does not mutate CLIENT_* presets. French journal 50 kW and thesis 55 kW stay distinct.
 * Sigma A = printed formula; Sigma B = graph-fit with explicit EXTRAPOLATED zones.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC239_CASE_IDENTITIES_H
#define ELECTROMAGNETIC_OPHELIE_DOC239_CASE_IDENTITIES_H

#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kDoc239CaseClientNcBaseline = "CLIENT_NC_BASELINE";
inline constexpr const char *kDoc239CaseClientStirringBaseline = "CLIENT_STIRRING_BASELINE";
inline constexpr const char *kDoc239CaseFrenchPaperNcBaseline = "FRENCH_PAPER_NC_BASELINE";
inline constexpr const char *kDoc239CaseJacCh3NoMarangoni = "JAC_CH3_55KW_NO_MARANGONI";
inline constexpr const char *kDoc239CaseJacCh3WithMarangoni = "JAC_CH3_55KW_WITH_MARANGONI";

inline constexpr const char *kDoc239SigmaPrintedFormulaK = "JAC_PRINTED_FORMULA_K";
inline constexpr const char *kDoc239SigmaGraphFitExtrapolated = "JAC_GRAPH_FIT_EXTRAPOLATED";
inline constexpr const char *kDoc239SigmaGraphTraceOnly = "JAC_GRAPH_TRACE_ONLY";

/** Printed thesis formula (III.6 / C2), T in kelvin. */
inline Real doc239SigmaPrintedFormulaK(Real temperature_k)
{
    return std::pow(Real(10), Real(3.7921) - Real(3179.8) / temperature_k);
}

/**
 * Doc239 graph-fit candidate (NOT original author formula).
 * Fit range 1100–1550 K; outside = model extrapolation.
 */
inline Real doc239SigmaGraphFitExtrapolated(Real temperature_k)
{
    return std::pow(Real(10), Real(3.76180257442893) - Real(3869.111098088437) / temperature_k);
}

inline bool doc239SigmaBIsExtrapolated(Real temperature_k)
{
    return temperature_k < Real(1100) || temperature_k > Real(1550);
}

inline OphelieTemperatureLaw makeDoc239SigmaPrintedFormulaLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1550), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        sigma.setPoint(i, t, doc239SigmaPrintedFormulaK(t));
    }
    return sigma;
}

inline OphelieTemperatureLaw makeDoc239SigmaGraphFitExtrapolatedLaw()
{
    OphelieTemperatureLaw sigma;
    sigma.logarithmic_value_ = true;
    const Real temperatures[12] = {Real(900),  Real(1000), Real(1100), Real(1200), Real(1300), Real(1400),
                                   Real(1473), Real(1500), Real(1550), Real(1600), Real(1650), Real(1700)};
    for (size_t i = 0; i < 12; ++i)
    {
        const Real t = temperatures[i];
        sigma.setPoint(i, t, doc239SigmaGraphFitExtrapolated(t));
    }
    return sigma;
}

/**
 * Install Jacoutot Ch.3 thermal/rheology from Doc239 sparse solver CSVs
 * (docs/ophelie/reference/doc239_runs/solver_tables/).
 * Falls back to CES Fig.2 digitization if CSV missing (and reports fallback).
 */
inline ThesisMaterialCurveLoadReport installDoc239JacCh3ThermalRheologyFromSolverCsv(
    OphelieFrenchGlassMaterialLaws &laws, const std::string &solver_tables_dir)
{
    ThesisMaterialCurveLoadReport r;
    r.k_path = solver_tables_dir + "/solver_jac_III8_k_T.csv";
    r.cp_path = solver_tables_dir + "/solver_jac_III9_cp_T.csv";
    r.mu_path = solver_tables_dir + "/solver_jac_III10_mu_T.csv";
    r.endpoint_policy =
        "Doc239: endpoint hold outside knots; mu figure cold platform ~8e4 retained; "
        "text mu upper limit 1e5 is a separate config clamp if enabled";
    r.interpolation = "k,cp linear; mu log-interp";
    const bool k_ok = loadOphelieTemperatureLawFromCsv(laws.k, r.k_path, false, r.k_msg, true, Real(303));
    const bool cp_ok = loadOphelieTemperatureLawFromCsv(laws.cp, r.cp_path, false, r.cp_msg, true, Real(303));
    const bool mu_ok = loadOphelieTemperatureLawFromCsv(laws.mu, r.mu_path, true, r.mu_msg, true, Real(303));
    r.ok = k_ok && cp_ok && mu_ok;
    if (r.ok)
    {
        laws.fig2_thermal_rheology_digitized_ = false;
        laws.paper_digitized_ = true;
        laws.k_source_ = "Doc239 Jac III.8 solver table (from figure trace; not re-extracted)";
        laws.cp_source_ = "Doc239 Jac III.9 solver table (from figure trace; not re-extracted)";
        laws.mu_source_ =
            "Doc239 Jac III.10 solver table (figure platform ~8e4; text 1e5 separate constraint)";
        r.k_1473 = laws.k.evaluate(Real(1473));
        r.cp_1473 = laws.cp.evaluate(Real(1473));
        r.mu_1473 = laws.mu.evaluate(Real(1473));
    }
    return r;
}

struct Doc239CaseApplyReport
{
    bool ok = false;
    std::string case_identity;
    std::string sigma_identity;
    std::string material_note;
    std::string thermal_msg;
    bool marangoni_ready = false;
    Real target_power_w = 0;
    Real rho = 0;
    Real beta = 0;
    Real sigma_at_T0 = 0;
    Real T0 = 1473;
};

/**
 * Apply Doc239 French identities onto literature + french geometry knobs.
 * Does not touch client STL / client material paths.
 */
inline Doc239CaseApplyReport applyDoc239CaseIdentity(const std::string &case_identity,
                                                     const std::string &sigma_identity,
                                                     FrenchLiteratureParameters &lit,
                                                     OphelieFrenchReducedCaseParams &french,
                                                     const std::string &solver_tables_dir,
                                                     bool apply_mu_text_upper_1e5)
{
    Doc239CaseApplyReport rep;
    rep.case_identity = case_identity;
    rep.sigma_identity = sigma_identity;
    rep.T0 = Real(1473);

    if (case_identity == kDoc239CaseJacCh3WithMarangoni)
    {
        rep.ok = false;
        rep.marangoni_ready = false;
        rep.material_note = "JAC_CH3_55KW_WITH_MARANGONI marked NOT_READY (no trusted Marangoni implementation)";
        std::cerr << "[ophelie][doc239] FATAL: " << rep.material_note << "\n";
        return rep;
    }

    const bool is_jac_ch3 = (case_identity == kDoc239CaseJacCh3NoMarangoni);
    const bool is_paper = (case_identity == kDoc239CaseFrenchPaperNcBaseline);

    if (!is_jac_ch3 && !is_paper)
    {
        rep.ok = false;
        rep.material_note = "unsupported case_identity for applyDoc239CaseIdentity: " + case_identity;
        std::cerr << "[ophelie][doc239] FATAL: " << rep.material_note << "\n";
        return rep;
    }

    // Preserve CLI coupling knobs already parsed into lit (--zero-q / --zero-buoyancy / ...).
    // Resetting literature must not silently re-enable Q or Boussinesq (Doc239 P3 OFF bug).
    // Also keep periodic refresh interval / sigma UR / em_control (Doc239 P6: interval=0.5 was
    // wiped back to default 10 s → periodic_refreshes=0 within end_time=2).
    const bool keep_zero_q = lit.zero_q;
    const bool keep_boussinesq = lit.enable_boussinesq;
    const bool keep_constant_sigma = lit.constant_sigma;
    const bool keep_thermal_diffusion = lit.enable_thermal_diffusion;
    const bool keep_constant_thermal_rheology = lit.constant_thermal_rheology;
    const ThermoEMCouplingMode keep_coupling = lit.coupling;
    const Real keep_em_update_interval_s = lit.em_update_interval_s;
    const Real keep_sigma_under_relaxation = lit.sigma_under_relaxation;
    const EMControlMode keep_em_control = lit.em_control;

    // Geometry shared French cylinder defaults already in frozen_q main; enforce Ch.3 numbers.
    french.glass_radius = Real(0.25);
    french.glass_half_height = Real(0.5) * Real(0.185);
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    french.coil.loop_radius = Real(0.285); // 570 mm diameter centerline approx; OPEN coil section
    french.coil.num_loops = 7;
    french.coil.z_min = Real(-0.0225);
    french.coil.z_max = Real(0.2075);
    french.auto_coil_z = false;
    french.frequency_hz = Real(282000);

    lit = makeFrenchNaturalLiteratureParameters();
    lit.frequency_hz = Real(282000);
    lit.t_initial_k = rep.T0;
    lit.zero_q = keep_zero_q;
    lit.enable_boussinesq = keep_boussinesq;
    lit.constant_sigma = keep_constant_sigma;
    lit.enable_thermal_diffusion = keep_thermal_diffusion;
    lit.constant_thermal_rheology = keep_constant_thermal_rheology;
    lit.coupling = keep_coupling;
    lit.em_update_interval_s = keep_em_update_interval_s;
    lit.sigma_under_relaxation = keep_sigma_under_relaxation;
    lit.em_control = keep_em_control;

    if (is_jac_ch3)
    {
        lit.target_glass_absorbed_power_w = Real(55000);
        french.target_joule_power = Real(55000);
        // Thesis Ch.3 constants
        lit.material.rho.clear();
        lit.material.beta.clear();
        lit.material.rho.setPoint(0, rep.T0, Real(2800));
        lit.material.beta.setPoint(0, rep.T0, Real(3.0e-5));
        auto thermal = installDoc239JacCh3ThermalRheologyFromSolverCsv(lit.material, solver_tables_dir);
        rep.thermal_msg = thermal.ok ? ("thermal_csv_ok " + thermal.mu_msg)
                                     : ("thermal_csv_FAIL " + thermal.k_msg + "|" + thermal.cp_msg + "|" +
                                        thermal.mu_msg);
        if (!thermal.ok)
        {
            // Doc245: missing thesis solver CSV must FAIL — never silent CES Fig.2 fallback.
            rep.ok = false;
            rep.material_note =
                "FATAL: Doc239 Jac Ch.3 solver CSV missing under " + solver_tables_dir +
                " (refusing CES Fig.2 silent fallback). Restore docs/ophelie/reference/doc239_runs/solver_tables/";
            std::cerr << "[ophelie][doc239] " << rep.material_note << "\n"
                      << "[ophelie][doc239] " << rep.thermal_msg << "\n";
            return rep;
        }
        if (apply_mu_text_upper_1e5)
        {
            // Config clamp to text upper limit; figure platform evidence remains in raw/solver CSV.
            for (size_t i = 0; i < lit.material.mu.point_count_; ++i)
            {
                if (lit.material.mu.value_[i] > Real(1.0e5))
                    lit.material.mu.value_[i] = Real(1.0e5);
            }
            lit.material.mu_source_ =
                "Doc239 Jac III.10 + text upper clamp mu<=1e5 (figure platform evidence kept in CSV)";
        }
        rep.material_note = "JAC_CH3_55KW_NO_MARANGONI: rho2800 beta3e-5 H0.185 f282kHz P=55kW gamma=0";
        rep.marangoni_ready = false;
    }
    else
    {
        // Journal paper baseline identity (50 kW); keep existing CEP material unless sigma override.
        lit.target_glass_absorbed_power_w = Real(50000);
        french.target_joule_power = Real(50000);
        rep.material_note = "FRENCH_PAPER_NC_BASELINE: journal 50kW identity; do not mix with thesis 55kW metrics";
        rep.thermal_msg = "paper_baseline_uses_existing_CEP_material_unless_sigma_override";
    }

    // Sigma identity (A/B). First A/B round must only change sigma.
    if (sigma_identity.empty() || sigma_identity == "default")
    {
        if (is_jac_ch3)
            rep.sigma_identity = kDoc239SigmaGraphFitExtrapolated;
        else
            rep.sigma_identity = "CEP2008_journal_anchored_sigma16_shape";
    }

    if (rep.sigma_identity == kDoc239SigmaPrintedFormulaK || rep.sigma_identity == "A")
    {
        lit.material.sigma = makeDoc239SigmaPrintedFormulaLaw();
        lit.material.sigma_source_ = "Doc239 JAC_PRINTED_FORMULA_K: 10**(3.7921-3179.8/T)";
        lit.material.paper_digitized_ = true;
        rep.sigma_identity = kDoc239SigmaPrintedFormulaK;
    }
    else if (rep.sigma_identity == kDoc239SigmaGraphFitExtrapolated || rep.sigma_identity == "B" ||
             rep.sigma_identity == "B_fit")
    {
        lit.material.sigma = makeDoc239SigmaGraphFitExtrapolatedLaw();
        lit.material.sigma_source_ =
            "Doc239 JAC_GRAPH_FIT_EXTRAPOLATED: 10**(3.76180-3869.11/T); outside 1100-1550K = EXTRAPOLATED";
        lit.material.paper_digitized_ = true;
        rep.sigma_identity = kDoc239SigmaGraphFitExtrapolated;
    }
    else if (rep.sigma_identity == kDoc239SigmaGraphTraceOnly)
    {
        std::string msg;
        const std::string path = solver_tables_dir + "/solver_jac_III12_sigma_trace_T.csv";
        if (!loadOphelieTemperatureLawFromCsv(lit.material.sigma, path, true, msg, false))
        {
            rep.ok = false;
            rep.material_note = "cannot load graph trace sigma: " + msg;
            std::cerr << "[ophelie][doc239] FATAL: " << rep.material_note << "\n";
            return rep;
        }
        lit.material.sigma_source_ = "Doc239 JAC_GRAPH_TRACE_ONLY (readable figure range; no low-T invent)";
        lit.material.paper_digitized_ = true;
    }
    else if (is_paper && (rep.sigma_identity == "CEP2008_journal_anchored_sigma16_shape" ||
                          rep.sigma_identity == "default"))
    {
        // keep makeFrenchNaturalLiteratureParameters sigma
        rep.sigma_identity = "CEP2008_journal_anchored_sigma16_shape";
    }
    else
    {
        rep.ok = false;
        rep.material_note = "unknown sigma_identity: " + rep.sigma_identity;
        std::cerr << "[ophelie][doc239] FATAL: " << rep.material_note << "\n";
        return rep;
    }

    french.sigma_glass = lit.material.sigma.evaluate(rep.T0);
    french.frequency_hz = lit.frequency_hz;
    rep.rho = lit.material.rho.evaluate(rep.T0);
    rep.beta = lit.material.beta.evaluate(rep.T0);
    rep.sigma_at_T0 = french.sigma_glass;
    rep.target_power_w = lit.target_glass_absorbed_power_w;
    rep.ok = true;

    std::cout << "[ophelie][doc239] case_identity=" << rep.case_identity
              << " sigma_identity=" << rep.sigma_identity << " T0=" << rep.T0
              << " sigma(T0)=" << rep.sigma_at_T0 << " rho=" << rep.rho << " beta=" << rep.beta
              << " P_target=" << rep.target_power_w
              << " marangoni_ready=" << (rep.marangoni_ready ? 1 : 0)
              << " zero_q=" << (lit.zero_q ? 1 : 0)
              << " boussinesq=" << (lit.enable_boussinesq ? 1 : 0) << "\n"
              << "[ophelie][doc239] material_note=" << rep.material_note << "\n"
              << "[ophelie][doc239] thermal=" << rep.thermal_msg << "\n";
    if (rep.sigma_identity == kDoc239SigmaGraphFitExtrapolated)
    {
        std::cout << "[ophelie][doc239] sigma_B status_at_T0="
                  << (doc239SigmaBIsExtrapolated(rep.T0) ? "EXTRAPOLATED" : "FIT_IN_RANGE") << "\n";
    }
    return rep;
}

/**
 * Doc239 thermal BC packet (no include of OphelieThermalDiffusionOneWayOptions).
 * Jac Ch.3 (正文62–67): h_side=h_bottom=300, h_free=20, ε=0.9, Tcool=Tamb=303, Trad=423.
 * Journal paper baseline keeps prior production mix (h_bottom=35, ε=0.8, T*=300).
 * Bare-wall at T0=1473 still loses ≫55 kW; skull/crust series is a separate OPEN item.
 */
struct Doc239ThermalBcPacket
{
    bool set = false;
    Real h_side = 0;
    Real h_bottom = 0;
    Real h_free = 0;
    Real emissivity = 0;
    Real t_cool = 0;
    Real t_ambient = 0;
    Real t_rad_ambient = 0;
    std::string note;
};

inline Doc239ThermalBcPacket doc239ThermalBcForCase(const std::string &case_identity)
{
    Doc239ThermalBcPacket p;
    if (case_identity == kDoc239CaseJacCh3NoMarangoni)
    {
        p.set = true;
        p.h_side = Real(300);
        p.h_bottom = Real(300);
        p.h_free = Real(20);
        p.emissivity = Real(0.9);
        p.t_cool = Real(303);
        p.t_ambient = Real(303);
        p.t_rad_ambient = Real(423);
        p.note = "JAC_CH3 thesis BC h_side=h_bottom=300 h_free=20 eps=0.9 Tcool=Tamb=303 Trad=423 "
                 "(bare-wall T0 loss still >>55kW; skull/crust OPEN)";
        return p;
    }
    if (case_identity == kDoc239CaseFrenchPaperNcBaseline)
    {
        p.set = true;
        p.h_side = Real(300);
        p.h_bottom = Real(35);
        p.h_free = Real(20);
        p.emissivity = Real(0.8);
        p.t_cool = Real(300);
        p.t_ambient = Real(300);
        p.t_rad_ambient = Real(300);
        p.note = "FRENCH_PAPER_NC_BASELINE production BC h_bottom=35 eps=0.8 T*=300 (do not mix with Jac Ch.3)";
        return p;
    }
    return p;
}

inline bool isDoc239CliOption(const char *arg)
{
    return std::strncmp(arg, "--doc239-case-identity=", 23) == 0 ||
           std::strncmp(arg, "--doc239-sigma-identity=", 24) == 0 ||
           std::strncmp(arg, "--doc239-solver-tables-dir=", 27) == 0 ||
           std::strcmp(arg, "--doc239-mu-text-upper-1e5") == 0;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC239_CASE_IDENTITIES_H
