/**
 * @file electromagnetic_ophelie_doc245_material_dataset.h
 * @brief Doc245 independent material-dataset switch (μ/k/cp only).
 *
 * Decouples property tables from case identity (50 kW / β / ρ / BC / Q).
 * Missing thesis CSVs must FAIL — never silent CES Fig.2 fallback on JAC_THESIS_*.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC245_MATERIAL_DATASET_H
#define ELECTROMAGNETIC_OPHELIE_DOC245_MATERIAL_DATASET_H

#include "electromagnetic_ophelie_french_material_laws.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kDoc245MaterialLegacyApprox = "LEGACY_APPROX_TABLE";
inline constexpr const char *kDoc245MaterialJacThesisRaw = "JAC_THESIS_RAW_BORROWED";
inline constexpr const char *kDoc245MaterialJacThesisCepAnchored = "JAC_THESIS_SHAPE_CEP_ANCHORED";

struct Doc245MaterialApplyReport
{
    bool ok = false;
    std::string dataset_id;
    std::string tables_dir;
    std::string note;
    std::string k_msg;
    std::string cp_msg;
    std::string mu_msg;
    std::string k_path;
    std::string cp_path;
    std::string mu_path;
    int k_nodes = 0;
    int cp_nodes = 0;
    int mu_nodes = 0;
    Real k_1473 = 0;
    Real cp_1473 = 0;
    Real mu_1473 = 0;
    /** Snapshot of non-material knobs that must remain unchanged by this switch. */
    Real rho_1473_before = 0;
    Real beta_1473_before = 0;
    Real sigma_1473_before = 0;
    Real rho_1473_after = 0;
    Real beta_1473_after = 0;
    Real sigma_1473_after = 0;
    Real power_w_before = 0;
    Real power_w_after = 0;
};

inline std::string resolveDoc245SolverTablesDir(const std::string &requested)
{
    const std::vector<std::string> candidates = {
        requested,
        std::string("docs/ophelie/reference/doc239_runs/solver_tables"),
        std::string("../../../../../../docs/ophelie/reference/doc239_runs/solver_tables"),
        std::string("/home/ubuntu/SPHinXsys/docs/ophelie/reference/doc239_runs/solver_tables"),
    };
    for (const auto &c : candidates)
    {
        std::ifstream probe(c + "/solver_jac_III8_k_T.csv");
        if (probe.good())
            return c;
    }
    return requested;
}

/**
 * Apply material dataset to laws only. Does not touch rho/beta/sigma/power/geometry.
 * For JAC_THESIS_*: missing CSV → ok=false (no CES fallback).
 */
inline Doc245MaterialApplyReport applyDoc245MaterialDataset(const std::string &dataset_id,
                                                            OphelieFrenchGlassMaterialLaws &laws,
                                                            Real target_power_w,
                                                            const std::string &tables_dir_in,
                                                            bool apply_cep1473_anchors = false)
{
    Doc245MaterialApplyReport rep;
    rep.dataset_id = dataset_id.empty() ? kDoc245MaterialLegacyApprox : dataset_id;
    rep.power_w_before = target_power_w;
    rep.power_w_after = target_power_w;
    rep.rho_1473_before = laws.rho.evaluate(Real(1473));
    rep.beta_1473_before = laws.beta.evaluate(Real(1473));
    rep.sigma_1473_before = laws.sigma.evaluate(Real(1473));

    if (rep.dataset_id == kDoc245MaterialLegacyApprox || rep.dataset_id == "legacy" ||
        rep.dataset_id == "CES_FIG2")
    {
        installJacoutotCes2008Fig2ThermalRheology(laws);
        laws.fig2_thermal_rheology_digitized_ = true;
        laws.paper_digitized_ = false;
        laws.mu_source_ = "LEGACY_APPROX_TABLE: CES 2008 Fig.2 NUM visual digitization @Table-1 1473K";
        laws.cp_source_ = laws.mu_source_;
        laws.k_source_ = laws.mu_source_;
        rep.dataset_id = kDoc245MaterialLegacyApprox;
        rep.k_nodes = static_cast<int>(laws.k.point_count_);
        rep.cp_nodes = static_cast<int>(laws.cp.point_count_);
        rep.mu_nodes = static_cast<int>(laws.mu.point_count_);
        rep.k_1473 = laws.k.evaluate(Real(1473));
        rep.cp_1473 = laws.cp.evaluate(Real(1473));
        rep.mu_1473 = laws.mu.evaluate(Real(1473));
        rep.ok = true;
        rep.note = "legacy CES Fig.2 sparse approx; NOT Doc239 thesis digitization";
    }
    else if (rep.dataset_id == kDoc245MaterialJacThesisRaw ||
             rep.dataset_id == kDoc245MaterialJacThesisCepAnchored ||
             rep.dataset_id == "jac_thesis" || rep.dataset_id == "thesis_raw")
    {
        const bool want_anchor = (rep.dataset_id == kDoc245MaterialJacThesisCepAnchored) ||
                                 apply_cep1473_anchors;
        rep.dataset_id =
            want_anchor ? kDoc245MaterialJacThesisCepAnchored : kDoc245MaterialJacThesisRaw;
        rep.tables_dir = resolveDoc245SolverTablesDir(tables_dir_in);
        ThesisMaterialCurveLoadReport thermal;
        thermal.k_path = rep.tables_dir + "/solver_jac_III8_k_T.csv";
        thermal.cp_path = rep.tables_dir + "/solver_jac_III9_cp_T.csv";
        thermal.mu_path = rep.tables_dir + "/solver_jac_III10_mu_T.csv";
        rep.k_path = thermal.k_path;
        rep.cp_path = thermal.cp_path;
        rep.mu_path = thermal.mu_path;
        const bool k_ok =
            loadOphelieTemperatureLawFromCsv(laws.k, thermal.k_path, false, thermal.k_msg, true, Real(303));
        const bool cp_ok =
            loadOphelieTemperatureLawFromCsv(laws.cp, thermal.cp_path, false, thermal.cp_msg, true, Real(303));
        const bool mu_ok =
            loadOphelieTemperatureLawFromCsv(laws.mu, thermal.mu_path, true, thermal.mu_msg, true, Real(303));
        rep.k_msg = thermal.k_msg;
        rep.cp_msg = thermal.cp_msg;
        rep.mu_msg = thermal.mu_msg;
        if (!(k_ok && cp_ok && mu_ok))
        {
            rep.ok = false;
            rep.note = "FATAL: JAC thesis solver CSV missing/unreadable; refusing CES silent fallback. "
                       "Need solver_jac_III8_k_T.csv / III9_cp / III10_mu under " +
                       rep.tables_dir;
            std::cerr << "[ophelie][doc245][material] " << rep.note << "\n"
                      << "  k: " << rep.k_msg << "\n"
                      << "  cp: " << rep.cp_msg << "\n"
                      << "  mu: " << rep.mu_msg << "\n";
            return rep;
        }
        laws.fig2_thermal_rheology_digitized_ = false;
        laws.paper_digitized_ = true;
        if (want_anchor)
        {
            const Real k_raw = laws.k.evaluate(Real(1473));
            const Real cp_raw = laws.cp.evaluate(Real(1473));
            const Real mu_raw = laws.mu.evaluate(Real(1473));
            if (!(k_raw > TinyReal && cp_raw > TinyReal && mu_raw > TinyReal))
            {
                rep.ok = false;
                rep.note = "FATAL: CEP anchor failed (nonpositive thesis value at 1473 K)";
                return rep;
            }
            laws.k.scaleValues(Real(4) / k_raw);
            laws.cp.scaleValues(Real(1150) / cp_raw);
            laws.mu.scaleValues(Real(4) / mu_raw);
            laws.k_source_ =
                "JAC_THESIS_SHAPE_CEP_ANCHORED: Doc239 III.8 shape × CEP1473 k=4 (NOT raw thesis)";
            laws.cp_source_ =
                "JAC_THESIS_SHAPE_CEP_ANCHORED: Doc239 III.9 shape × CEP1473 cp=1150 (NOT raw thesis)";
            laws.mu_source_ =
                "JAC_THESIS_SHAPE_CEP_ANCHORED: Doc239 III.10 shape × CEP1473 mu=4 (NOT raw thesis)";
            rep.note = "shape retained; values scaled to journal Table-1 anchors at 1473 K — "
                       "do not call this raw thesis";
        }
        else
        {
            laws.k_source_ =
                "JAC_THESIS_RAW_BORROWED: Doc239 Jac III.8 solver table (figure trace downsample)";
            laws.cp_source_ =
                "JAC_THESIS_RAW_BORROWED: Doc239 Jac III.9 solver table (figure trace downsample)";
            laws.mu_source_ =
                "JAC_THESIS_RAW_BORROWED: Doc239 Jac III.10 solver table (log-interp; cold platform ~8e4)";
            rep.note = "borrowed thesis figure-derived curves for journal 50 kW NC; "
                       "rho/beta/sigma/power unchanged by this switch";
        }
        rep.k_nodes = static_cast<int>(laws.k.point_count_);
        rep.cp_nodes = static_cast<int>(laws.cp.point_count_);
        rep.mu_nodes = static_cast<int>(laws.mu.point_count_);
        rep.k_1473 = laws.k.evaluate(Real(1473));
        rep.cp_1473 = laws.cp.evaluate(Real(1473));
        rep.mu_1473 = laws.mu.evaluate(Real(1473));
        rep.ok = true;
    }
    else
    {
        rep.ok = false;
        rep.note = "unknown material_dataset: " + rep.dataset_id +
                   " (expected LEGACY_APPROX_TABLE | JAC_THESIS_RAW_BORROWED | "
                   "JAC_THESIS_SHAPE_CEP_ANCHORED)";
        std::cerr << "[ophelie][doc245][material] FATAL: " << rep.note << "\n";
        return rep;
    }

    rep.rho_1473_after = laws.rho.evaluate(Real(1473));
    rep.beta_1473_after = laws.beta.evaluate(Real(1473));
    rep.sigma_1473_after = laws.sigma.evaluate(Real(1473));
    const Real eps = Real(1.0e-12);
    if (std::abs(rep.rho_1473_after - rep.rho_1473_before) > eps ||
        std::abs(rep.beta_1473_after - rep.beta_1473_before) > eps ||
        std::abs(rep.sigma_1473_after - rep.sigma_1473_before) > eps ||
        std::abs(rep.power_w_after - rep.power_w_before) > eps)
    {
        rep.ok = false;
        rep.note = "FATAL: material_dataset mutated rho/beta/sigma/power (decoupling violated)";
        std::cerr << "[ophelie][doc245][material] " << rep.note << "\n";
        return rep;
    }

    std::cout << "[ophelie][doc245][material] dataset=" << rep.dataset_id
              << " tables_dir=" << (rep.tables_dir.empty() ? "(in-code)" : rep.tables_dir)
              << " k_nodes=" << rep.k_nodes << " cp_nodes=" << rep.cp_nodes
              << " mu_nodes=" << rep.mu_nodes << " @1473K: k=" << rep.k_1473 << " cp=" << rep.cp_1473
              << " mu=" << rep.mu_1473 << " rho=" << rep.rho_1473_after << " beta=" << rep.beta_1473_after
              << " sigma=" << rep.sigma_1473_after << " P=" << rep.power_w_after << "\n"
              << "[ophelie][doc245][material] note=" << rep.note << "\n"
              << "[ophelie][doc245][material] sources mu=[" << laws.mu_source_ << "] cp=["
              << laws.cp_source_ << "] k=[" << laws.k_source_ << "]\n";
    return rep;
}

inline bool writeDoc245MaterialProbeCsv(const std::string &path, const std::string &dataset_id,
                                        const OphelieFrenchGlassMaterialLaws &laws,
                                        const OphelieFrenchGlassMaterialLaws *legacy_for_compare = nullptr)
{
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "dataset_id,property,T_K,table_interp_value,actual_return_value,unit,handling,"
           "clamped_or_extrapolated,legacy_CES_Fig2_value\n";
    const Real probes[] = {Real(750),  Real(1000), Real(1200), Real(1400), Real(1473),
                           Real(1500), Real(303),  Real(1800)};
    struct Prop
    {
        const char *name;
        const char *unit;
        const OphelieTemperatureLaw *law;
        const OphelieTemperatureLaw *legacy;
        bool log_interp;
    };
    const OphelieTemperatureLaw *leg_mu = legacy_for_compare ? &legacy_for_compare->mu : nullptr;
    const OphelieTemperatureLaw *leg_cp = legacy_for_compare ? &legacy_for_compare->cp : nullptr;
    const OphelieTemperatureLaw *leg_k = legacy_for_compare ? &legacy_for_compare->k : nullptr;
    const Prop props[] = {
        {"mu", "Pa_s", &laws.mu, leg_mu, true},
        {"cp", "J_per_kg_K", &laws.cp, leg_cp, false},
        {"k", "W_per_m_K", &laws.k, leg_k, false},
        {"sigma", "S_per_m", &laws.sigma, nullptr, true},
    };
    for (const auto &p : props)
    {
        if (p.law->point_count_ == 0)
            continue;
        const Real t_lo = p.law->temperature_k_[0];
        const Real t_hi = p.law->temperature_k_[p.law->point_count_ - 1];
        for (Real T : probes)
        {
            const Real v = p.law->evaluate(T);
            const bool outside = (T < t_lo - Real(1e-9)) || (T > t_hi + Real(1e-9));
            const char *handling = outside ? "endpoint_hold" : (p.log_interp ? "log_linear_interp" : "linear_interp");
            out << dataset_id << "," << p.name << "," << T << "," << v << "," << v << "," << p.unit << ","
                << handling << "," << (outside ? 1 : 0) << ",";
            if (p.legacy && p.legacy->point_count_ > 0)
                out << p.legacy->evaluate(T);
            else
                out << "";
            out << "\n";
        }
        // Always emit table min/max nodes
        out << dataset_id << "," << p.name << "," << t_lo << "," << p.law->value_[0] << ","
            << p.law->evaluate(t_lo) << "," << p.unit << ",table_min_node,0,";
        if (p.legacy && p.legacy->point_count_ > 0)
            out << p.legacy->evaluate(t_lo);
        out << "\n";
        out << dataset_id << "," << p.name << "," << t_hi << "," << p.law->value_[p.law->point_count_ - 1]
            << "," << p.law->evaluate(t_hi) << "," << p.unit << ",table_max_node,0,";
        if (p.legacy && p.legacy->point_count_ > 0)
            out << p.legacy->evaluate(t_hi);
        out << "\n";
    }
    return true;
}

inline bool writeDoc245MaterialEffectiveJson(const std::string &path, const Doc245MaterialApplyReport &rep,
                                             const OphelieFrenchGlassMaterialLaws &laws, Real t_ref_k)
{
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "{\n";
    out << "  \"dataset_id\": \"" << rep.dataset_id << "\",\n";
    out << "  \"tables_dir\": \"" << rep.tables_dir << "\",\n";
    out << "  \"note\": \"" << rep.note << "\",\n";
    out << "  \"ok\": " << (rep.ok ? "true" : "false") << ",\n";
    out << "  \"paths\": {\n";
    out << "    \"k\": \"" << rep.k_path << "\",\n";
    out << "    \"cp\": \"" << rep.cp_path << "\",\n";
    out << "    \"mu\": \"" << rep.mu_path << "\"\n";
    out << "  },\n";
    out << "  \"load_messages\": {\n";
    out << "    \"k\": \"" << rep.k_msg << "\",\n";
    out << "    \"cp\": \"" << rep.cp_msg << "\",\n";
    out << "    \"mu\": \"" << rep.mu_msg << "\"\n";
    out << "  },\n";
    out << "  \"nodes\": {\"k\": " << rep.k_nodes << ", \"cp\": " << rep.cp_nodes
        << ", \"mu\": " << rep.mu_nodes << "},\n";
    out << "  \"at_1473K\": {\"k\": " << rep.k_1473 << ", \"cp\": " << rep.cp_1473
        << ", \"mu\": " << rep.mu_1473 << "},\n";
    out << "  \"unchanged_knobs\": {\n";
    out << "    \"rho_before\": " << rep.rho_1473_before << ", \"rho_after\": " << rep.rho_1473_after
        << ",\n";
    out << "    \"beta_before\": " << rep.beta_1473_before << ", \"beta_after\": " << rep.beta_1473_after
        << ",\n";
    out << "    \"sigma_before\": " << rep.sigma_1473_before
        << ", \"sigma_after\": " << rep.sigma_1473_after << ",\n";
    out << "    \"power_W_before\": " << rep.power_w_before << ", \"power_W_after\": " << rep.power_w_after
        << "\n";
    out << "  },\n";
    out << "  \"sources\": {\n";
    out << "    \"mu\": \"" << laws.mu_source_ << "\",\n";
    out << "    \"cp\": \"" << laws.cp_source_ << "\",\n";
    out << "    \"k\": \"" << laws.k_source_ << "\",\n";
    out << "    \"sigma\": \"" << laws.sigma_source_ << "\"\n";
    out << "  },\n";
    out << "  \"evaluate_at_T_ref_K\": " << t_ref_k << ",\n";
    out << "  \"evaluate\": {\n";
    out << "    \"mu\": " << laws.mu.evaluate(t_ref_k) << ",\n";
    out << "    \"cp\": " << laws.cp.evaluate(t_ref_k) << ",\n";
    out << "    \"k\": " << laws.k.evaluate(t_ref_k) << ",\n";
    out << "    \"sigma\": " << laws.sigma.evaluate(t_ref_k) << ",\n";
    out << "    \"rho\": " << laws.rho.evaluate(t_ref_k) << ",\n";
    out << "    \"beta\": " << laws.beta.evaluate(t_ref_k) << "\n";
    out << "  },\n";
    out << "  \"sigma_feedback_to_frozen_Q\": false,\n";
    out << "  \"sigma_feedback_note\": \"frozen_q mainline: sigma table may exist but does not update frozen heat source this round\"\n";
    out << "}\n";
    return true;
}

/** Curve CSV for plotting legacy vs thesis (T, mu/cp/k). */
inline bool writeDoc245MaterialCurveCompareCsv(const std::string &path,
                                               const OphelieFrenchGlassMaterialLaws &active,
                                               const OphelieFrenchGlassMaterialLaws &legacy)
{
    std::ofstream out(path);
    if (!out)
        return false;
    out << std::setprecision(10);
    out << "T_K,mu_active,mu_legacy,cp_active,cp_legacy,k_active,k_legacy\n";
    for (int T = 400; T <= 1800; T += 25)
    {
        const Real t = static_cast<Real>(T);
        out << t << "," << active.mu.evaluate(t) << "," << legacy.mu.evaluate(t) << ","
            << active.cp.evaluate(t) << "," << legacy.cp.evaluate(t) << "," << active.k.evaluate(t)
            << "," << legacy.k.evaluate(t) << "\n";
    }
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC245_MATERIAL_DATASET_H
