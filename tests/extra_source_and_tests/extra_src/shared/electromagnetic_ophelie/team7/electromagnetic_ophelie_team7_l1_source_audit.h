#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_L1_SOURCE_AUDIT_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_L1_SOURCE_AUDIT_H

#include "electromagnetic_ophelie_team7_coil_path_source.h"
#include "electromagnetic_ophelie_team7_native_geometry.h"
#include "electromagnetic_ophelie_team7_probe.h"
#include "electromagnetic_ophelie_team7_validation.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct Team7L1SourceAuditRecord
{
    std::string category;
    std::string test_name;
    Real metric_value = 0.0;
    Real threshold = 0.0;
    bool passed = false;
    std::string notes;
};

struct Team7L1SourceAuditSummary
{
    OphelieCoilSourceModel source_model = OphelieCoilSourceModel::VolumeRacetrack;
    Team7BzProbeLine probe_line = Team7BzProbeLine::A1B1;
    Real coil_source_scale = 1.0;
    Real frequency_hz = 50.0;
    std::string reference_csv_kind = "phase0_total_muFEM";
    size_t n_probes = 0;
    Real bz_source_rms_abs_mT = 0.0;
    Real bz_source_rms_rel = 0.0;
    Real peak_sim_mT = 0.0;
    Real peak_ref_mT = 0.0;
    Real peak_x_sim_mm = 0.0;
    Real peak_x_ref_mm = 0.0;
    Real correlation = 0.0;
    Real best_fit_scale_l2 = 1.0;
    Real best_fit_scale_peak = 1.0;
    /** RMS_rel between volume-racetrack coil Bz and filament-racetrack at probes (diagnostic). */
    Real volume_vs_filament_rms_rel = 0.0;
    /** RMS_rel between ref phase0 and (coil_sim + ind_imag_skin) when L2 decomp available. */
    Real ref_total_vs_coil_plus_skin_rms_rel = 0.0;
    bool coil_path_audit_ok = false;
    bool reference_loaded = false;
};

inline Real team7L1PearsonCorrelation(const StdVec<Real> &sim_mT, const StdVec<Real> &ref_mT)
{
    const size_t n = std::min(sim_mT.size(), ref_mT.size());
    if (n < 2)
    {
        return 0.0;
    }
    Real mean_sim = 0.0;
    Real mean_ref = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        mean_sim += sim_mT[i];
        mean_ref += ref_mT[i];
    }
    mean_sim /= static_cast<Real>(n);
    mean_ref /= static_cast<Real>(n);
    Real cov = 0.0;
    Real var_sim = 0.0;
    Real var_ref = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real ds = sim_mT[i] - mean_sim;
        const Real dr = ref_mT[i] - mean_ref;
        cov += ds * dr;
        var_sim += ds * ds;
        var_ref += dr * dr;
    }
    return cov / (std::sqrt(var_sim * var_ref) + TinyReal);
}

inline Real team7L1BestFitScaleL2(const StdVec<Real> &sim_mT, const StdVec<Real> &ref_mT)
{
    const size_t n = std::min(sim_mT.size(), ref_mT.size());
    Real dot_sr = 0.0;
    Real dot_ss = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        dot_sr += sim_mT[i] * ref_mT[i];
        dot_ss += sim_mT[i] * sim_mT[i];
    }
    return dot_ss > TinyReal ? dot_sr / dot_ss : Real(1);
}

inline void team7L1ExtractSimRefVectors(const StdVec<Team7BzProbePoint> &coil_probes, StdVec<Real> &sim_mT,
                                        StdVec<Real> &ref_mT)
{
    sim_mT.clear();
    ref_mT.clear();
    sim_mT.reserve(coil_probes.size());
    ref_mT.reserve(coil_probes.size());
    for (const Team7BzProbePoint &probe : coil_probes)
    {
        sim_mT.push_back(probe.bz_sim_mT);
        ref_mT.push_back(probe.bz_ref_phase0_mT);
    }
}

inline void team7L1AddGeometryAuditRecords(const OphelieTeam7NativeDerivedGeometry &derived,
                                           const OphelieTeam7NativeMesh &mesh, Team7BzProbeLine probe_line,
                                           Real mm_to_m, StdVec<Team7L1SourceAuditRecord> &records)
{
    const Real inv_mm = 1.0 / mm_to_m;
    const Vec3d coil_lo = vec3FromVecdMm(derived.coil_bbox_.lower_, inv_mm);
    const Vec3d coil_hi = vec3FromVecdMm(derived.coil_bbox_.upper_, inv_mm);
    const Vec3d plate_lo = vec3FromVecdMm(derived.plate_bbox_.lower_, inv_mm);
    const Vec3d plate_hi = vec3FromVecdMm(derived.plate_bbox_.upper_, inv_mm);
    const Team7BzProbeLineSpec spec = team7BzProbeLineSpec(probe_line);

    const auto add = [&](const std::string &name, Real metric, Real threshold, bool pass, const std::string &notes)
    {
        records.push_back({"geometry", name, metric, threshold, pass, notes});
    };

    const Real dz_above_plate_top_mm = spec.z_mm - plate_hi[2];
    const Real dz_below_coil_bottom_mm = coil_lo[2] - spec.z_mm;
    const bool probe_in_air_gap = spec.z_mm > plate_hi[2] && spec.z_mm < coil_lo[2];
    add("probe_z_above_plate_top_mm", dz_above_plate_top_mm, 50.0, dz_above_plate_top_mm > Real(0),
        "A1-B1 z should be above plate top (air gap)");
    add("probe_z_below_coil_bottom_mm", dz_below_coil_bottom_mm, 50.0, dz_below_coil_bottom_mm > Real(0),
        "A1-B1 z should be below coil bottom");
    add("probe_z_in_air_gap", probe_in_air_gap ? Real(1) : Real(0), Real(0.5), probe_in_air_gap,
        "probe between plate top and coil bottom");
    add("probe_y_in_coil_y_span", (spec.y_mm >= coil_lo[1] && spec.y_mm <= coil_hi[1]) ? Real(1) : Real(0),
        Real(0.5), spec.y_mm >= coil_lo[1] && spec.y_mm <= coil_hi[1], "probe y within coil y bbox");
    add("coil_x_span_mm", coil_hi[0] - coil_lo[0], 400.0, coil_hi[0] > coil_lo[0],
        "coil x extent for peak placement audit");
    add("plate_x_span_mm", plate_hi[0] - plate_lo[0], 400.0, plate_hi[0] > plate_lo[0], "plate x extent");
    add("stl_scale_to_meter", mesh.stl_scale_to_meter_, 1.0, mesh.stl_scale_to_meter_ > TinyReal, "mm→m scale");
    (void)plate_lo;
}

inline void team7L1AddReferenceMetadataRecords(bool reference_loaded, Team7BzProbeLine probe_line,
                                             const StdVec<Team7BzProbePoint> &reference_probes,
                                             StdVec<Team7L1SourceAuditRecord> &records)
{
    const auto add = [&](const std::string &name, Real metric, Real threshold, bool pass, const std::string &notes)
    {
        records.push_back({"reference", name, metric, threshold, pass, notes});
    };
    add("reference_csv_loaded", reference_loaded ? Real(1) : Real(0), Real(0.5), reference_loaded,
        team7BzProbeLineSpec(probe_line).csv_basename);
    add("reference_is_phase0_total_not_pure_source", Real(1), Real(0.5), false,
        "TEAM7_Bz_*_reference_mT.csv is phase0 total probe; NOT f=0 source-only");
    if (!reference_probes.empty())
    {
        const Team7BzProbeLineSpec spec = team7BzProbeLineSpec(probe_line);
        const Real y_err_mm = std::abs(reference_probes.front().position_m[1] / 1.0e-3 - spec.y_mm);
        const Real z_err_mm = std::abs(reference_probes.front().position_m[2] / 1.0e-3 - spec.z_mm);
        add("probe_line_y_mm_match", y_err_mm, 0.01, y_err_mm < Real(0.01), "loaded probe y vs CSV spec");
        add("probe_line_z_mm_match", z_err_mm, 0.01, z_err_mm < Real(0.01), "loaded probe z vs CSV spec");
        add("n_reference_probes", static_cast<Real>(reference_probes.size()), 1.0,
            reference_probes.size() >= 10, "probe count along x");
    }
}

inline Team7L1SourceAuditSummary computeTeam7L1SourceAuditSummary(
    const StdVec<Team7BzProbePoint> &coil_probes, OphelieCoilSourceModel source_model, Team7BzProbeLine probe_line,
    Real coil_source_scale, Real frequency_hz, bool coil_path_audit_ok, bool reference_loaded,
    const Team7BzCompareMetrics &coil_metrics)
{
    Team7L1SourceAuditSummary summary;
    summary.source_model = source_model;
    summary.probe_line = probe_line;
    summary.coil_source_scale = coil_source_scale;
    summary.frequency_hz = frequency_hz;
    summary.coil_path_audit_ok = coil_path_audit_ok;
    summary.reference_loaded = reference_loaded;
    summary.n_probes = coil_probes.size();
    summary.bz_source_rms_abs_mT = coil_metrics.rms_abs_error_mT;
    summary.bz_source_rms_rel = coil_metrics.rms_rel_error;
    summary.peak_sim_mT = coil_metrics.peak_sim_mT;
    summary.peak_ref_mT = coil_metrics.peak_ref_mT;
    if (coil_metrics.n_probes > 0 && !coil_probes.empty())
    {
        summary.peak_x_sim_mm = coil_probes[coil_metrics.peak_sim_index].x_mm;
        summary.peak_x_ref_mm = coil_probes[coil_metrics.peak_ref_index].x_mm;
    }
    StdVec<Real> sim_mT;
    StdVec<Real> ref_mT;
    team7L1ExtractSimRefVectors(coil_probes, sim_mT, ref_mT);
    summary.correlation = team7L1PearsonCorrelation(sim_mT, ref_mT);
    summary.best_fit_scale_l2 = team7L1BestFitScaleL2(sim_mT, ref_mT);
    if (std::abs(summary.peak_sim_mT) > TinyReal)
    {
        summary.best_fit_scale_peak = summary.peak_ref_mT / summary.peak_sim_mT;
    }
    return summary;
}

inline void team7L1FilamentCrossCheckAtProbes(const BoundingBoxd &coil_bbox, const OphelieTeam7CoilPathSourceSpec &spec,
                                              Real ds_m, Real mu0, Real softening_length,
                                              const StdVec<Team7BzProbePoint> &reference_probes,
                                              const StdVec<Team7BzProbePoint> &volume_coil_probes,
                                              Real &volume_vs_filament_rms_rel)
{
    volume_vs_filament_rms_rel = 0.0;
    if (reference_probes.empty() || volume_coil_probes.empty())
    {
        return;
    }
    StdVec<OphelieCurrentMomentSample> filament_moments;
    buildOphelieTeam7FilamentMomentsFromCoilPath(coil_bbox, spec, ds_m, filament_moments);
    if (filament_moments.empty())
    {
        return;
    }
    StdVec<Team7BzProbePoint> filament_probes;
    evaluateTeam7FilamentBiotSavartBzAtProbes(filament_moments, mu0, softening_length, reference_probes,
                                              filament_probes);
    Real sum_err_sq = 0.0;
    Real sum_vol_sq = 0.0;
    const size_t n = std::min(volume_coil_probes.size(), filament_probes.size());
    for (size_t i = 0; i < n; ++i)
    {
        const Real vol = volume_coil_probes[i].bz_sim_mT;
        const Real fil = filament_probes[i].bz_sim_mT;
        const Real err = vol - fil;
        sum_err_sq += err * err;
        sum_vol_sq += vol * vol;
    }
    volume_vs_filament_rms_rel = std::sqrt(sum_err_sq) / (std::sqrt(sum_vol_sq) + TinyReal);
}

inline Real team7L1RefTotalVsCoilPlusSkinRmsRel(const StdVec<Team7BzProbePoint> &reference_probes,
                                                const StdVec<Team7BzProbePoint> &coil_probes,
                                                const StdVec<Team7ProbeBzDecomposition> &decomp)
{
    if (reference_probes.empty() || coil_probes.size() != decomp.size())
    {
        return 0.0;
    }
    Real sum_err_sq = 0.0;
    Real sum_ref_sq = 0.0;
    for (size_t i = 0; i < reference_probes.size(); ++i)
    {
        const Real ref = reference_probes[i].bz_ref_phase0_mT;
        const Real hybrid = coil_probes[i].bz_sim_mT + decomp[i].bz_ind_imag_skin_mT;
        const Real err = hybrid - ref;
        sum_err_sq += err * err;
        sum_ref_sq += ref * ref;
    }
    return std::sqrt(sum_err_sq) / (std::sqrt(sum_ref_sq) + TinyReal);
}

inline StdVec<Team7L1SourceAuditRecord> computeTeam7L1SourceAudit(
    const StdVec<Team7BzProbePoint> &reference_probes, const StdVec<Team7BzProbePoint> &coil_probes,
    const StdVec<Team7ProbeBzDecomposition> *decomp, bool reference_ok, OphelieCoilSourceModel source_model,
    Real coil_source_scale, Team7BzProbeLine probe_line, const OphelieTeam7NativeDerivedGeometry &derived,
    const OphelieTeam7NativeMesh &mesh, const OphelieTeam7CoilPathSourceSpec &coil_path_spec,
    const OphelieParameters &params, bool coil_path_audit_ok, const Team7BzCompareMetrics &coil_metrics,
    Real filament_ds_m, Team7L1SourceAuditSummary &summary_out)
{
    StdVec<Team7L1SourceAuditRecord> records;
    summary_out = computeTeam7L1SourceAuditSummary(coil_probes, source_model, probe_line, coil_source_scale,
                                                   params.frequency_, coil_path_audit_ok, reference_ok, coil_metrics);

    team7L1AddReferenceMetadataRecords(reference_ok, probe_line, reference_probes, records);
    team7L1AddGeometryAuditRecords(derived, mesh, probe_line, 1.0e-3, records);

    const auto add_source = [&](const std::string &name, Real metric, Real threshold, bool pass,
                                const std::string &notes)
    {
        records.push_back({"source_bz", name, metric, threshold, pass, notes});
    };

    add_source("bz_source_rms_abs_mT", summary_out.bz_source_rms_abs_mT, 2.0,
               summary_out.bz_source_rms_abs_mT < Real(2.0), "coil-only Bz vs phase0-total ref (diagnostic)");
    add_source("bz_source_rms_rel", summary_out.bz_source_rms_rel, 0.70, summary_out.bz_source_rms_rel < Real(0.70),
               "relative RMS; L1 smoke threshold 0.70");
    add_source("peak_sim_over_ref", summary_out.peak_sim_mT / (summary_out.peak_ref_mT + TinyReal), 1.5,
               std::abs(summary_out.peak_sim_mT / (summary_out.peak_ref_mT + TinyReal) - Real(1)) < Real(0.35),
               "peak amplitude ratio");
    add_source("peak_x_offset_mm", summary_out.peak_x_sim_mm - summary_out.peak_x_ref_mm, 100.0,
               std::abs(summary_out.peak_x_sim_mm - summary_out.peak_x_ref_mm) < Real(100.0),
               "peak x sim minus ref (TEAM7 often ~+70–90 mm)");
    add_source("correlation", summary_out.correlation, 0.90, summary_out.correlation > Real(0.90),
               "Pearson(coil_sim, ref_phase0)");
    add_source("best_fit_scale_l2", summary_out.best_fit_scale_l2, 1.0,
               summary_out.best_fit_scale_l2 > Real(0.5) && summary_out.best_fit_scale_l2 < Real(1.5),
               "argmin ||scale*sim-ref||_2");
    add_source("best_fit_scale_peak", summary_out.best_fit_scale_peak, 1.0,
               summary_out.best_fit_scale_peak > Real(0.5) && summary_out.best_fit_scale_peak < Real(1.5),
               "peak_ref/peak_sim");
    add_source("coil_path_audit", coil_path_audit_ok ? Real(1) : Real(0), Real(0.5), coil_path_audit_ok,
               "NI / integrated J_src audit");

    if (source_model == OphelieCoilSourceModel::VolumeRacetrack && reference_ok)
    {
        team7L1FilamentCrossCheckAtProbes(derived.coil_bbox_, coil_path_spec, filament_ds_m, params.mu0_,
                                          params.softening_length_, reference_probes, coil_probes,
                                          summary_out.volume_vs_filament_rms_rel);
        records.push_back({"cross_check", "volume_vs_filament_rms_rel", summary_out.volume_vs_filament_rms_rel, 0.30,
                           summary_out.volume_vs_filament_rms_rel < Real(0.30),
                           "volume-racetrack vs filament-racetrack at probes"});
    }

    if (decomp != nullptr && !decomp->empty() && reference_ok)
    {
        summary_out.ref_total_vs_coil_plus_skin_rms_rel =
            team7L1RefTotalVsCoilPlusSkinRmsRel(reference_probes, coil_probes, *decomp);
        records.push_back(
            {"cross_check", "ref_phase0_vs_coil_plus_ind_skin_rms_rel",
             summary_out.ref_total_vs_coil_plus_skin_rms_rel, summary_out.bz_source_rms_rel,
             summary_out.ref_total_vs_coil_plus_skin_rms_rel < summary_out.bz_source_rms_rel,
             "if ref≈coil+skin then phase0 ref may be coil-dominated at probes"});
    }

    return records;
}

inline void printTeam7L1SourceAuditReport(const Team7L1SourceAuditSummary &summary,
                                          const StdVec<Team7L1SourceAuditRecord> &records)
{
    size_t n_pass = 0;
    size_t n_fail = 0;
    std::cout << "[team7] P6a L1 source/reference/probe audit:" << std::endl;
    std::cout << "[team7]   model=" << ophelieCoilSourceModelName(summary.source_model)
              << " line=" << team7BzProbeLineName(summary.probe_line) << " scale=" << summary.coil_source_scale
              << " ref_kind=" << summary.reference_csv_kind << std::endl;
    std::cout << "[team7]   Bz_source RMS_abs=" << summary.bz_source_rms_abs_mT
              << " RMS_rel=" << summary.bz_source_rms_rel << " corr=" << summary.correlation
              << " best_fit_L2=" << summary.best_fit_scale_l2 << " best_fit_peak=" << summary.best_fit_scale_peak
              << std::endl;
    std::cout << "[team7]   peak sim/ref=" << summary.peak_sim_mT << "/" << summary.peak_ref_mT
              << " peak_x sim/ref=" << summary.peak_x_sim_mm << "/" << summary.peak_x_ref_mm << std::endl;
    if (summary.volume_vs_filament_rms_rel > TinyReal)
    {
        std::cout << "[team7]   volume_vs_filament_rms_rel=" << summary.volume_vs_filament_rms_rel << std::endl;
    }
    if (summary.ref_total_vs_coil_plus_skin_rms_rel > TinyReal)
    {
        std::cout << "[team7]   ref_vs_coil+ind_skin_rms_rel=" << summary.ref_total_vs_coil_plus_skin_rms_rel
                  << std::endl;
    }
    for (const Team7L1SourceAuditRecord &record : records)
    {
        if (record.passed)
        {
            ++n_pass;
        }
        else
        {
            ++n_fail;
        }
        std::cout << "  " << (record.passed ? "PASS" : "FAIL") << " [" << record.category << "] " << record.test_name
                  << " metric=" << record.metric_value << " threshold=" << record.threshold;
        if (!record.notes.empty())
        {
            std::cout << " (" << record.notes << ")";
        }
        std::cout << std::endl;
    }
    std::cout << "[team7] P6a L1 audit summary: pass=" << n_pass << " fail=" << n_fail
              << " (diagnostic-only; phase0 ref is NOT pure source-only)" << std::endl;
}

inline void writeTeam7L1SourceAuditCsv(const std::string &detail_path, const std::string &summary_path,
                                       const Team7L1SourceAuditSummary &summary,
                                       const StdVec<Team7L1SourceAuditRecord> &records)
{
    namespace fs = std::filesystem;
    for (const std::string &path : {detail_path, summary_path})
    {
        const fs::path parent = fs::path(path).parent_path();
        if (!parent.empty())
        {
            fs::create_directories(parent);
        }
    }
    {
        std::ofstream out(detail_path);
        out << "category,test_name,metric,threshold,passed,notes\n";
        for (const Team7L1SourceAuditRecord &record : records)
        {
            out << record.category << "," << record.test_name << "," << record.metric_value << "," << record.threshold
                << "," << (record.passed ? 1 : 0) << ",\"" << record.notes << "\"\n";
        }
    }
    {
        std::ofstream out(summary_path);
        out << "source_model,probe_line,coil_source_scale,frequency_hz,reference_csv_kind,n_probes,"
               "bz_source_rms_abs_mT,bz_source_rms_rel,peak_sim_mT,peak_ref_mT,peak_x_sim_mm,peak_x_ref_mm,"
               "correlation,best_fit_scale_l2,best_fit_scale_peak,volume_vs_filament_rms_rel,"
               "ref_total_vs_coil_plus_skin_rms_rel,coil_path_audit_ok,reference_loaded\n";
        out << ophelieCoilSourceModelName(summary.source_model) << "," << team7BzProbeLineName(summary.probe_line)
            << "," << summary.coil_source_scale << "," << summary.frequency_hz << "," << summary.reference_csv_kind
            << "," << summary.n_probes << "," << summary.bz_source_rms_abs_mT << "," << summary.bz_source_rms_rel << ","
            << summary.peak_sim_mT << "," << summary.peak_ref_mT << "," << summary.peak_x_sim_mm << ","
            << summary.peak_x_ref_mm << "," << summary.correlation << "," << summary.best_fit_scale_l2 << ","
            << summary.best_fit_scale_peak << "," << summary.volume_vs_filament_rms_rel << ","
            << summary.ref_total_vs_coil_plus_skin_rms_rel << "," << (summary.coil_path_audit_ok ? 1 : 0) << ","
            << (summary.reference_loaded ? 1 : 0) << "\n";
    }
    std::cout << "[team7] P6a L1 source audit CSV: " << detail_path << " ; summary: " << summary_path << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_L1_SOURCE_AUDIT_H
