#include "gpu_paper_recording.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

namespace fs = std::filesystem;
namespace sphinxsys_gpu_paper
{

namespace
{
void ensureDir(const std::string &dir)
{
    fs::create_directories(dir);
}

std::string getenvOr(const char *k, const std::string &fallback)
{
    const char *v = std::getenv(k);
    return (v && *v) ? std::string(v) : fallback;
}
} // namespace

RunRecorder::RunRecorder(std::string run_dir, std::string run_id, std::string case_id, std::string mode)
    : run_dir_(std::move(run_dir)), run_id_(std::move(run_id)), case_id_(std::move(case_id)), mode_(std::move(mode))
{
    ensureDir(run_dir_);
    ensureDir(run_dir_ + "/snapshots");
    meta_["case_id"] = case_id_;
    meta_["mode"] = mode_;
    meta_["run_id"] = run_id_;
    meta_["record_kind"] = "measured";
    meta_["cuda_visible_devices"] = getenvOr("CUDA_VISIBLE_DEVICES", "");
}

void RunRecorder::setMeta(const std::string &key, const std::string &value) { meta_[key] = value; }

void RunRecorder::registerMissing(MissingMetric m) { missing_.push_back(std::move(m)); }

void RunRecorder::seedDefaultToolMissings()
{
    auto add_if_absent = [&](const std::string &metric, const std::string &status,
                             const std::string &reason, const std::string &next)
    {
        for (const auto &m : missing_)
            if (m.metric == metric)
                return;
        missing_.push_back({metric, status, reason, "predeclared_before_suite", next});
    };

    add_if_absent("kernel_profile_nsight", "not_collected",
                  "Nsight/kernel counters collected only on B2/B3 representative short profile window",
                  "run supplemental profile mode; fill kernel_profile.csv or keep missing");
    add_if_absent("roofline", "not_applicable",
                  "Requires trusted op/byte counts; not fabricated",
                  "optional representative profile only");
    add_if_absent("gpu_energy_integral_J", "not_collected",
                  "Derived only if telemetry power coverage is valid over matched window",
                  "enable telemetry sampling; derive in postprocess when coverage_ratio OK");
    add_if_absent("cpu_package_energy_J", "not_collected",
                  "Needs RAPL/perf permission; not assumed available",
                  "attempt on CPU baseline host; else permission_denied");
    add_if_absent("device_event_timing_full_queue", "not_collected",
                  "Event timing only when queue profiling covers all related work",
                  "enable on profile short window; never replace wall_solver");
    add_if_absent("multi_gpu_scaling", "not_applicable",
                  "Single-GPU paper scope", "none");
    add_if_absent("cross_vendor_gpu", "not_applicable",
                  "Out of scope for this campaign", "none");
    add_if_absent("cross_sph_software_benchmark", "not_applicable",
                  "Out of scope for this campaign", "none");
}

void RunRecorder::writeHeaders()
{
    if (headers_written_)
        return;
    probes_ofs_.open(path("probes.csv"));
    probes_ofs_ << "run_id,physical_time_s,probe_id,body,quantity,value,unit,x,y,z,requested_time_s\n";
    diag_ofs_.open(path("diagnostics.csv"));
    diag_ofs_ << "run_id,physical_time_s,body,quantity,value,unit,measurement_status\n";
    step_ofs_.open(path("step_stats.csv"));
    step_ofs_ << "run_id,step_type,calls,dt_min_s,dt_mean_s,dt_max_s,sort_calls,relation_inner_calls,"
                 "relation_contact_calls,dt_limiter,status\n";
    body_ofs_.open(path("body_stats.csv"));
    body_ofs_ << "run_id,physical_time_s,body,N_active,N_capacity,neighbor_kind,directed_pairs,"
                 "mean_neighbors,max_neighbors,occupied_cells,cell_occupancy_mean,cell_occupancy_max,"
                 "collection_mode\n";
    throughput_ofs_.open(path("throughput.csv"));
    throughput_ofs_ << "run_id,module,update_count,neighbor_visit_count,denominator_wall_s,MPPS,MNIPS,"
                       "step_type,count_method,pair_definition,measurement_status\n";
    memory_ofs_.open(path("memory_stats.csv"));
    memory_ofs_ << "run_id,scope,metric,value_bytes,N_denominator,measurement_method,measurement_status,"
                   "missing_reason\n";
    transfer_ofs_.open(path("transfer_timing.csv"));
    transfer_ofs_ << "run_id,parent_scope,category,direction,calls,bytes,host_wall_s,device_event_s,"
                     "measurement_status,missing_reason,evidence\n";
    sync_ofs_.open(path("sync_timing.csv"));
    sync_ofs_ << "run_id,parent_scope,operation,calls,host_wall_s,overlaps_device_work,measurement_status,notes\n";
    telemetry_ofs_.open(path("telemetry.csv"));
    telemetry_ofs_ << "run_id,wall_timestamp_s,gpu_power_W,gpu_temperature_C,sm_clock_MHz,memory_clock_MHz,"
                      "gpu_utilization_pct,device_memory_used_bytes,process_rss_bytes,throttle_reason,"
                      "measurement_status\n";
    kernel_ofs_.open(path("kernel_profile.csv"));
    kernel_ofs_ << "run_id,kernel,scope,tool_version,native_metric_name,value,unit,semantic_metric,"
                   "replay_passes,measurement_status,missing_reason,raw_report\n";
    energy_ofs_.open(path("energy_summary.csv"));
    energy_ofs_ << "run_id,scope,energy_J,time_coverage_ratio,integration_method,sample_count,"
                   "window_wall_s,measurement_status,missing_reason\n";

    // Placeholders with open schema for later fill
    {
        std::ofstream ofs(path("timing_scope.csv"));
        ofs << "run_id,scope_id,parent_scope_id,label,clock,aggregation,exclusive_siblings,wall_s,"
               "completion_boundary,notes\n";
    }
    {
        std::ofstream ofs(path("accuracy_cost.csv"));
        ofs << "run_id,reference_hash,quantity,error_norm,error_value,normalization,wall_solver_s,dp,N,"
               "measurement_status\n";
    }
    {
        std::ofstream ofs(path("instrumentation_check.csv"));
        ofs << "pair_id,case_id,baseline_run_id,instrumented_run_id,run_key_match,baseline_s,"
               "instrumented_s,overhead_ratio,instrumentation,record_kind\n";
    }
    {
        std::ofstream ofs(path("run_aggregate.csv"));
        ofs << "comparison_key,case_id,platform,precision,group,n,median_s,min_s,max_s,sample_sd_s,cv,"
               "cpu_reference_key,solver_speedup,record_kind\n";
    }

    seedDefaultToolMissings();
    headers_written_ = true;
}

void RunRecorder::addProbe(const ProbeSample &p)
{
    if (!headers_written_)
        writeHeaders();
    probes_ofs_ << csvEscape(run_id_) << "," << p.physical_time_s << "," << csvEscape(p.probe_id) << ","
                << csvEscape(p.body) << "," << csvEscape(p.quantity) << "," << p.value << ","
                << csvEscape(p.unit) << "," << p.x << "," << p.y << "," << p.z << "," << p.requested_time_s
                << "\n";
}

void RunRecorder::addDiagnostic(double t, const std::string &body, const std::string &quantity, double value,
                                const std::string &unit, const std::string &status)
{
    if (!headers_written_)
        writeHeaders();
    diag_ofs_ << csvEscape(run_id_) << "," << t << "," << csvEscape(body) << "," << csvEscape(quantity) << ","
              << value << "," << csvEscape(unit) << "," << csvEscape(status) << "\n";
}

void RunRecorder::addStepBlock(const std::string &step_type, int calls, double dt_min, double dt_mean,
                               double dt_max, int sort_calls, int rel_inner, int rel_contact,
                               const std::string &dt_limiter, const std::string &status)
{
    if (!headers_written_)
        writeHeaders();
    step_ofs_ << csvEscape(run_id_) << "," << csvEscape(step_type) << "," << calls << "," << dt_min << ","
              << dt_mean << "," << dt_max << "," << sort_calls << "," << rel_inner << "," << rel_contact << ","
              << csvEscape(dt_limiter) << "," << csvEscape(status) << "\n";
}

void RunRecorder::addBodyStats(double t, const std::string &body, int n_active, int n_capacity,
                               const std::string &neighbor_kind, double directed_pairs, double mean_n,
                               double max_n, double occupied_cells, double cell_mean, double cell_max,
                               const std::string &mode)
{
    if (!headers_written_)
        writeHeaders();
    body_ofs_ << csvEscape(run_id_) << "," << t << "," << csvEscape(body) << "," << n_active << ","
              << n_capacity << "," << csvEscape(neighbor_kind) << "," << directed_pairs << "," << mean_n << ","
              << max_n << "," << occupied_cells << "," << cell_mean << "," << cell_max << ","
              << csvEscape(mode) << "\n";
}

void RunRecorder::addThroughput(const std::string &module, double update_count, double neighbor_visits,
                                double denom_wall_s, const std::string &step_type,
                                const std::string &count_method, const std::string &pair_def,
                                const std::string &status)
{
    if (!headers_written_)
        writeHeaders();
    double mpps = (denom_wall_s > 0.0) ? (update_count / 1.0e6 / denom_wall_s) : 0.0;
    double mnips = (denom_wall_s > 0.0 && neighbor_visits > 0.0) ? (neighbor_visits / 1.0e6 / denom_wall_s) : 0.0;
    throughput_ofs_ << csvEscape(run_id_) << "," << csvEscape(module) << "," << update_count << ","
                    << neighbor_visits << "," << denom_wall_s << "," << mpps << "," << mnips << ","
                    << csvEscape(step_type) << "," << csvEscape(count_method) << "," << csvEscape(pair_def)
                    << "," << csvEscape(status) << "\n";
}

void RunRecorder::addMemory(const std::string &scope, const std::string &metric, double value_bytes,
                            double n_denom, const std::string &method, const std::string &status,
                            const std::string &missing_reason)
{
    if (!headers_written_)
        writeHeaders();
    memory_ofs_ << csvEscape(run_id_) << "," << csvEscape(scope) << "," << csvEscape(metric) << ","
                << value_bytes << "," << n_denom << "," << csvEscape(method) << "," << csvEscape(status) << ","
                << csvEscape(missing_reason) << "\n";
}

void RunRecorder::addTransfer(const std::string &parent, const std::string &category,
                              const std::string &direction, int calls, double bytes, double host_wall,
                              double device_event, const std::string &status,
                              const std::string &missing_reason, const std::string &evidence)
{
    if (!headers_written_)
        writeHeaders();
    transfer_ofs_ << csvEscape(run_id_) << "," << csvEscape(parent) << "," << csvEscape(category) << ","
                  << csvEscape(direction) << "," << calls << "," << bytes << "," << host_wall << ","
                  << device_event << "," << csvEscape(status) << "," << csvEscape(missing_reason) << ","
                  << csvEscape(evidence) << "\n";
}

void RunRecorder::addSync(const std::string &parent, const std::string &operation, int calls,
                          double host_wall, bool overlaps, const std::string &status,
                          const std::string &notes)
{
    if (!headers_written_)
        writeHeaders();
    sync_ofs_ << csvEscape(run_id_) << "," << csvEscape(parent) << "," << csvEscape(operation) << ","
              << calls << "," << host_wall << "," << (overlaps ? "true" : "false") << ","
              << csvEscape(status) << "," << csvEscape(notes) << "\n";
}

void RunRecorder::addTelemetry(double wall_ts, double power, double temp, double sm_clock, double mem_clock,
                               double util, double mem_used, double rss, const std::string &throttle,
                               const std::string &status)
{
    if (!headers_written_)
        writeHeaders();
    telemetry_ofs_ << csvEscape(run_id_) << "," << wall_ts << "," << power << "," << temp << "," << sm_clock
                   << "," << mem_clock << "," << util << "," << mem_used << "," << rss << ","
                   << csvEscape(throttle) << "," << csvEscape(status) << "\n";
}

void RunRecorder::addKernelProfilePlaceholder(const std::string &kernel, const std::string &scope,
                                              const std::string &status, const std::string &missing_reason)
{
    if (!headers_written_)
        writeHeaders();
    kernel_ofs_ << csvEscape(run_id_) << "," << csvEscape(kernel) << "," << csvEscape(scope) << ",,,,,"
                << ",," << csvEscape(status) << "," << csvEscape(missing_reason) << ",\n";
}

void RunRecorder::addEnergyPlaceholder(const std::string &scope, const std::string &status,
                                       const std::string &reason)
{
    if (!headers_written_)
        writeHeaders();
    energy_ofs_ << csvEscape(run_id_) << "," << csvEscape(scope) << ",,,,,"
                << "," << csvEscape(status) << "," << csvEscape(reason) << "\n";
}

void RunRecorder::writeManifestJson() const
{
    std::ofstream ofs(path("run_manifest.json"));
    ofs << "{\n";
    ofs << "  \"run_id\": \"" << run_id_ << "\",\n";
    ofs << "  \"case_id\": \"" << case_id_ << "\",\n";
    ofs << "  \"mode\": \"" << mode_ << "\",\n";
    ofs << "  \"record_kind\": \"measured\",\n";
    for (const auto &kv : meta_)
    {
        if (kv.first == "run_id" || kv.first == "case_id" || kv.first == "mode")
            continue;
        ofs << "  \"" << kv.first << "\": \"" << kv.second << "\",\n";
    }
    ofs << "  \"_notice\": \"Fill remaining V4 fields via runner/preflight; nulls not written as 0.\"\n";
    ofs << "}\n";
}

void RunRecorder::finalize(const std::string &status, double t_start, double actual_t_end,
                           double wall_process, double wall_init, double wall_warmup, double wall_prewindow,
                           double wall_solver, double wall_solver_inclusive, double wall_diagnostics,
                           double wall_vtp, double wall_restart, double wall_validation, double wall_finalize,
                           double peak_vram_sampled, double device_alloc_peak, double peak_rss)
{
    if (!headers_written_)
        writeHeaders();

    // stage_timing
    {
        std::ofstream ofs(path("stage_timing.csv"));
        ofs << "run_id,stage,parent_scope,calls,wall_s,device_event_s,inclusive_or_exclusive,"
               "timing_method,concurrent,notes\n";
        for (const auto &kv : stages_.walls())
        {
            int c = 0;
            auto it = stages_.calls().find(kv.first);
            if (it != stages_.calls().end())
                c = it->second;
            ofs << csvEscape(run_id_) << "," << csvEscape(kv.first) << ",solver," << c << "," << kv.second
                << ",,exclusive,host_steady_clock,false,\n";
        }
    }

    double accounted = wall_init + wall_warmup + wall_prewindow + wall_solver + wall_diagnostics + wall_vtp +
                       wall_restart + wall_validation + wall_finalize;
    double unaccounted = wall_process - accounted;

    {
        std::ofstream ofs(path("run_summary.csv"));
        ofs << "run_id,record_kind,status,t_start,actual_t_end,wall_process_s,wall_init_s,wall_warmup_s,"
               "wall_prewindow_s,wall_solver_s,wall_solver_window_inclusive_s,wall_diagnostics_s,"
               "wall_vtp_output_s,wall_restart_output_s,wall_validation_postprocess_s,wall_finalize_s,"
               "wall_unaccounted_s,peak_vram_sampled_bytes,device_allocated_peak_bytes,peak_rss_bytes\n";
        ofs << csvEscape(run_id_) << ",measured," << csvEscape(status) << "," << t_start << ","
            << actual_t_end << "," << wall_process << "," << wall_init << "," << wall_warmup << ","
            << wall_prewindow << "," << wall_solver << "," << wall_solver_inclusive << ","
            << wall_diagnostics << "," << wall_vtp << "," << wall_restart << "," << wall_validation << ","
            << wall_finalize << "," << unaccounted << "," << peak_vram_sampled << "," << device_alloc_peak
            << "," << peak_rss << "\n";
    }

    {
        std::ofstream ofs(path("missing_metrics.csv"));
        ofs << "run_id,metric,measurement_status,reason,evidence,next_action\n";
        for (const auto &m : missing_)
        {
            ofs << csvEscape(run_id_) << "," << csvEscape(m.metric) << "," << csvEscape(m.status) << ","
                << csvEscape(m.reason) << "," << csvEscape(m.evidence) << "," << csvEscape(m.next_action)
                << "\n";
        }
    }

    writeManifestJson();

    {
        std::ofstream ofs(path("exit_status.json"));
        ofs << "{\"run_id\":\"" << run_id_ << "\",\"status\":\"" << status << "\"}\n";
    }

    probes_ofs_.flush();
    diag_ofs_.flush();
    step_ofs_.flush();
    body_ofs_.flush();
    throughput_ofs_.flush();
    memory_ofs_.flush();
    transfer_ofs_.flush();
    sync_ofs_.flush();
    telemetry_ofs_.flush();
    kernel_ofs_.flush();
    energy_ofs_.flush();

    std::cout << "[gpu_paper] finalized run_dir=" << run_dir_ << " status=" << status << std::endl;
}

static std::vector<double> parseJsonNumberArray(const std::string &text, const std::string &key)
{
    std::vector<double> out;
    const std::string pat = "\"" + key + "\"";
    auto pos = text.find(pat);
    if (pos == std::string::npos)
        return out;
    pos = text.find('[', pos);
    if (pos == std::string::npos)
        return out;
    auto end = text.find(']', pos);
    if (end == std::string::npos)
        return out;
    std::string body = text.substr(pos + 1, end - pos - 1);
    std::stringstream ss(body);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        try
        {
            out.push_back(std::stod(item));
        }
        catch (...)
        {
        }
    }
    return out;
}

std::vector<double> loadProbeTimesFromBoundSchedule(const std::string &json_path)
{
    std::ifstream ifs(json_path);
    if (!ifs)
        return {};
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    auto v = parseJsonNumberArray(buffer.str(), "probe_targets");
    if (v.empty())
        v = parseJsonNumberArray(buffer.str(), "probe_times");
    if (v.empty())
        v = parseJsonNumberArray(buffer.str(), "probes");
    return v;
}

std::vector<double> loadFieldTimesFromBoundSchedule(const std::string &json_path)
{
    std::ifstream ifs(json_path);
    if (!ifs)
        return {};
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    auto v = parseJsonNumberArray(buffer.str(), "field_targets");
    if (v.empty())
        v = parseJsonNumberArray(buffer.str(), "field_times");
    if (v.empty())
        v = parseJsonNumberArray(buffer.str(), "fields");
    return v;
}

bool dueForSample(double physical_time, const std::vector<double> &times, size_t &cursor,
                  double *requested_out)
{
    if (cursor >= times.size())
        return false;
    if (physical_time + 1.0e-15 < times[cursor])
        return false;
    if (requested_out)
        *requested_out = times[cursor];
    ++cursor;
    // If one step crosses multiple targets, only keep the first actual state (V4 rule).
    while (cursor < times.size() && physical_time + 1.0e-15 >= times[cursor])
        ++cursor;
    return true;
}

} // namespace sphinxsys_gpu_paper
