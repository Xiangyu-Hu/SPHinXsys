/**
 * @file gpu_paper_recording.h
 * @brief Shared run recording for GPU paper experiments (V4 templates).
 * @note GPU paper experiment (docs/gpu_paper_discussion). Not a production solver API.
 */
#ifndef SPHINXSYS_GPU_PAPER_RECORDING_H
#define SPHINXSYS_GPU_PAPER_RECORDING_H

#include <chrono>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sphinxsys_gpu_paper
{

inline std::string csvEscape(const std::string &s)
{
    bool need = s.find_first_of(",\"\n\r") != std::string::npos;
    if (!need)
        return s;
    std::string o = "\"";
    for (char c : s)
    {
        if (c == '"')
            o += "\"\"";
        else
            o += c;
    }
    o += "\"";
    return o;
}

/** Wall-clock stage accumulator (exclusive siblings under one parent). */
class StageClock
{
  public:
    void start(const std::string &stage)
    {
        active_ = stage;
        t0_ = std::chrono::steady_clock::now();
    }
    void stop()
    {
        if (active_.empty())
            return;
        auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
        wall_s_[active_] += dt;
        calls_[active_] += 1;
        active_.clear();
    }
    void add(const std::string &stage, double seconds, int calls = 1)
    {
        wall_s_[stage] += seconds;
        calls_[stage] += calls;
    }
    const std::map<std::string, double> &walls() const { return wall_s_; }
    const std::map<std::string, int> &calls() const { return calls_; }

  private:
    std::string active_;
    std::chrono::steady_clock::time_point t0_{};
    std::map<std::string, double> wall_s_;
    std::map<std::string, int> calls_;
};

struct MissingMetric
{
    std::string metric;
    std::string status; // unsupported|permission_denied|not_collected|not_applicable|...
    std::string reason;
    std::string evidence;
    std::string next_action;
};

struct ProbeSample
{
    double physical_time_s = 0.0;
    double requested_time_s = 0.0;
    std::string probe_id;
    std::string body;
    std::string quantity;
    double value = 0.0;
    std::string unit;
    double x = 0.0, y = 0.0, z = 0.0;
};

/**
 * Creates a run directory and writes V4-aligned CSV/JSON skeletons.
 * Metrics not yet measured are registered into missing_metrics at finalize.
 */
class RunRecorder
{
  public:
    RunRecorder(std::string run_dir, std::string run_id, std::string case_id, std::string mode);

    const std::string &runDir() const { return run_dir_; }
    const std::string &runId() const { return run_id_; }

    void setMeta(const std::string &key, const std::string &value);
    void registerMissing(MissingMetric m);
    void addProbe(const ProbeSample &p);
    void addDiagnostic(double t, const std::string &body, const std::string &quantity,
                       double value, const std::string &unit, const std::string &status = "measured");
    void addStepBlock(const std::string &step_type, int calls, double dt_min, double dt_mean, double dt_max,
                      int sort_calls, int rel_inner, int rel_contact, const std::string &dt_limiter,
                      const std::string &status = "measured");
    void addBodyStats(double t, const std::string &body, int n_active, int n_capacity,
                      const std::string &neighbor_kind, double directed_pairs, double mean_n, double max_n,
                      double occupied_cells, double cell_mean, double cell_max, const std::string &mode);
    void addThroughput(const std::string &module, double update_count, double neighbor_visits,
                       double denom_wall_s, const std::string &step_type, const std::string &count_method,
                       const std::string &pair_def, const std::string &status);
    void addMemory(const std::string &scope, const std::string &metric, double value_bytes,
                   double n_denom, const std::string &method, const std::string &status,
                   const std::string &missing_reason = "");
    void addTransfer(const std::string &parent, const std::string &category, const std::string &direction,
                     int calls, double bytes, double host_wall, double device_event,
                     const std::string &status, const std::string &missing_reason, const std::string &evidence);
    void addSync(const std::string &parent, const std::string &operation, int calls, double host_wall,
                 bool overlaps, const std::string &status, const std::string &notes);
    void addTelemetry(double wall_ts, double power, double temp, double sm_clock, double mem_clock,
                      double util, double mem_used, double rss, const std::string &throttle,
                      const std::string &status);
    void addKernelProfilePlaceholder(const std::string &kernel, const std::string &scope,
                                     const std::string &status, const std::string &missing_reason);
    void addEnergyPlaceholder(const std::string &scope, const std::string &status, const std::string &reason);

    StageClock &stages() { return stages_; }

    /** Write headers / empty applicable files once at start. */
    void writeHeaders();

    /** Flush all buffered rows and summary; apply default missings for unused tool metrics. */
    void finalize(const std::string &status, double t_start, double actual_t_end,
                  double wall_process, double wall_init, double wall_warmup, double wall_prewindow,
                  double wall_solver, double wall_solver_inclusive, double wall_diagnostics,
                  double wall_vtp, double wall_restart, double wall_validation, double wall_finalize,
                  double peak_vram_sampled, double device_alloc_peak, double peak_rss);

    /** Ensure required tool-attempt missings exist (call before first formal case suite). */
    void seedDefaultToolMissings();

  private:
    std::string path(const std::string &name) const { return run_dir_ + "/" + name; }
    void writeManifestJson() const;

    std::string run_dir_;
    std::string run_id_;
    std::string case_id_;
    std::string mode_;
    std::map<std::string, std::string> meta_;
    std::vector<MissingMetric> missing_;
    std::vector<ProbeSample> probes_;
    StageClock stages_;

    std::ofstream probes_ofs_;
    std::ofstream diag_ofs_;
    std::ofstream step_ofs_;
    std::ofstream body_ofs_;
    std::ofstream throughput_ofs_;
    std::ofstream memory_ofs_;
    std::ofstream transfer_ofs_;
    std::ofstream sync_ofs_;
    std::ofstream telemetry_ofs_;
    std::ofstream kernel_ofs_;
    std::ofstream energy_ofs_;
    bool headers_written_ = false;
};

/** Parse simple schedule JSON produced by freeze_output_schedule.py (probe_times array). */
std::vector<double> loadProbeTimesFromBoundSchedule(const std::string &json_path);
std::vector<double> loadFieldTimesFromBoundSchedule(const std::string &json_path);

/** Advance schedule cursor: return true if physical_time reached next requested sample. */
bool dueForSample(double physical_time, const std::vector<double> &times, size_t &cursor,
                  double *requested_out = nullptr);

} // namespace sphinxsys_gpu_paper

#endif
