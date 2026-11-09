#ifndef PAPER_BENCHMARK_RECORDER_H
#define PAPER_BENCHMARK_RECORDER_H

#include "benchmark_config.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef SPHINXSYS_BENCHMARK_GIT_COMMIT
#define SPHINXSYS_BENCHMARK_GIT_COMMIT "unknown"
#endif
#ifndef SPHINXSYS_BENCHMARK_BUILD
#define SPHINXSYS_BENCHMARK_BUILD "unknown"
#endif
#ifndef SPHINXSYS_BENCHMARK_COMPILER
#define SPHINXSYS_BENCHMARK_COMPILER "unknown"
#endif

namespace paper_bench
{
inline std::string formatMetric(double value)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10)
           << value;
    return stream.str();
}

inline std::string peakRssKb()
{
#if defined(__linux__)
    std::ifstream status("/proc/self/status");
    if (!status)
    {
        return "unavailable";
    }
    std::string line;
    while (std::getline(status, line))
    {
        if (line.rfind("VmHWM:", 0) != 0)
        {
            continue;
        }
        std::istringstream parser(line.substr(6));
        long kilobytes = -1;
        std::string unit;
        parser >> kilobytes >> unit;
        if (!parser || kilobytes < 0)
        {
            return "unavailable";
        }
        return std::to_string(kilobytes);
    }
#endif
    return "unavailable";
}

struct BenchmarkSummary
{
    double dp = 0.0;
    std::size_t initial_particle_count = 0;
    std::size_t particle_count = 0;
    double physical_time = 0.0;
    std::size_t outer_steps = 0;
    std::size_t advection_steps = 0;
    std::size_t acoustic_steps = 0;
    std::size_t solid_steps = 0;
    double wall_seconds = 0.0;
    double compute_seconds = 0.0;
    double io_seconds = 0.0;
    double init_seconds = 0.0;
    double verification_seconds = 0.0;
    double time_per_outer_step = 0.0;
    std::size_t particle_updates = 0;
    std::string particle_update_definition = "unset";
    double particle_updates_per_second = 0.0;
    double mpps = 0.0;
    std::string neighbor_interactions;
    std::string mnips;
    std::string gpips;
    std::string peak_rss_kb = "unavailable";
    std::string peak_gpu_memory_kb = "unavailable";
    std::string sorting_seconds;
    std::string cll_seconds;
    std::string configuration_seconds;
    std::string advection_seconds;
    std::string acoustic_component_seconds;
    std::string pressure_seconds;
    std::string density_seconds;
    std::string viscous_seconds;
    std::string solid_component_seconds;
    std::string contact_seconds;
    std::string fsi_seconds;
    std::string diffusion_seconds;
    std::string reduction_seconds;
    std::string status = "completed";
    std::vector<std::pair<std::string, std::string>> extra_fields;

    void setParticleUpdates(std::size_t updates,
                            const std::string &definition)
    {
        particle_updates = updates;
        particle_update_definition = definition;
        if (compute_seconds > 0.0 && updates > 0)
        {
            particle_updates_per_second =
                static_cast<double>(updates) / compute_seconds;
            mpps = particle_updates_per_second / 1.0e6;
        }
        else
        {
            particle_updates_per_second = 0.0;
            mpps = 0.0;
        }
    }

    void captureHostMemory()
    {
        peak_rss_kb = peakRssKb();
    }
};

class BenchmarkRecorder
{
  public:
    BenchmarkRecorder(const std::string &case_name,
                      const BenchmarkConfig &config)
        : case_name_(safeComponent("case name", case_name)),
          run_id_(safeComponent("run id", config.run_id)),
          original_working_directory_(std::filesystem::current_path()),
          config_(config),
          git_commit_(environment("SPH_BENCH_GIT_COMMIT",
                                  SPHINXSYS_BENCHMARK_GIT_COMMIT)),
          build_(environment("SPH_BENCH_BUILD", SPHINXSYS_BENCHMARK_BUILD)),
          backend_(environment("SPH_BENCH_BACKEND", compiledBackend())),
          device_(environment("SPH_BENCH_DEVICE")),
          precision_(environment("SPH_BENCH_PRECISION", compiledPrecision()))
    {
        std::error_code error;
        const std::filesystem::path result_root =
            std::filesystem::absolute(config.result_directory);
        const std::filesystem::path case_directory =
            result_root / case_name_;
        std::filesystem::create_directories(case_directory, error);
        if (error)
        {
            throw std::runtime_error(
                "Unable to create benchmark result directory '" +
                case_directory.string() + "': " + error.message());
        }

        run_directory_ = case_directory / run_id_;
        if (!std::filesystem::create_directory(run_directory_, error))
        {
            const std::string reason =
                error ? error.message() : "run directory already exists";
            throw std::runtime_error("Refusing to overwrite benchmark run '" +
                                     run_directory_.string() + "': " + reason);
        }
        writeEnvironment();
    }

    const std::filesystem::path &runDirectory() const
    {
        return run_directory_;
    }

    const std::filesystem::path &originalWorkingDirectory() const
    {
        return original_working_directory_;
    }

    void activateRunDirectory() const
    {
        std::error_code error;
        std::filesystem::current_path(run_directory_, error);
        if (error)
        {
            throw std::runtime_error(
                "Unable to activate benchmark run directory '" +
                run_directory_.string() + "': " + error.message());
        }
    }

    void stageInputAssets() const
    {
        stageAssetDirectory("input");
        stageAssetDirectory("reload");
    }

    void writeSummary(BenchmarkSummary summary) const
    {
        if (summary.peak_rss_kb == "unavailable" ||
            summary.peak_rss_kb.empty())
        {
            summary.captureHostMemory();
        }

        const std::filesystem::path path = run_directory_ / "summary.csv";
        std::ofstream stream(path);
        if (!stream)
        {
            throw std::runtime_error("Unable to open '" + path.string() + "'");
        }

        stream << "case,run,git,build,precision,backend,device,dp,"
                  "initial_particle_count,particle_count,physical_time,"
                  "outer_steps,advection_steps,acoustic_steps,solid_steps,"
                  "wall_seconds,compute_seconds,io_seconds,init_seconds,"
                  "time_per_outer_step,status,requested_end_time,"
                  "benchmark_mode,output_enabled,output_interval,resolution,"
                  "verification_seconds,particle_updates,"
                  "particle_update_definition,particle_updates_per_second,"
                  "mpps,neighbor_interactions,mnips,gpips,peak_rss_kb,"
                  "peak_gpu_memory_kb,sorting_seconds,cll_seconds,"
                  "configuration_seconds,advection_seconds,"
                  "acoustic_component_seconds,pressure_seconds,"
                  "density_seconds,viscous_seconds,solid_component_seconds,"
                  "contact_seconds,fsi_seconds,diffusion_seconds,"
                  "reduction_seconds";
        for (const auto &field : summary.extra_fields)
        {
            stream << ',' << csv(field.first);
        }
        stream << '\n';
        stream << csv(case_name_) << ',' << csv(run_id_) << ','
               << csv(git_commit_) << ',' << csv(build_) << ','
               << csv(precision_) << ',' << csv(backend_) << ','
               << csv(device_) << ','
               << std::setprecision(std::numeric_limits<double>::max_digits10)
               << summary.dp << ',' << summary.initial_particle_count << ','
               << summary.particle_count << ',' << summary.physical_time << ','
               << summary.outer_steps << ',' << summary.advection_steps << ','
               << summary.acoustic_steps << ',' << summary.solid_steps << ','
               << summary.wall_seconds << ',' << summary.compute_seconds << ','
               << summary.io_seconds << ',' << summary.init_seconds << ','
               << summary.time_per_outer_step << ',' << csv(summary.status)
               << ',' << config_.end_time << ','
               << (config_.benchmark_mode ? "true" : "false") << ','
               << (config_.output_enabled ? "true" : "false") << ','
               << config_.output_interval << ',' << csv(config_.resolution)
               << ',' << summary.verification_seconds << ','
               << summary.particle_updates << ','
               << csv(summary.particle_update_definition) << ','
               << summary.particle_updates_per_second << ',' << summary.mpps
               << ',' << csv(summary.neighbor_interactions) << ','
               << csv(summary.mnips) << ',' << csv(summary.gpips) << ','
               << csv(summary.peak_rss_kb) << ','
               << csv(summary.peak_gpu_memory_kb) << ','
               << csv(summary.sorting_seconds) << ','
               << csv(summary.cll_seconds) << ','
               << csv(summary.configuration_seconds) << ','
               << csv(summary.advection_seconds) << ','
               << csv(summary.acoustic_component_seconds) << ','
               << csv(summary.pressure_seconds) << ','
               << csv(summary.density_seconds) << ','
               << csv(summary.viscous_seconds) << ','
               << csv(summary.solid_component_seconds) << ','
               << csv(summary.contact_seconds) << ','
               << csv(summary.fsi_seconds) << ','
               << csv(summary.diffusion_seconds) << ','
               << csv(summary.reduction_seconds);
        for (const auto &field : summary.extra_fields)
        {
            stream << ',' << csv(field.second);
        }
        stream << '\n';
        stream.close();
        if (!stream)
        {
            throw std::runtime_error("Unable to write '" + path.string() + "'");
        }
    }

  private:
    std::string case_name_;
    std::string run_id_;
    std::filesystem::path original_working_directory_;
    BenchmarkConfig config_;
    std::string git_commit_;
    std::string build_;
    std::string backend_;
    std::string device_;
    std::string precision_;
    std::filesystem::path run_directory_;

    static std::string environment(const char *name,
                                   const std::string &fallback = "unknown")
    {
        const char *value = std::getenv(name);
        return value != nullptr && *value != '\0' ? value : fallback;
    }

    static std::string compiledBackend()
    {
#if SPHINXSYS_USE_SYCL
        return "sycl";
#else
        return "host_tbb";
#endif
    }

    static std::string compiledPrecision()
    {
#if SPHINXSYS_USE_FLOAT
        return "float";
#else
        return "double";
#endif
    }

    static std::string safeComponent(const std::string &label,
                                     const std::string &value)
    {
        if (value.empty() || value == "." || value == ".." ||
            value.find('/') != std::string::npos ||
            value.find('\\') != std::string::npos)
        {
            throw std::invalid_argument(label + " must be one safe path component");
        }
        return value;
    }

    static std::string csv(const std::string &value)
    {
        if (value.find_first_of(",\"\r\n") == std::string::npos)
        {
            return value;
        }
        std::string escaped = "\"";
        for (const char character : value)
        {
            escaped += character == '"' ? "\"\"" : std::string(1, character);
        }
        return escaped + '"';
    }

    void stageAssetDirectory(const std::filesystem::path &name) const
    {
        const std::filesystem::path source =
            original_working_directory_ / name;
        const std::filesystem::path destination = run_directory_ / name;
        std::error_code error;
        if (!std::filesystem::exists(source, error))
        {
            if (error)
            {
                throw std::runtime_error(
                    "Unable to inspect benchmark input asset directory '" +
                    source.string() + "': " + error.message());
            }
            return;
        }
        if (!std::filesystem::is_directory(source, error))
        {
            if (error)
            {
                throw std::runtime_error(
                    "Unable to inspect benchmark input asset directory '" +
                    source.string() + "': " + error.message());
            }
            throw std::runtime_error(
                "Benchmark input asset path is not a directory: '" +
                source.string() + "'");
        }

        std::filesystem::create_directories(destination, error);
        if (error)
        {
            throw std::runtime_error(
                "Unable to create staged benchmark input asset directory '" +
                destination.string() + "': " + error.message());
        }
        std::filesystem::copy(
            source, destination,
            std::filesystem::copy_options::recursive |
                std::filesystem::copy_options::overwrite_existing,
            error);
        if (error)
        {
            throw std::runtime_error(
                "Unable to stage benchmark input assets from '" +
                source.string() + "' to '" + destination.string() +
                "': " + error.message());
        }
    }

    static std::string json(const std::string &value)
    {
        std::ostringstream escaped;
        for (const unsigned char character : value)
        {
            switch (character)
            {
            case '"':
                escaped << "\\\"";
                break;
            case '\\':
                escaped << "\\\\";
                break;
            case '\b':
                escaped << "\\b";
                break;
            case '\f':
                escaped << "\\f";
                break;
            case '\n':
                escaped << "\\n";
                break;
            case '\r':
                escaped << "\\r";
                break;
            case '\t':
                escaped << "\\t";
                break;
            default:
                if (character < 0x20)
                {
                    escaped << "\\u" << std::hex << std::setw(4)
                            << std::setfill('0')
                            << static_cast<unsigned int>(character)
                            << std::dec << std::setfill(' ');
                }
                else
                {
                    escaped << static_cast<char>(character);
                }
            }
        }
        return escaped.str();
    }

    void writeEnvironment() const
    {
        const std::filesystem::path path = run_directory_ / "environment.json";
        std::ofstream stream(path);
        if (!stream)
        {
            throw std::runtime_error("Unable to open '" + path.string() + "'");
        }
        stream << "{\n"
               << "  \"case\": \"" << json(case_name_) << "\",\n"
               << "  \"run\": \"" << json(run_id_) << "\",\n"
               << "  \"git\": \"" << json(git_commit_) << "\",\n"
               << "  \"build\": \"" << json(build_) << "\",\n"
               << "  \"backend\": \"" << json(backend_) << "\",\n"
               << "  \"device\": \"" << json(device_) << "\",\n"
               << "  \"precision\": \"" << json(precision_) << "\",\n"
               << "  \"host\": \"" << json(environment("SPH_BENCH_HOST"))
               << "\",\n"
               << "  \"os\": \"" << json(environment("SPH_BENCH_OS"))
               << "\",\n"
               << "  \"compiler\": \""
               << json(environment("SPH_BENCH_COMPILER",
                                   SPHINXSYS_BENCHMARK_COMPILER))
               << "\"\n"
               << "}\n";
        stream.close();
        if (!stream)
        {
            throw std::runtime_error("Unable to write '" + path.string() + "'");
        }
    }
};
} // namespace paper_bench

#endif // PAPER_BENCHMARK_RECORDER_H
