#ifndef PAPER_BENCHMARK_CONFIG_H
#define PAPER_BENCHMARK_CONFIG_H

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef SPHINXSYS_BENCHMARK_RESULT_ROOT
#define SPHINXSYS_BENCHMARK_RESULT_ROOT "benchmark_results"
#endif

namespace paper_bench
{
inline std::filesystem::path defaultResultDirectory()
{
    return std::filesystem::path(SPHINXSYS_BENCHMARK_RESULT_ROOT);
}
struct BenchmarkDefaults
{
    double dp = 0.05;
    double end_time = 20.0;
    double output_interval = 0.0;
    std::vector<std::pair<std::string, double>> resolutions{
        {"standard", 0.050}, {"s1", 0.050}, {"s2", 0.040},
        {"s3", 0.032},      {"s4", 0.025}, {"s5", 0.020},
        {"s6", 0.016}};
};

struct BenchmarkConfig
{
    double dp = 0.05;
    std::string resolution = "standard";
    double end_time = 20.0;
    bool output_enabled = true;
    double output_interval = 0.0;
    bool benchmark_mode = false;
    std::filesystem::path result_directory = defaultResultDirectory();
    std::string run_id;

    static BenchmarkConfig parse(int &argc, char *argv[])
    {
        return parse(argc, argv, BenchmarkDefaults{});
    }

    static BenchmarkConfig parse(int &argc, char *argv[],
                                 const BenchmarkDefaults &defaults)
    {
        BenchmarkConfig config;
        config.dp = defaults.dp;
        config.end_time = defaults.end_time;
        config.output_interval = defaults.output_interval;
        config.result_directory = defaultResultDirectory();
        bool dp_was_set = false;
        bool resolution_was_set = false;
        bool output_was_set = false;
        int output_index = 1;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument(argv[index]);
            const std::size_t equals = argument.find('=');
            const std::string option = argument.substr(
                0, equals == std::string::npos ? argument.size() : equals);
            if (!isBenchmarkOption(option))
            {
                argv[output_index++] = argv[index];
                continue;
            }

            if (option == "--benchmark" || option == "--output")
            {
                bool enabled = true;
                if (equals != std::string::npos)
                {
                    enabled = boolean(option, argument.substr(equals + 1));
                }
                else if (index + 1 < argc && isBoolean(argv[index + 1]))
                {
                    enabled = boolean(option, argv[++index]);
                }
                if (option == "--benchmark")
                {
                    config.benchmark_mode = enabled;
                }
                else
                {
                    config.output_enabled = enabled;
                    output_was_set = true;
                }
                continue;
            }

            const std::string value =
                requiredValue(option, argument, equals, index, argc, argv);
            if (option == "--dp")
            {
                config.dp = positiveDouble(option, value);
                dp_was_set = true;
            }
            else if (option == "--resolution")
            {
                config.resolution = resolutionLabel(value, defaults);
                resolution_was_set = true;
            }
            else if (option == "--end-time")
            {
                config.end_time = positiveDouble(option, value);
            }
            else if (option == "--output-interval")
            {
                config.output_interval = positiveDouble(option, value);
            }
            else if (option == "--result-dir")
            {
                if (value.empty())
                {
                    throw std::invalid_argument(option + " must not be empty");
                }
                config.result_directory = value;
            }
            else if (option == "--run-id")
            {
                validateRunId(value);
                config.run_id = value;
            }
        }
        argc = output_index;
        argv[output_index] = nullptr;

        if (!dp_was_set && resolution_was_set)
        {
            config.dp = resolutionDp(config.resolution, defaults);
        }
        if (config.benchmark_mode && !output_was_set)
        {
            config.output_enabled = false;
        }
        if (config.run_id.empty())
        {
            config.run_id = defaultRunId();
        }
        return config;
    }

    static double resolutionDp(const std::string &label)
    {
        return resolutionDp(label, BenchmarkDefaults{});
    }

    static double resolutionDp(const std::string &label,
                               const BenchmarkDefaults &defaults)
    {
        for (const auto &resolution : defaults.resolutions)
        {
            if (resolution.first == label)
            {
                return resolution.second;
            }
        }
        throw std::invalid_argument("--resolution has no dp mapping for '" +
                                    label + "'");
    }

  private:
    static bool isBenchmarkOption(const std::string &option)
    {
        return option == "--benchmark" || option == "--dp" ||
               option == "--resolution" || option == "--end-time" ||
               option == "--output" || option == "--output-interval" ||
               option == "--result-dir" || option == "--run-id";
    }

    static bool isBoolean(const std::string &value)
    {
        return value == "true" || value == "false" || value == "1" ||
               value == "0" || value == "on" || value == "off";
    }

    static std::string requiredValue(const std::string &option,
                                     const std::string &argument,
                                     std::size_t equals, int &index,
                                     int argc, char *argv[])
    {
        if (equals != std::string::npos)
        {
            const std::string value = argument.substr(equals + 1);
            if (!value.empty())
            {
                return value;
            }
        }
        else if (index + 1 < argc)
        {
            return argv[++index];
        }
        throw std::invalid_argument("Missing value for " + option);
    }

    static double positiveDouble(const std::string &option,
                                 const std::string &value)
    {
        char *end = nullptr;
        errno = 0;
        const double parsed = std::strtod(value.c_str(), &end);
        if (errno != 0 || end == value.c_str() || *end != '\0' ||
            !std::isfinite(parsed) || parsed <= 0.0)
        {
            throw std::invalid_argument(option + " must be a positive number");
        }
        return parsed;
    }

    static bool boolean(const std::string &option, const std::string &value)
    {
        if (value == "true" || value == "1" || value == "on")
        {
            return true;
        }
        if (value == "false" || value == "0" || value == "off")
        {
            return false;
        }
        throw std::invalid_argument(option + " must be true or false");
    }

    static std::string resolutionLabel(std::string value,
                                       const BenchmarkDefaults &defaults)
    {
        for (char &character : value)
        {
            if (character >= 'A' && character <= 'Z')
            {
                character = static_cast<char>(character - 'A' + 'a');
            }
        }
        for (const auto &resolution : defaults.resolutions)
        {
            if (resolution.first == value)
            {
                return value;
            }
        }
        throw std::invalid_argument("--resolution has no dp mapping for '" +
                                    value + "'");
    }

    static void validateRunId(const std::string &value)
    {
        if (value.empty() || value == "." || value == ".." ||
            value.find('/') != std::string::npos ||
            value.find('\\') != std::string::npos)
        {
            throw std::invalid_argument(
                "--run-id must be one safe path component");
        }
    }

    static std::string defaultRunId()
    {
        const std::time_t now = std::time(nullptr);
        std::tm local_time{};
#if defined(_WIN32)
        localtime_s(&local_time, &now);
#else
        localtime_r(&now, &local_time);
#endif
        char buffer[32]{};
        if (std::strftime(buffer, sizeof(buffer), "%Y%m%d-%H%M%S",
                          &local_time) == 0)
        {
            throw std::runtime_error("Unable to create benchmark run id");
        }
        return buffer;
    }
};
} // namespace paper_bench

#endif // PAPER_BENCHMARK_CONFIG_H
