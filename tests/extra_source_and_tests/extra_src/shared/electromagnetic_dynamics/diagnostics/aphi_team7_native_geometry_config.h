#ifndef APHI_TEAM7_NATIVE_GEOMETRY_CONFIG_H
#define APHI_TEAM7_NATIVE_GEOMETRY_CONFIG_H

#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** SI geometry / discretization for native TEAM7 particle generation and reload tests. */
struct AphiTeam7NativeGeometryConfig
{
    /** Minimum particle spacing [m] (finest resolution near surfaces). */
    Real dp_reference = 0.006;
    Real stl_scale = 0.001;
    Vec3d air_box_lower{-0.050, -0.050, -0.050};
    Vec3d air_box_upper{0.350, 0.350, 0.200};
    std::string air_preset = "small";
    /**
     * Air-only dyadic coarsening away from coil/plate (0 = uniform air at dp_reference).
     * Coil/plate always use local_refinement_levels=0 (single layer at dp_reference).
     */
    int local_refinement_levels = 0;
    /** "uniform" | "multires" (empty = inferred from local_refinement_levels). */
    std::string discretization;
    /** If set, load/stage particles from ./reload_cases/<id>/ before solve. */
    std::string reload_case_id;
    /** Last successful `--team7-case=` preset name (empty if not set via CLI). */
    std::string team7_case_preset;
};

inline std::string team7NativeCliSuffix(const std::string &arg, const std::string &prefix)
{
    if (arg.rfind(prefix, 0) != 0)
    {
        return {};
    }
    return arg.substr(prefix.size());
}

/** Coil/plate: uniform lattice at dp_reference (no dyadic levels). */
inline Real team7NativeSolidAdaptiveGlobalResolution(const AphiTeam7NativeGeometryConfig &config)
{
    return config.dp_reference;
}

inline int team7NativeSolidRefinementLevels(const AphiTeam7NativeGeometryConfig &config)
{
    (void)config;
    return 0;
}

/** Air adaptive global_resolution / spacing_ref (coarsest air lattice spacing). */
inline Real team7NativeAirAdaptiveGlobalResolution(const AphiTeam7NativeGeometryConfig &config)
{
    return config.dp_reference *
           std::pow(2.0, static_cast<Real>(config.local_refinement_levels));
}

inline int team7NativeAirRefinementLevels(const AphiTeam7NativeGeometryConfig &config)
{
    return config.local_refinement_levels;
}

/** Domain SPHSystem GlobalResolution: finest spacing (coil/plate + air near conductors). */
inline Real team7NativeSphReferenceDp(const AphiTeam7NativeGeometryConfig &config)
{
    return config.dp_reference;
}

/** @deprecated alias; use team7NativeSphReferenceDp for SPHSystem ctor. */
inline Real team7NativeAdaptiveGlobalResolution(const AphiTeam7NativeGeometryConfig &config)
{
    return team7NativeSphReferenceDp(config);
}

inline bool team7NativeRefinementActive(const AphiTeam7NativeGeometryConfig &config)
{
    return config.local_refinement_levels > 0;
}

/** Finest SmoothingLengthRatio (= h_ref / h_local) on the most refined particles. */
inline Real team7NativeFinestSmoothingLengthRatio(const AphiTeam7NativeGeometryConfig &config)
{
    return std::pow(2.0, static_cast<Real>(config.local_refinement_levels));
}

inline std::string team7NativeDiscretizationLabel(const AphiTeam7NativeGeometryConfig &config)
{
    if (!config.discretization.empty())
    {
        return config.discretization;
    }
    return team7NativeRefinementActive(config) ? "multires" : "uniform";
}

inline void finalizeTeam7NativeGeometryConfig(AphiTeam7NativeGeometryConfig &config)
{
    if (config.local_refinement_levels < 0)
    {
        config.local_refinement_levels = 0;
    }
    if (config.discretization == "uniform")
    {
        config.local_refinement_levels = 0;
    }
    else if (config.discretization == "multires" && config.local_refinement_levels == 0)
    {
        std::cerr << "warning: TEAM7 multires discretization requires --team7-air-refinement-level=N > 0\n";
    }
    else if (config.discretization.empty())
    {
        config.discretization = team7NativeRefinementActive(config) ? "multires" : "uniform";
    }
}

inline std::string team7NativeReloadCasesRoot()
{
    return "./reload_cases";
}

inline std::string team7NativeWorkingReloadDirectory()
{
    return "./reload";
}

/** Canonical snapshot id, e.g. uniform_legacy_dp6 or multires_small_dp6_l2. */
inline std::string team7NativeReloadCaseId(const AphiTeam7NativeGeometryConfig &config)
{
    if (!config.reload_case_id.empty())
    {
        return config.reload_case_id;
    }
    const std::string disc = team7NativeDiscretizationLabel(config);
    const int dp_mm = static_cast<int>(std::lround(config.dp_reference * 1000.0));
    std::ostringstream id;
    id << disc << "_" << config.air_preset << "_dp" << dp_mm;
    if (disc == "multires" && config.local_refinement_levels > 0)
    {
        id << "_l" << config.local_refinement_levels;
    }
    return id.str();
}

inline std::string team7NativeReloadCaseDirectory(const std::string &case_id)
{
    return team7NativeReloadCasesRoot() + "/" + case_id;
}

inline bool copyTeam7NativeReloadArtifacts(const std::string &from_dir, const std::string &to_dir)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(to_dir, ec);
    const StdVec<std::string> files = {"Coil_rld.xml", "Plate_rld.xml", "Air_rld.xml", "team7_native_geometry.txt"};
    bool all_ok = true;
    for (const std::string &file : files)
    {
        const fs::path src = fs::path(from_dir) / file;
        const fs::path dst = fs::path(to_dir) / file;
        if (!fs::exists(src))
        {
            all_ok = false;
            continue;
        }
        fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
        all_ok = all_ok && !ec;
    }
    return all_ok;
}

inline bool publishTeam7NativeReloadCaseFromWorkingDirectory(const AphiTeam7NativeGeometryConfig &config)
{
    const std::string case_id = team7NativeReloadCaseId(config);
    return copyTeam7NativeReloadArtifacts(team7NativeWorkingReloadDirectory(),
                                          team7NativeReloadCaseDirectory(case_id));
}

/** Legacy snapshot folder names kept for existing ./reload_cases/ trees. */
inline StdVec<std::string> team7NativeReloadCaseLookupOrder(const std::string &case_id)
{
    StdVec<std::string> order;
    if (case_id.empty())
    {
        return order;
    }
    order.push_back(case_id);
    if (case_id == "uniform_legacy_dp6")
    {
        order.push_back("team7_native_uniform_legacy_dp6_l0");
    }
    else if (case_id == "team7_native_uniform_legacy_dp6_l0")
    {
        order.push_back("uniform_legacy_dp6");
    }
    else if (case_id == "multires_small_dp6_l2")
    {
        order.push_back("team7_native_multires_small_dp6_l2");
    }
    else if (case_id == "team7_native_multires_small_dp6_l2")
    {
        order.push_back("multires_small_dp6_l2");
    }
    else if (case_id == "multires_small_dp6_l1" || case_id == "team7_native_multires_small_dp6_l1")
    {
        if (case_id == "multires_small_dp6_l1")
        {
            order.push_back("team7_native_multires_small_dp6_l1");
        }
        else
        {
            order.push_back("multires_small_dp6_l1");
        }
    }
    return order;
}

inline bool stageTeam7NativeReloadCaseToWorkingDirectory(const std::string &case_id)
{
    for (const std::string &candidate : team7NativeReloadCaseLookupOrder(case_id))
    {
        if (copyTeam7NativeReloadArtifacts(team7NativeReloadCaseDirectory(candidate),
                                          team7NativeWorkingReloadDirectory()))
        {
            if (candidate != case_id)
            {
                std::cout << "team7_reload_case_staged_alias=" << candidate << " requested=" << case_id << "\n";
            }
            return true;
        }
    }
    return false;
}

inline std::string parseTeam7ReloadCaseIdFromCli(int ac, char *av[])
{
    static constexpr const char *k_prefix = "--team7-reload-case=";
    for (int i = 1; i < ac; ++i)
    {
        const std::string value = team7NativeCliSuffix(av[i], k_prefix);
        if (!value.empty())
        {
            return value;
        }
    }
    return {};
}

inline bool stageTeam7NativeReloadCaseFromCli(int ac, char *av[])
{
    const std::string case_id = parseTeam7ReloadCaseIdFromCli(ac, av);
    if (case_id.empty())
    {
        return true;
    }
    const bool ok = stageTeam7NativeReloadCaseToWorkingDirectory(case_id);
    if (!ok)
    {
        std::cerr << "error: failed to stage reload case \"" << case_id << "\" from "
                  << team7NativeReloadCaseDirectory(case_id) << " to " << team7NativeWorkingReloadDirectory()
                  << "\n";
    }
    else
    {
        std::cout << "team7_reload_case_staged=" << case_id << " from=" << team7NativeReloadCaseDirectory(case_id)
                  << "\n";
    }
    return ok;
}

inline AphiTeam7NativeGeometryConfig defaultTeam7NativeGeometryConfig()
{
    return AphiTeam7NativeGeometryConfig{};
}

inline void applyTeam7NativeAirPreset(const std::string &preset, AphiTeam7NativeGeometryConfig &config)
{
    if (preset == "small")
    {
        config.air_preset = "small";
        config.air_box_lower = Vec3d(-0.050, -0.050, -0.050);
        config.air_box_upper = Vec3d(0.350, 0.350, 0.200);
    }
    else if (preset == "big")
    {
        config.air_preset = "big";
        config.air_box_lower = Vec3d(-1.353, -1.353, -0.300);
        config.air_box_upper = Vec3d(1.647, 1.647, 0.449);
    }
    else if (preset == "legacy")
    {
        config.air_preset = "legacy";
        config.air_box_lower = Vec3d(-0.200, -0.200, -0.200);
        config.air_box_upper = Vec3d(0.500, 0.500, 0.300);
    }
}

inline bool applyTeam7NativeCasePreset(const std::string &preset, AphiTeam7NativeGeometryConfig &config)
{
    if (preset == "uniform-legacy-dp6")
    {
        config.dp_reference = 0.006;
        config.discretization = "uniform";
        config.local_refinement_levels = 0;
        applyTeam7NativeAirPreset("legacy", config);
        return true;
    }
    if (preset == "multires-small-dp6-l1")
    {
        config.dp_reference = 0.006;
        config.discretization = "multires";
        config.local_refinement_levels = 1;
        applyTeam7NativeAirPreset("small", config);
        return true;
    }
    if (preset == "multires-small-dp6-l2")
    {
        config.dp_reference = 0.006;
        config.discretization = "multires";
        config.local_refinement_levels = 2;
        applyTeam7NativeAirPreset("small", config);
        return true;
    }
    if (preset == "uniform-small-dp6")
    {
        config.dp_reference = 0.006;
        config.discretization = "uniform";
        config.local_refinement_levels = 0;
        applyTeam7NativeAirPreset("small", config);
        return true;
    }
    if (preset == "uniform")
    {
        config.discretization = "uniform";
        config.local_refinement_levels = 0;
        return true;
    }
    if (preset == "multires")
    {
        config.discretization = "multires";
        return true;
    }
    return false;
}

inline Real parseTeam7EnvReal(const char *name, Real fallback)
{
    if (const char *value = std::getenv(name))
    {
        return static_cast<Real>(std::atof(value));
    }
    return fallback;
}

inline std::string parseTeam7EnvString(const char *name, const std::string &fallback)
{
    if (const char *value = std::getenv(name))
    {
        return std::string(value);
    }
    return fallback;
}

inline Vec3d parseTeam7CommaVec3(const std::string &text, const Vec3d &fallback)
{
    std::stringstream stream(text);
    std::string item;
    StdVec<Real> values;
    while (std::getline(stream, item, ',') && values.size() < 3)
    {
        values.push_back(static_cast<Real>(std::atof(item.c_str())));
    }
    if (values.size() != 3)
    {
        return fallback;
    }
    return Vec3d(values[0], values[1], values[2]);
}

inline void applyTeam7NativeGeometryEnvironment(AphiTeam7NativeGeometryConfig &config)
{
    const std::string preset = parseTeam7EnvString("TEAM7_AIR_PRESET", config.air_preset);
    if (!preset.empty())
    {
        applyTeam7NativeAirPreset(preset, config);
    }
    config.dp_reference = parseTeam7EnvReal("TEAM7_DP_REFERENCE",
                                            parseTeam7EnvReal("TEAM7_DP_0", config.dp_reference));
    config.stl_scale = parseTeam7EnvReal("TEAM7_STL_SCALE", config.stl_scale);
    config.air_box_lower[0] = parseTeam7EnvReal("TEAM7_AIR_BOX_LOWER_X", config.air_box_lower[0]);
    config.air_box_lower[1] = parseTeam7EnvReal("TEAM7_AIR_BOX_LOWER_Y", config.air_box_lower[1]);
    config.air_box_lower[2] = parseTeam7EnvReal("TEAM7_AIR_BOX_LOWER_Z", config.air_box_lower[2]);
    config.air_box_upper[0] = parseTeam7EnvReal("TEAM7_AIR_BOX_UPPER_X", config.air_box_upper[0]);
    config.air_box_upper[1] = parseTeam7EnvReal("TEAM7_AIR_BOX_UPPER_Y", config.air_box_upper[1]);
    config.air_box_upper[2] = parseTeam7EnvReal("TEAM7_AIR_BOX_UPPER_Z", config.air_box_upper[2]);
    config.local_refinement_levels = static_cast<int>(parseTeam7EnvReal(
        "TEAM7_LOCAL_REFINEMENT_LEVEL",
        parseTeam7EnvReal("TEAM7_AIR_MR_LEVELS", static_cast<Real>(config.local_refinement_levels))));
    if (parseTeam7EnvString("TEAM7_AIR_MULTIRES", "0") == "1" && config.local_refinement_levels == 0)
    {
        config.local_refinement_levels = 1;
    }
}

inline void printTeam7NativeGeometryConfigHelp()
{
    std::cout << "TEAM7 native geometry options (SI meters):\n"
              << "  --team7-case=<preset>       uniform-legacy-dp6 | uniform-small-dp6 | multires-small-dp6-l1 | multires-small-dp6-l2 | uniform | multires\n"
              << "  --team7-discretization=uniform|multires  uniform forces level 0; multires needs air refinement level\n"
              << "  --team7-reload-case=<id>    stage ./reload_cases/<id>/ into ./reload/ before solve\n"
              << "  --team7-dp=<m>              minimum particle spacing dp_reference (default 0.006)\n"
              << "  --team7-air=<preset>        small | big | legacy\n"
              << "  --team7-air-lower=x,y,z     explicit air box lower [m]\n"
              << "  --team7-air-upper=x,y,z     explicit air box upper [m]\n"
              << "  --team7-write-vtp           write ./output/*_ite_0000000000.vtp after EM solve (P3 etc.)\n"
              << "  --team7-local-refinement-level=N  air only: 0=uniform; N>0 coarsens air by 2^N (coil/plate stay uniform)\n"
              << "  --team7-air-refinement-level=N     alias (preferred name for air-only MR)\n"
              << "  --team7-refinement-level=N         alias for --team7-local-refinement-level\n"
              << "  --team7-air-mr              alias: set local_refinement_level=1 if still 0\n"
              << "  --team7-air-mr-levels=N     alias for --team7-local-refinement-level\n"
              << "  --team7-help                print this help\n"
              << "Environment: TEAM7_DP_REFERENCE, TEAM7_DP_0, TEAM7_AIR_PRESET, TEAM7_AIR_BOX_* , TEAM7_WRITE_VTP=1\n"
              << "Air presets (SI [m], outer box minus coil/plate STL):\n"
              << "  small   (-0.05,-0.05,-0.05) .. (0.35,0.35,0.20)   ~0.04 m^3  fast dev\n"
              << "  legacy  (-0.20,-0.20,-0.20) .. (0.50,0.50,0.30)   ~0.25 m^3  slightly larger than small\n"
              << "  big     (-1.353,-1.353,-0.30) .. (1.647,1.647,0.449) ~6.7 m^3  COMSOL TEAM7 far field\n"
              << "Reload snapshots: ./reload_cases/uniform_legacy_dp6/ | ./reload_cases/multires_small_dp6_l1/ | ./reload_cases/multires_small_dp6_l2/\n";
}

inline bool team7NativeCliRequestsWriteVtp(int ac, char *av[])
{
    if (parseTeam7EnvString("TEAM7_WRITE_VTP", "0") == "1")
    {
        return true;
    }
    for (int i = 1; i < ac; ++i)
    {
        if (std::string(av[i]) == "--team7-write-vtp")
        {
            return true;
        }
    }
    return false;
}

inline bool parseTeam7LegacyAirRefinementArg(const std::string &arg, int &out_levels)
{
    const StdVec<std::string> prefixes = {"air_refinement_levels=", "air_local_refinement_levels="};
    for (const std::string &prefix : prefixes)
    {
        if (arg.rfind(prefix, 0) == 0)
        {
            out_levels = std::atoi(arg.substr(prefix.size()).c_str());
            return true;
        }
    }
    return false;
}

inline bool parseTeam7NativeGeometryCli(int ac, char *av[], AphiTeam7NativeGeometryConfig &config)
{
    bool show_help = false;
    for (int i = 1; i < ac; ++i)
    {
        const std::string arg = av[i];
        int legacy_air_levels = 0;
        if (parseTeam7LegacyAirRefinementArg(arg, legacy_air_levels))
        {
            std::cerr << "warning: \"" << arg
                      << "\" is not parsed (missing --team7- prefix). Use --team7-air-refinement-level="
                      << legacy_air_levels << " or --team7-local-refinement-level=" << legacy_air_levels << "\n";
            continue;
        }
        if (arg == "--team7-help" || arg == "--help-team7")
        {
            show_help = true;
        }
        else if (const std::string preset = team7NativeCliSuffix(arg, "--team7-case="); !preset.empty())
        {
            if (applyTeam7NativeCasePreset(preset, config))
            {
                config.team7_case_preset = preset;
            }
            else
            {
                std::cerr << "warning: unknown --team7-case preset: " << preset << "\n";
            }
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-discretization="); !value.empty())
        {
            config.discretization = value;
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-reload-case="); !value.empty())
        {
            config.reload_case_id = value;
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-dp="); !value.empty())
        {
            config.dp_reference = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-air="); !value.empty())
        {
            applyTeam7NativeAirPreset(value, config);
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-air-lower="); !value.empty())
        {
            config.air_box_lower = parseTeam7CommaVec3(value, config.air_box_lower);
            config.air_preset = "custom";
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-air-upper="); !value.empty())
        {
            config.air_box_upper = parseTeam7CommaVec3(value, config.air_box_upper);
            config.air_preset = "custom";
        }
        else if (arg == "--team7-air-mr")
        {
            if (config.local_refinement_levels == 0)
            {
                config.local_refinement_levels = 1;
            }
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-local-refinement-level="); !value.empty())
        {
            config.local_refinement_levels = std::atoi(value.c_str());
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-air-refinement-level="); !value.empty())
        {
            config.local_refinement_levels = std::atoi(value.c_str());
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-refinement-level="); !value.empty())
        {
            config.local_refinement_levels = std::atoi(value.c_str());
        }
        else if (const std::string value = team7NativeCliSuffix(arg, "--team7-air-mr-levels="); !value.empty())
        {
            config.local_refinement_levels = std::atoi(value.c_str());
        }
        else if (arg.rfind("--team7-dp-air-finest=", 0) == 0)
        {
            std::cerr << "warning: --team7-dp-air-finest is deprecated; use --team7-dp and --team7-local-refinement-level\n";
        }
    }
    finalizeTeam7NativeGeometryConfig(config);
    if (show_help)
    {
        printTeam7NativeGeometryConfigHelp();
    }
    return show_help;
}

inline bool writeTeam7NativeGeometryMetadata(const std::string &path, const AphiTeam7NativeGeometryConfig &config)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "dp_reference_m," << config.dp_reference << "\n";
    output << "dp_0_m," << config.dp_reference << "\n";
    output << "stl_scale," << config.stl_scale << "\n";
    output << "air_preset," << config.air_preset << "\n";
    output << "air_lower_x_m," << config.air_box_lower[0] << "\n";
    output << "air_lower_y_m," << config.air_box_lower[1] << "\n";
    output << "air_lower_z_m," << config.air_box_lower[2] << "\n";
    output << "air_upper_x_m," << config.air_box_upper[0] << "\n";
    output << "air_upper_y_m," << config.air_box_upper[1] << "\n";
    output << "air_upper_z_m," << config.air_box_upper[2] << "\n";
    output << "local_refinement_levels," << config.local_refinement_levels << "\n";
    output << "air_local_refinement_levels," << config.local_refinement_levels << "\n";
    output << "air_refinement_levels," << config.local_refinement_levels << "\n";
    output << "air_multiresolution," << (team7NativeRefinementActive(config) ? 1 : 0) << "\n";
    output << "dp_air_coarsest_m," << team7NativeAirAdaptiveGlobalResolution(config) << "\n";
    output << "solid_refinement_levels,0\n";
    output << "dp_air_finest_m," << config.dp_reference << "\n";
    output << "sph_reference_dp_m," << team7NativeSphReferenceDp(config) << "\n";
    output << "sph_system_resolution_m," << team7NativeSphReferenceDp(config) << "\n";
    output << "discretization," << team7NativeDiscretizationLabel(config) << "\n";
    output << "reload_case_id," << team7NativeReloadCaseId(config) << "\n";
    return true;
}

inline bool loadTeam7NativeGeometryMetadata(const std::string &path, AphiTeam7NativeGeometryConfig &config)
{
    std::ifstream input(path);
    if (!input)
    {
        return false;
    }
    std::string line;
    while (std::getline(input, line))
    {
        const size_t comma = line.find(',');
        if (comma == std::string::npos)
        {
            continue;
        }
        const std::string key = line.substr(0, comma);
        const std::string value = line.substr(comma + 1);
        if (key == "dp_reference_m" || key == "dp_0_m")
        {
            config.dp_reference = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "stl_scale")
        {
            config.stl_scale = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_preset")
        {
            config.air_preset = value;
        }
        else if (key == "air_lower_x_m")
        {
            config.air_box_lower[0] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_lower_y_m")
        {
            config.air_box_lower[1] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_lower_z_m")
        {
            config.air_box_lower[2] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_upper_x_m")
        {
            config.air_box_upper[0] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_upper_y_m")
        {
            config.air_box_upper[1] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_upper_z_m")
        {
            config.air_box_upper[2] = static_cast<Real>(std::atof(value.c_str()));
        }
        else if (key == "air_multiresolution")
        {
            if (std::atoi(value.c_str()) != 0 && config.local_refinement_levels == 0)
            {
                config.local_refinement_levels = 1;
            }
        }
        else if (key == "dp_air_finest_m")
        {
            /** Legacy field; dp_reference_m is authoritative (do not override uniform dp=6 mm). */
            (void)value;
        }
        else if (key == "local_refinement_levels" || key == "air_local_refinement_levels" ||
                 key == "air_refinement_levels")
        {
            config.local_refinement_levels = std::atoi(value.c_str());
        }
        else if (key == "discretization")
        {
            config.discretization = value;
        }
        else if (key == "reload_case_id")
        {
            config.reload_case_id = value;
        }
    }
    finalizeTeam7NativeGeometryConfig(config);
    return true;
}

inline bool team7NativeReloadCaseIdsCompatible(const std::string &requested, const std::string &loaded)
{
    if (requested.empty() || loaded.empty())
    {
        return true;
    }
    if (requested == loaded)
    {
        return true;
    }
    for (const std::string &a : team7NativeReloadCaseLookupOrder(requested))
    {
        for (const std::string &b : team7NativeReloadCaseLookupOrder(loaded))
        {
            if (a == b)
            {
                return true;
            }
        }
    }
    return false;
}

inline bool team7NativeAuditReloadMetadata(const AphiTeam7NativeGeometryConfig &resolved,
                                           const std::string &requested_reload_case_id,
                                           const std::string &metadata_path = "./reload/team7_native_geometry.txt")
{
    AphiTeam7NativeGeometryConfig from_metadata;
    if (!loadTeam7NativeGeometryMetadata(metadata_path, from_metadata))
    {
        if (!requested_reload_case_id.empty())
        {
            std::cerr << "warning: --team7-reload-case=" << requested_reload_case_id << " but no metadata at "
                      << metadata_path << "\n";
        }
        return true;
    }
    std::cout << "team7_reload_metadata: reload_case_id=" << from_metadata.reload_case_id
              << " dp_reference_m=" << from_metadata.dp_reference << " air_preset=" << from_metadata.air_preset
              << " local_refinement_levels=" << from_metadata.local_refinement_levels
              << " discretization=" << team7NativeDiscretizationLabel(from_metadata) << "\n";
    if (!requested_reload_case_id.empty() && !from_metadata.reload_case_id.empty() &&
        !team7NativeReloadCaseIdsCompatible(requested_reload_case_id, from_metadata.reload_case_id))
    {
        std::cerr << "error: requested reload case \"" << requested_reload_case_id
                  << "\" does not match metadata reload_case_id=\"" << from_metadata.reload_case_id << "\"\n";
        return false;
    }
    if (!requested_reload_case_id.empty())
    {
        constexpr Real tol = static_cast<Real>(1.0e-9);
        if (std::abs(resolved.dp_reference - from_metadata.dp_reference) > tol)
        {
            std::cout << "note: dp_reference differs from reload metadata: resolved=" << resolved.dp_reference
                      << " metadata=" << from_metadata.dp_reference << "\n";
        }
        if (resolved.air_preset != from_metadata.air_preset && resolved.air_preset != "custom")
        {
            std::cout << "note: air_preset differs from reload metadata: resolved=" << resolved.air_preset
                      << " metadata=" << from_metadata.air_preset << "\n";
        }
        if (resolved.local_refinement_levels != from_metadata.local_refinement_levels)
        {
            std::cout << "note: local_refinement_levels differs from reload metadata: resolved="
                      << resolved.local_refinement_levels << " metadata=" << from_metadata.local_refinement_levels
                      << "\n";
        }
    }
    return true;
}

inline void filterSphinxsysCommandlineArgs(int ac, char *av[], StdVec<std::string> &storage, StdVec<char *> &filtered)
{
    storage.clear();
    filtered.clear();
    if (ac <= 0 || av == nullptr)
    {
        return;
    }
    storage.emplace_back(av[0]);
    filtered.push_back(storage.back().data());
    for (int i = 1; i < ac; ++i)
    {
        const std::string arg = av[i];
        if (arg.rfind("--team7-", 0) == 0)
        {
            continue;
        }
        storage.push_back(arg);
        filtered.push_back(storage.back().data());
    }
}

/** Defaults -> reload metadata (if present) -> environment -> CLI overrides. */
inline AphiTeam7NativeGeometryConfig resolveTeam7NativeGeometryConfig(int ac, char *av[],
                                                                      const std::string &metadata_path =
                                                                          "./reload/team7_native_geometry.txt")
{
    const std::string requested_reload_case_id = parseTeam7ReloadCaseIdFromCli(ac, av);
    if (!stageTeam7NativeReloadCaseFromCli(ac, av))
    {
        std::cerr << "error: --team7-reload-case staging failed; ensure ./reload_cases/<id>/ exists\n";
        std::exit(1);
    }
    AphiTeam7NativeGeometryConfig config = defaultTeam7NativeGeometryConfig();
    loadTeam7NativeGeometryMetadata(metadata_path, config);
    applyTeam7NativeGeometryEnvironment(config);
    parseTeam7NativeGeometryCli(ac, av, config);
    finalizeTeam7NativeGeometryConfig(config);
    if (!team7NativeAuditReloadMetadata(config, requested_reload_case_id, metadata_path))
    {
        std::exit(1);
    }
    return config;
}

inline void printTeam7NativeGeometryConfig(const AphiTeam7NativeGeometryConfig &config)
{
    std::cout << "TEAM7 native geometry config:"
              << " team7_case=" << (config.team7_case_preset.empty() ? "(none)" : config.team7_case_preset)
              << " discretization=" << team7NativeDiscretizationLabel(config)
              << " reload_case_id=" << team7NativeReloadCaseId(config)
              << " dp_reference_m=" << config.dp_reference << " stl_scale=" << config.stl_scale
              << " air_preset=" << config.air_preset << " air_lower_m=(" << config.air_box_lower.transpose()
              << ") air_upper_m=(" << config.air_box_upper.transpose() << ")"
              << " air_local_refinement_levels=" << config.local_refinement_levels
              << " coil_plate_refinement_levels=0"
              << " air_coarsest_spacing_m=" << team7NativeAirAdaptiveGlobalResolution(config)
              << " sph_reference_dp_m=" << team7NativeSphReferenceDp(config) << "\n";
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_GEOMETRY_CONFIG_H
