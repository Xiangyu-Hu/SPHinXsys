#ifndef ELECTROMAGNETIC_OPHELIE_CLI_H
#define ELECTROMAGNETIC_OPHELIE_CLI_H

#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_racetrack_source.h"

#include <cstdlib>
#include <cstring>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline const char *phiSolverKindName(OpheliePhiSolverKind kind)
{
    switch (kind)
    {
    case OpheliePhiSolverKind::GMRES:
        return "GMRES";
    case OpheliePhiSolverKind::PCG:
        return "PCG";
    default:
        return "Jacobi";
    }
}

struct OphelieTestCliOptions
{
    bool compare_level0 = false;
    /** TEAM7 integration smoke: do not fail on phi residual / divJ (fields must still be finite and positive). */
    bool ophelie_smoke = false;
    /** Use TEAM7 STL geometry (particle_generation_em / COMSOL mm mesh, scaled to m). */
    bool native_stl = false;
    bool native_small_air_box = true;
    bool native_relax_air = true;
    Real native_dp_mm = 6.0;
    bool sigma_user_set = false;
    bool enable_self_induction = false;
    bool skip_relaxation = false;
    size_t relaxation_steps = 0;
    size_t relaxation_log_every = 50;
    /** Default 100 (official relax tests); --relax-vtp-every=0 disables relax VTP dumps. */
    size_t relaxation_vtp_every = 100;
    /** Non-empty: use this folder for Reload.xml (see test_3d_ophelie_team7/reload/README.md). */
    std::string reload_dir;
    bool no_power_scaling = false;
    bool compare_team7_bz = false;
    bool compare_team7_bz_loop = false;
    bool compare_team7_bz_rect_loop = false;
    std::string team7_reference_dir;
    Real team7_bz_rms_smoke_threshold = 0.5;
    Real team7_coil_turns = 0.0;
    OphelieCoilSourceModel coil_source_model = OphelieCoilSourceModel::VolumeETheta;
    Real racetrack_inset_mm = 0.0;
    Real racetrack_z_mm = 99.0;
    Real racetrack_ds_mm = 2.0;
    bool racetrack_sweep = false;
};

inline bool isOphelieTestCommandLineOption(const char *arg)
{
    return std::strcmp(arg, "--no-phi") == 0 || std::strcmp(arg, "--ophelie-compare-level0") == 0 ||
           std::strcmp(arg, "--ophelie-smoke") == 0 || std::strcmp(arg, "--native-stl") == 0 ||
           std::strcmp(arg, "--native-standard-air") == 0 || std::strcmp(arg, "--native-no-air-relax") == 0 ||
           std::strncmp(arg, "--native-dp-mm=", 15) == 0 || std::strcmp(arg, "--no-power-scaling") == 0 ||
           std::strncmp(arg, "--target-power=", 15) == 0 || std::strcmp(arg, "--compare-team7-bz") == 0 ||
           std::strcmp(arg, "--compare-team7-bz-loop") == 0 || std::strcmp(arg, "--compare-team7-bz-rect-loop") == 0 || std::strncmp(arg, "--team7-coil-turns=", 19) == 0 ||
           std::strcmp(arg, "--coil-j-outer-shell") == 0 || std::strncmp(arg, "--coil-j-outer-shell=", 21) == 0 ||
           std::strncmp(arg, "--team7-reference-dir=", 22) == 0 ||
           std::strncmp(arg, "--team7-bz-rms-threshold=", 25) == 0 ||
           std::strcmp(arg, "--self-induction") == 0 || std::strncmp(arg, "--phi-solver=", 13) == 0 ||
           std::strncmp(arg, "--self-induction-relax=", 23) == 0 ||
           std::strncmp(arg, "--self-induction-max-iter=", 26) == 0 || std::strcmp(arg, "--skip-relax") == 0 ||
           std::strncmp(arg, "--relax-steps=", 14) == 0 || std::strncmp(arg, "--relax-log-every=", 18) == 0 ||
           std::strncmp(arg, "--relax-vtp-every=", 18) == 0 || std::strncmp(arg, "--reload-dir=", 13) == 0 ||
           std::strncmp(arg, "--sigma=", 8) == 0 || std::strncmp(arg, "--coil-source-model=", 20) == 0 ||
           std::strncmp(arg, "--racetrack-inset-mm=", 21) == 0 || std::strncmp(arg, "--racetrack-z-mm=", 17) == 0 ||
           std::strncmp(arg, "--racetrack-ds-mm=", 18) == 0 || std::strcmp(arg, "--racetrack-sweep") == 0;
}

inline void applyOphelieTestCommandLineOption(const char *arg, OphelieParameters &params, OphelieTestCliOptions &cli_options)
{
    if (std::strcmp(arg, "--no-phi") == 0)
    {
        params.enable_phi_correction_ = false;
    }
    else if (std::strcmp(arg, "--ophelie-compare-level0") == 0)
    {
        cli_options.compare_level0 = true;
    }
    else if (std::strcmp(arg, "--ophelie-smoke") == 0)
    {
        cli_options.ophelie_smoke = true;
    }
    else if (std::strcmp(arg, "--native-stl") == 0)
    {
        cli_options.native_stl = true;
    }
    else if (std::strcmp(arg, "--native-standard-air") == 0)
    {
        cli_options.native_small_air_box = false;
    }
    else if (std::strcmp(arg, "--native-no-air-relax") == 0)
    {
        cli_options.native_relax_air = false;
    }
    else if (std::strncmp(arg, "--native-dp-mm=", 15) == 0)
    {
        cli_options.native_dp_mm = static_cast<Real>(std::atof(arg + 15));
    }
    else if (std::strcmp(arg, "--no-power-scaling") == 0)
    {
        cli_options.no_power_scaling = true;
        params.enable_power_scaling_ = false;
    }
    else if (std::strncmp(arg, "--target-power=", 15) == 0)
    {
        params.target_joule_power_ = static_cast<Real>(std::atof(arg + 15));
        if (params.target_joule_power_ <= TinyReal)
        {
            params.enable_power_scaling_ = false;
        }
    }
    else if (std::strcmp(arg, "--compare-team7-bz") == 0)
    {
        cli_options.compare_team7_bz = true;
    }
    else if (std::strcmp(arg, "--compare-team7-bz-loop") == 0)
    {
        cli_options.compare_team7_bz = true;
        cli_options.compare_team7_bz_loop = true;
    }
    else if (std::strcmp(arg, "--compare-team7-bz-rect-loop") == 0)
    {
        cli_options.compare_team7_bz = true;
        cli_options.compare_team7_bz_rect_loop = true;
        cli_options.coil_source_model = OphelieCoilSourceModel::FilamentRacetrack;
    }
    else if (std::strncmp(arg, "--coil-source-model=", 20) == 0)
    {
        cli_options.coil_source_model = parseOphelieCoilSourceModel(std::string(arg + 20));
        if (cli_options.coil_source_model == OphelieCoilSourceModel::FilamentRacetrack)
        {
            cli_options.compare_team7_bz = true;
        }
    }
    else if (std::strncmp(arg, "--racetrack-inset-mm=", 21) == 0)
    {
        cli_options.racetrack_inset_mm = static_cast<Real>(std::atof(arg + 21));
    }
    else if (std::strncmp(arg, "--racetrack-z-mm=", 17) == 0)
    {
        cli_options.racetrack_z_mm = static_cast<Real>(std::atof(arg + 17));
    }
    else if (std::strncmp(arg, "--racetrack-ds-mm=", 18) == 0)
    {
        cli_options.racetrack_ds_mm = static_cast<Real>(std::atof(arg + 18));
    }
    else if (std::strcmp(arg, "--racetrack-sweep") == 0)
    {
        cli_options.racetrack_sweep = true;
        cli_options.compare_team7_bz = true;
        cli_options.coil_source_model = OphelieCoilSourceModel::FilamentRacetrack;
    }
    else if (std::strncmp(arg, "--team7-coil-turns=", 19) == 0)
    {
        cli_options.team7_coil_turns = static_cast<Real>(std::atof(arg + 19));
    }
    else if (std::strcmp(arg, "--coil-j-outer-shell") == 0)
    {
        params.coil_j_outer_shell_only_ = true;
        params.coil_j_outer_shell_radius_fraction_ = 0.85;
    }
    else if (std::strncmp(arg, "--coil-j-outer-shell=", 21) == 0)
    {
        params.coil_j_outer_shell_only_ = true;
        params.coil_j_outer_shell_radius_fraction_ = static_cast<Real>(std::atof(arg + 21));
    }
    else if (std::strncmp(arg, "--team7-reference-dir=", 22) == 0)
    {
        cli_options.team7_reference_dir = std::string(arg + 22);
    }
    else if (std::strncmp(arg, "--team7-bz-rms-threshold=", 25) == 0)
    {
        cli_options.team7_bz_rms_smoke_threshold = static_cast<Real>(std::atof(arg + 25));
    }
    else if (std::strcmp(arg, "--self-induction") == 0)
    {
        cli_options.enable_self_induction = true;
        params.enable_self_induction_ = true;
        params.enable_phi_correction_ = true;
        std::cout << "[ophelie] WARNING: --self-induction is EXPERIMENTAL (phase convention under review). "
                     "Do not use for validation reports.\n";
    }
    else if (std::strncmp(arg, "--phi-solver=", 13) == 0)
    {
        const std::string solver_name(arg + 13);
        if (solver_name == "PCG")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::PCG;
        }
        else if (solver_name == "GMRES")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
        }
        else if (solver_name == "Jacobi")
        {
            params.phi_solver_kind_ = OpheliePhiSolverKind::Jacobi;
        }
    }
    else if (std::strncmp(arg, "--self-induction-relax=", 23) == 0)
    {
        params.self_induction_relaxation_factor_ = static_cast<Real>(std::atof(arg + 23));
    }
    else if (std::strncmp(arg, "--self-induction-max-iter=", 26) == 0)
    {
        params.self_induction_max_iterations_ = static_cast<size_t>(std::atoi(arg + 26));
    }
    else if (std::strcmp(arg, "--skip-relax") == 0)
    {
        cli_options.skip_relaxation = true;
    }
    else if (std::strncmp(arg, "--relax-steps=", 14) == 0)
    {
        cli_options.relaxation_steps = static_cast<size_t>(std::atoi(arg + 14));
    }
    else if (std::strncmp(arg, "--relax-log-every=", 18) == 0)
    {
        cli_options.relaxation_log_every = static_cast<size_t>(std::atoi(arg + 18));
    }
    else if (std::strncmp(arg, "--relax-vtp-every=", 18) == 0)
    {
        cli_options.relaxation_vtp_every = static_cast<size_t>(std::atoi(arg + 18));
    }
    else if (std::strncmp(arg, "--reload-dir=", 13) == 0)
    {
        cli_options.reload_dir = std::string(arg + 13);
    }
    else if (std::strncmp(arg, "--sigma=", 8) == 0)
    {
        params.sigma_glass_ = static_cast<Real>(std::atof(arg + 8));
        cli_options.sigma_user_set = true;
    }
}

inline StdVec<std::string> filterOphelieTestCommandLine(int ac, char *av[], OphelieParameters &params,
                                                        OphelieTestCliOptions &cli_options)
{
    StdVec<std::string> filtered_arguments;
    filtered_arguments.emplace_back(av[0]);
    for (int arg_index = 1; arg_index < ac; ++arg_index)
    {
        const char *arg = av[arg_index];
        if (isOphelieTestCommandLineOption(arg))
        {
            applyOphelieTestCommandLineOption(arg, params, cli_options);
            continue;
        }
        filtered_arguments.emplace_back(arg);
    }
    return filtered_arguments;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLI_H
