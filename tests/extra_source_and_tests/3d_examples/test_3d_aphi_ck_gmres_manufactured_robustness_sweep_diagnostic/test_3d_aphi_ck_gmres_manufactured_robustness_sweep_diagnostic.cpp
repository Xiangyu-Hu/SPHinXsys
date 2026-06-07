/**
 * Stage 9B: full A-phi manufactured-solution GMRES robustness sweep (informational).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_gmres_robustness_sweep_helpers.h"

#include <array>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

void printSweepCase(const char *sweep_name, const char *value_label, Real value,
                    const AphiManufacturedRobustnessResult &result)
{
    std::cout << " case=" << sweep_name << " " << value_label << "=" << value
              << " particles=" << result.total_real_particles
              << " init_norm=" << result.initial_residual_norm << " outer=" << result.outer_iteration_count
              << " arnoldi=" << result.arnoldi_step_count << " rel_res=" << result.final_relative_residual
              << " true_rel_res=" << result.final_true_relative_residual
              << " monotonic_outer=" << (result.monotonic_outer_residual ? 1 : 0)
              << " converged=" << (result.converged ? 1 : 0)
              << " breakdown_code=" << result.breakdown_code_name;
}

AphiManufacturedRobustnessResult runCase(int ac, char *av[], const AphiManufacturedRobustnessParams &params,
                                         Real tolerance, UnsignedInt restart_dimension,
                                         UnsignedInt max_outer_iterations)
{
    return runFullAphiManufacturedGMRES(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
}

} // namespace

int main(int ac, char *av[])
{
    const Real tolerance = 1.0e-5;
    const UnsignedInt restart_dimension = 30;
    const UnsignedInt max_outer_iterations = 30;

    const AphiManufacturedRobustnessParams baseline;
    const std::array<Real, 3> dp_values = {0.2, 0.1, 0.075};
    const std::array<Real, 3> omega_values = {0.1, 1.25, 10.0};
    const std::array<Real, 3> sigma_values = {0.0, 0.1, 2.0};
    const std::array<Real, 3> penalty_values = {1.0, 10.0, 100.0};

    UnsignedInt converged_count = 0;
    UnsignedInt total_count = 0;

    std::cout << "test_3d_aphi_ck_gmres_manufactured_robustness_sweep_diagnostic";

    {
        AphiManufacturedRobustnessParams params = baseline;
        const AphiManufacturedRobustnessResult result =
            runCase(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
        printSweepCase("baseline", "dp", params.dp_0, result);
        converged_count += result.converged ? 1 : 0;
        ++total_count;
    }

    for (const Real dp : dp_values)
    {
        AphiManufacturedRobustnessParams params = baseline;
        params.dp_0 = dp;
        const AphiManufacturedRobustnessResult result =
            runCase(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
        printSweepCase("dp", "dp", dp, result);
        converged_count += result.converged ? 1 : 0;
        ++total_count;
    }

    for (const Real omega : omega_values)
    {
        AphiManufacturedRobustnessParams params = baseline;
        params.omega = omega;
        const AphiManufacturedRobustnessResult result =
            runCase(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
        printSweepCase("omega", "omega", omega, result);
        converged_count += result.converged ? 1 : 0;
        ++total_count;
    }

    for (const Real sigma : sigma_values)
    {
        AphiManufacturedRobustnessParams params = baseline;
        params.sigma = sigma;
        const AphiManufacturedRobustnessResult result =
            runCase(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
        printSweepCase("sigma", "sigma", sigma, result);
        converged_count += result.converged ? 1 : 0;
        ++total_count;
    }

    for (const Real penalty : penalty_values)
    {
        AphiManufacturedRobustnessParams params = baseline;
        params.phi_gauge_penalty = penalty;
        const AphiManufacturedRobustnessResult result =
            runCase(ac, av, params, tolerance, restart_dimension, max_outer_iterations);
        printSweepCase("phi_gauge_penalty", "penalty", penalty, result);
        converged_count += result.converged ? 1 : 0;
        ++total_count;
    }

    std::cout << " converged_count=" << converged_count << " total_count=" << total_count << " passed=1" << std::endl;
    return 0;
}
