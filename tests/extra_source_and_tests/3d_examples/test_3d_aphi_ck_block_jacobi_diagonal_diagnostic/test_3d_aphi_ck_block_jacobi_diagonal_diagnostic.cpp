/**
 * Stage 8A: block-Jacobi diagonal min/max/non-positive count diagnostic.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <algorithm>
#include <iostream>
#include <limits>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

struct RealArrayStats
{
    Real min_value = std::numeric_limits<Real>::infinity();
    Real max_value = -std::numeric_limits<Real>::infinity();
    size_t non_positive_count = 0;
};

RealArrayStats hostRealArrayStats(BaseParticles &particles, const std::string &name, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, name);
    const Real *data = particles.getVariableDataByName<Real>(name);
    RealArrayStats stats;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        stats.min_value = std::min(stats.min_value, data[i]);
        stats.max_value = std::max(stats.max_value, data[i]);
        if (data[i] <= Real(0))
        {
            ++stats.non_positive_count;
        }
    }
    return stats;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real phi_gauge_penalty = 10.0;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);

    AphiLhsAssemblyOptions options;
    options.terms.laplace_a = true;
    options.terms.laplace_phi = true;
    options.terms.reaction = true;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal(
        DynamicsArgs(test_body.inner(), names.material, options.omega, options));

    initialize_aphi_variables.exec();
    set_material.exec();
    test_body.updateRelations();
    compute_jacobi_diagonal.exec();

    const size_t total_real_particles = test_body.body.getBaseParticles().TotalRealParticles();
    const AphiBlockJacobiDiagonalNames diag_names;
    const RealArrayStats laplace_a_stats =
        hostRealArrayStats(test_body.body.getBaseParticles(), diag_names.laplace_a_diag, total_real_particles);
    const RealArrayStats laplace_phi_stats =
        hostRealArrayStats(test_body.body.getBaseParticles(), diag_names.laplace_phi_diag, total_real_particles);

    Real d_phi_min = std::numeric_limits<Real>::infinity();
    Real d_phi_max = -std::numeric_limits<Real>::infinity();
    size_t d_phi_non_positive_count = 0;
    syncVariableToHost<Real>(test_body.body.getBaseParticles(), diag_names.laplace_phi_diag);
    const Real *laplace_phi_diag =
        test_body.body.getBaseParticles().getVariableDataByName<Real>(diag_names.laplace_phi_diag);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        const Real d_phi = laplace_phi_diag[i] + phi_gauge_penalty;
        d_phi_min = std::min(d_phi_min, d_phi);
        d_phi_max = std::max(d_phi_max, d_phi);
        if (d_phi <= Real(0))
        {
            ++d_phi_non_positive_count;
        }
    }

    const bool passed = std::isfinite(laplace_a_stats.min_value) && std::isfinite(laplace_a_stats.max_value) &&
                        std::isfinite(laplace_phi_stats.min_value) && std::isfinite(laplace_phi_stats.max_value);

    std::cout << "test_3d_aphi_ck_block_jacobi_diagonal_diagnostic"
              << " laplace_a_diag_min=" << laplace_a_stats.min_value
              << " laplace_a_diag_max=" << laplace_a_stats.max_value
              << " laplace_a_diag_non_positive_count=" << laplace_a_stats.non_positive_count
              << " laplace_phi_diag_min=" << laplace_phi_stats.min_value
              << " laplace_phi_diag_max=" << laplace_phi_stats.max_value
              << " laplace_phi_diag_non_positive_count=" << laplace_phi_stats.non_positive_count
              << " phi_penalty=" << phi_gauge_penalty << " d_phi_min=" << d_phi_min << " d_phi_max=" << d_phi_max
              << " d_phi_non_positive_count=" << d_phi_non_positive_count << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
