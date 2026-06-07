/**
 * Stage 8C: random Rayleigh/curvature diagnostic for scalar phi Laplace + penalty.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_scalar_phi_diagnostic_helpers.h"

#include <iostream>
#include <limits>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

struct CurvatureStats
{
    Real min_rayleigh = std::numeric_limits<Real>::max();
    Real max_rayleigh = std::numeric_limits<Real>::lowest();
    UnsignedInt negative_count = 0;
    UnsignedInt near_zero_count = 0;
};

CurvatureStats measureCurvatureStats(SPHBody &body, Inner<> &inner, AphiVariableNames &names,
                                     const AphiLhsAssemblyOptions &options, BaseParticles &particles,
                                     size_t total_real_particles, UnsignedInt seed, Real near_zero_tol)
{
    StateDynamics<MainExecutionPolicy, AssignRandomScalarPhiSingleBlockCK> assign_random(body, names.solution, seed);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_operator(body, inner, names.solution, names.lhs, names.material,
                                                                options.omega, options);
    assign_random.exec();
    apply_operator.exec();

    const Real pq = hostBlockDotProduct(particles, names.solution, names.lhs, total_real_particles);
    const Real pp = hostBlockDotProduct(particles, names.solution, names.solution, total_real_particles);
    const Real rayleigh = pq / (pp + TinyReal);

    CurvatureStats stats;
    stats.min_rayleigh = rayleigh;
    stats.max_rayleigh = rayleigh;
    if (rayleigh < Real(0))
    {
        ++stats.negative_count;
    }
    if (std::abs(rayleigh) <= near_zero_tol)
    {
        ++stats.near_zero_count;
    }
    return stats;
}

void mergeStats(CurvatureStats &aggregate, const CurvatureStats &sample)
{
    aggregate.min_rayleigh = std::min(aggregate.min_rayleigh, sample.min_rayleigh);
    aggregate.max_rayleigh = std::max(aggregate.max_rayleigh, sample.max_rayleigh);
    aggregate.negative_count += sample.negative_count;
    aggregate.near_zero_count += sample.near_zero_count;
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
    const Real near_zero_tol = 1.0e-12;
    const UnsignedInt seed_count = 20;
    const Real penalty_values[4] = {Real(0), Real(10), Real(100), Real(1000)};

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);

    initialize_aphi_variables.exec();
    set_material.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();

    std::cout << "test_3d_aphi_ck_scalar_phi_random_curvature_diagnostic"
              << " total_real_particles=" << total_real_particles;

    for (const Real penalty : penalty_values)
    {
        const AphiLhsAssemblyOptions options = scalarPhiLaplacePenaltyOptions(penalty, penalty > Real(0));
        CurvatureStats aggregate;
        aggregate.min_rayleigh = std::numeric_limits<Real>::max();
        aggregate.max_rayleigh = std::numeric_limits<Real>::lowest();

        for (UnsignedInt seed = 1; seed <= seed_count; ++seed)
        {
            mergeStats(aggregate, measureCurvatureStats(test_body.body, test_body.inner(), names, options, particles,
                                                        total_real_particles, seed, near_zero_tol));
        }

        std::cout << " penalty=" << penalty << " min_rayleigh=" << aggregate.min_rayleigh
                  << " max_rayleigh=" << aggregate.max_rayleigh << " negative_count=" << aggregate.negative_count
                  << " near_zero_count=" << aggregate.near_zero_count;
    }

    std::cout << " passed=1" << std::endl;
    return 0;
}
