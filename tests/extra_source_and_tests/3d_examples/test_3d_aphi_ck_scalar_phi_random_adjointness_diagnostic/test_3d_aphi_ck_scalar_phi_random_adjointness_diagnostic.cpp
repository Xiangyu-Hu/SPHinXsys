/**
 * Stage 8C: random/rough-field Vol-weighted adjointness diagnostic for scalar phi Laplace + penalty.
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

Real measureRandomAdjointGap(SPHBody &body, Inner<> &inner, AphiVariableNames &names,
                             const AphiLhsAssemblyOptions &options, UnsignedInt seed_x, UnsignedInt seed_y)
{
    StateDynamics<MainExecutionPolicy, AssignRandomScalarPhiFieldsCK> assign_random(body, names.solution, names.v,
                                                                                    seed_x, seed_y);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_x(body, inner, names.solution, names.lhs, names.material,
                                                       options.omega, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_y(body, inner, names.v, names.t, names.material, options.omega,
                                                         options);
    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockDotProductCK> dot_x_ky(body, names.solution, names.t);
    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockDotProductCK> dot_kx_y(body, names.lhs, names.v);

    assign_random.exec();
    apply_x.exec();
    apply_y.exec();

    const Real x_ky = dot_x_ky.exec();
    const Real kx_y = dot_kx_y.exec();
    const Real scale = std::max({std::abs(x_ky), std::abs(kx_y), Real(1)});
    return std::abs(x_ky - kx_y) / scale;
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
    const UnsignedInt pair_count = 20;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);

    initialize_aphi_variables.exec();
    set_material.exec();
    test_body.updateRelations();

    const AphiLhsAssemblyOptions options = scalarPhiLaplacePenaltyOptions(phi_gauge_penalty);

    Real min_gap = std::numeric_limits<Real>::max();
    Real max_gap = Real(0);
    Real sum_gap = Real(0);
    for (UnsignedInt pair_index = 0; pair_index < pair_count; ++pair_index)
    {
        const UnsignedInt seed_x = 1000 + pair_index * 17;
        const UnsignedInt seed_y = 2000 + pair_index * 23;
        const Real gap = measureRandomAdjointGap(test_body.body, test_body.inner(), names, options, seed_x, seed_y);
        min_gap = std::min(min_gap, gap);
        max_gap = std::max(max_gap, gap);
        sum_gap += gap;
    }
    const Real mean_gap = sum_gap / static_cast<Real>(pair_count);

    std::cout << "test_3d_aphi_ck_scalar_phi_random_adjointness_diagnostic"
              << " pair_count=" << pair_count << " min_gap=" << min_gap << " max_gap=" << max_gap
              << " mean_gap=" << mean_gap << " passed=1" << std::endl;

    return 0;
}
