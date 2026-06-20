/**
 * @file test_3d_ophelie_phi_solve_manufactured_A.cpp
 * @brief MMS: phi_exact cosine field; full chain RHS -> solve -> GradPhi -> E/J -> divJ.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "sphinxsys.h"

#include <cstring>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

struct ManufacturedCli
{
    OpheliePhiMmsSourceKind source_kind = OpheliePhiMmsSourceKind::DiscreteGrad;
    Real phi_gauge_penalty = 0.0;
};

inline ManufacturedCli parseManufacturedCli(int ac, char *av[])
{
    ManufacturedCli cli;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--mms-source=", 13) == 0)
        {
            cli.source_kind = parseOpheliePhiMmsSourceKind(std::string(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--phi-gauge-penalty=", 20) == 0)
        {
            cli.phi_gauge_penalty = static_cast<Real>(std::atof(av[i] + 20));
        }
    }
    return cli;
}

} // namespace

int main(int ac, char *av[])
{
    const ManufacturedCli cli = parseManufacturedCli(ac, av);

    OphelieParameters params;
    params.frequency_ = 300.0e3;
    params.sigma_glass_ = 16.0;
    params.enable_power_scaling_ = false;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_gauge_penalty_ = cli.phi_gauge_penalty;

    const Real dp = 0.08;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.24, 0.24, 0.24);
    const Real div_j_length = halfsize[0];
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    assignManufacturedASrcFromPhiExact<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, center,
                                                            halfsize, cli.source_kind);
    const OpheliePhiEquationResidualMetrics pre_solve_eq =
        evaluatePhiLhsRhsConsistency<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> level0(glass_body, glass_names, params);
    level0.exec();
    const OphelieDivJMetrics div_j_level0 =
        computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_length);
    syncGlassElectromagneticFieldsToHost(particles, glass_names);
    const Real e_level0_norm = std::sqrt(hostVolWeightedVecdNormSquared(particles, glass_names.e_imag, n));

    const Real phi_solver_rel_res = solvePhiImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    const OpheliePhiEquationResidualMetrics post_solve_eq =
        evaluatePhiLhsRhsConsistency<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        *glass_inner, glass_names);
    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, glass_names, params);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    compute_ejq_with_phi.exec();

    const OphelieDivJMetrics div_j_phi =
        computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_length);
    syncGlassElectromagneticFieldsToHost(particles, glass_names);
    const Real e_phi_norm = std::sqrt(hostVolWeightedVecdNormSquared(particles, glass_names.e_imag, n));
    const Real max_j_imag = hostVecdFieldMax(particles, glass_names.j_imag, n);

    const Real reconstructed_divj_red = div_j_level0.div_j_weighted_l2 / (div_j_phi.div_j_weighted_l2 + TinyReal);
    const char *source_label =
        cli.source_kind == OpheliePhiMmsSourceKind::DiscreteGrad ? "discrete-grad" : "continuous-grad";
    const bool phi_solver_passed = phi_solver_rel_res < 10.0 * params.phi_gmres_tolerance_;
    const bool eq_post_passed = post_solve_eq.eq_res_vol_l2 < 10.0 * params.phi_gmres_tolerance_;
    const bool divj_improved = reconstructed_divj_red > 1.25;
    const bool passed = n > 0 && phi_solver_passed && eq_post_passed && divj_improved;

    std::cout << "test_3d_ophelie_phi_solve_manufactured_A n=" << n << " mms_source=" << source_label
              << " phi_gauge_penalty=" << params.phi_gauge_penalty_ << " pre_eq_res_vol=" << pre_solve_eq.eq_res_vol_l2
              << " phi_solver_rel_res=" << phi_solver_rel_res << " post_eq_res_vol=" << post_solve_eq.eq_res_vol_l2
              << " E_L0=" << e_level0_norm << " E_phi=" << e_phi_norm << " max_JImag=" << max_j_imag
              << " divJ_L2_L0=" << div_j_level0.div_j_weighted_l2 << " divJ_L2_phi=" << div_j_phi.div_j_weighted_l2
              << " reconstructed_divJ_red=" << reconstructed_divj_red << " phi_solver_passed=" << (phi_solver_passed ? 1 : 0)
              << " divJ_improved=" << (divj_improved ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
