/**
 * @file test_3d_ophelie_phi_laplace_rhs_consistency.cpp
 * @brief Compare PairwiseLaplace(PhiExact) + gauge vs PhiRhsFromASrc before solving.
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

struct ConsistencyCli
{
    OpheliePhiMmsSourceKind source_kind = OpheliePhiMmsSourceKind::DiscreteGrad;
    Real phi_gauge_penalty = 0.0;
};

inline ConsistencyCli parseConsistencyCli(int ac, char *av[])
{
    ConsistencyCli cli;
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
    const ConsistencyCli cli = parseConsistencyCli(ac, av);

    OphelieParameters params;
    params.frequency_ = 300.0e3;
    params.sigma_glass_ = 16.0;
    params.phi_gauge_penalty_ = cli.phi_gauge_penalty;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;

    const Real dp = 0.08;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.24, 0.24, 0.24);
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
    const size_t n = glass_body.getBaseParticles().TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    assignManufacturedASrcFromPhiExact<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, center,
                                                            halfsize, cli.source_kind);
    const OpheliePhiEquationResidualMetrics eq_metrics =
        evaluatePhiLhsRhsConsistency<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    const char *source_label =
        cli.source_kind == OpheliePhiMmsSourceKind::DiscreteGrad ? "discrete-grad" : "continuous-grad";
    const bool passed = n > 0 && eq_metrics.eq_res_vol_l2 < 0.35;

    std::cout << "test_3d_ophelie_phi_laplace_rhs_consistency n=" << n << " mms_source=" << source_label
              << " phi_gauge_penalty=" << params.phi_gauge_penalty_ << " eq_res_linf=" << eq_metrics.eq_res_linf
              << " eq_res_vol=" << eq_metrics.eq_res_vol_l2 << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
