/**
 * @file test_3d_ophelie_phi_discrete_self_consistency.cpp
 * @brief Discrete MMS: phi_exact -> A=-G_h(phi)/omega -> compare L_h(phi) vs b_h.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_operator_diagnostics.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

int main(int, char *[])
{
    OphelieParameters params;
    params.frequency_ = 300.0e3;
    params.sigma_glass_ = 16.0;
    params.enable_power_scaling_ = false;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    params.phi_gauge_penalty_ = 0.0;

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

    const OpheliePhiDiscreteSelfConsistencyMetrics metrics =
        evaluateOpheliePhiDiscreteSelfConsistency<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params,
                                                                        center, halfsize);

    const bool passed = n > 0 && metrics.discrete_self_consistency_rel < Real(1.0e-6);

    std::cout << "test_3d_ophelie_phi_discrete_self_consistency n=" << n
              << " lhs_minus_rhs_vol_l2=" << metrics.lhs_minus_rhs_vol_l2
              << " rhs_vol_l2=" << metrics.rhs_vol_l2
              << " discrete_self_consistency_rel=" << metrics.discrete_self_consistency_rel
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
