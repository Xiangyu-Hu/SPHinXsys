/**
 * @file test_3d_ophelie_phi_operator_linear_consistency.cpp
 * @brief Verify DivSigmaGrad LHS repeat-apply, additivity, and scaling on random phi.
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
    params.phi_gauge_penalty_ = 0.0;
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

    const OpheliePhiOperatorLinearityMetrics metrics =
        evaluateOpheliePhiLhsLinearityMetrics<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    const Real repeat_tol = sizeof(Real) >= 8 ? Real(1.0e-10) : Real(1.0e-6);
    const Real linear_tol = sizeof(Real) >= 8 ? Real(1.0e-8) : Real(1.0e-6);
    const bool passed = n > 0 && metrics.repeat_apply_rel < repeat_tol && metrics.linearity_add_rel < linear_tol &&
                        metrics.linearity_scale_rel < linear_tol;

    std::cout << "test_3d_ophelie_phi_operator_linear_consistency n=" << n
              << " repeat_apply_rel=" << metrics.repeat_apply_rel
              << " linearity_add_rel=" << metrics.linearity_add_rel
              << " linearity_scale_rel=" << metrics.linearity_scale_rel << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
