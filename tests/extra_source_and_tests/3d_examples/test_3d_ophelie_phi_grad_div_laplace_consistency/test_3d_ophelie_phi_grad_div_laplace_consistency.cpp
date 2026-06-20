/**
 * @file test_3d_ophelie_phi_grad_div_laplace_consistency.cpp
 * @brief Compare div(sigma grad phi) [uncorrected] vs PairwiseLaplace(phi) [distance-weighted].
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
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
    assignManufacturedPhiExactCosine<MainExecutionPolicy>(glass_body, glass_names, center, halfsize);

    const OpheliePhiEquationResidualMetrics grad_div_vs_laplace =
        evaluateDivSigmaGradVsLaplaceConsistency<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    assignManufacturedASrcFromPhiExact<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, center,
                                                            halfsize, OpheliePhiMmsSourceKind::DiscreteGrad);
    const Real e_imag_norm =
        evaluateDiscreteManufacturedZeroEImagResidual<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    const bool passed = n > 0 && grad_div_vs_laplace.eq_res_vol_l2 < 0.15 && e_imag_norm < 1.0e-10;

    std::cout << "test_3d_ophelie_phi_grad_div_laplace_consistency n=" << n
              << " graddiv_vs_laplace_vol=" << grad_div_vs_laplace.eq_res_vol_l2
              << " graddiv_vs_laplace_linf=" << grad_div_vs_laplace.eq_res_linf
              << " discrete_mms_EImag_norm=" << e_imag_norm << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
