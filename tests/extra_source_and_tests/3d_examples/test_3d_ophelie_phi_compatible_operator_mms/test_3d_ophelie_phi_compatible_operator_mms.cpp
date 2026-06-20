/**
 * @file test_3d_ophelie_phi_compatible_operator_mms.cpp
 * @brief MMS box: uncorrected vs paired G_c/D_c linearity and discrete self-consistency.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_operator_diagnostics.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline bool passLinearity(const OpheliePhiOperatorLinearityMetrics &metrics)
{
    const Real repeat_tol = sizeof(Real) >= 8 ? Real(1.0e-10) : Real(1.0e-6);
    const Real linear_tol = sizeof(Real) >= 8 ? Real(1.0e-8) : Real(1.0e-6);
    return metrics.repeat_apply_rel < repeat_tol && metrics.linearity_add_rel < linear_tol &&
           metrics.linearity_scale_rel < linear_tol;
}

inline bool passSelfConsistency(const OpheliePhiDiscreteSelfConsistencyMetrics &metrics)
{
    return metrics.discrete_self_consistency_rel < Real(1.0e-5);
}

struct ModeResult
{
    const char *label = "";
    OpheliePhiOperatorLinearityMetrics linearity;
    OpheliePhiDiscreteSelfConsistencyMetrics self_consistency;
    OpheliePhiCompatibleOperatorMetrics operator_compare;
    bool linearity_passed = false;
    bool self_consistency_passed = false;
};

inline ModeResult evaluateMode(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                               OphelieParameters params, const Vecd &center, const Vecd &halfsize, const char *label)
{
    ModeResult result;
    result.label = label;
    result.linearity = evaluateOpheliePhiLhsLinearityMetrics<MainExecutionPolicy>(glass_body, inner, names, params);
    result.linearity_passed = passLinearity(result.linearity);

    assignManufacturedPhiExactCosine<MainExecutionPolicy>(glass_body, names, center, halfsize);
    assignManufacturedASrcFromPhiExact<MainExecutionPolicy>(glass_body, inner, names, params, center, halfsize,
                                                            OpheliePhiMmsSourceKind::DiscreteGrad);
    result.self_consistency =
        evaluateOpheliePhiDiscreteSelfConsistency<MainExecutionPolicy>(glass_body, inner, names, params, center,
                                                                         halfsize);
    result.self_consistency_passed = passSelfConsistency(result.self_consistency);
    result.operator_compare =
        evaluateOpheliePhiCompatibleVsUncorrectedOperators<MainExecutionPolicy>(glass_body, inner, names, params);
    return result;
}

} // namespace

int main(int, char *[])
{
    OphelieParameters params;
    params.frequency_ = 300.0e3;
    params.sigma_glass_ = 16.0;
    params.phi_gauge_penalty_ = 0.0;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;

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

    OphelieParameters unc_params = params;
    unc_params.phi_compatible_correction_ = false;
    const ModeResult unc =
        evaluateMode(glass_body, *glass_inner, glass_names, unc_params, center, halfsize, "uncorrected");

    OphelieParameters comp_params = params;
    comp_params.phi_compatible_correction_ = true;
    const ModeResult comp =
        evaluateMode(glass_body, *glass_inner, glass_names, comp_params, center, halfsize, "compatible");

    const bool passed = n > 0 && unc.linearity_passed && unc.self_consistency_passed && comp.linearity_passed &&
                        comp.self_consistency_passed;

    auto print_mode = [](const ModeResult &mode)
    {
        std::cout << " mode=" << mode.label << " linearity_passed=" << (mode.linearity_passed ? 1 : 0)
                  << " repeat_apply_rel=" << mode.linearity.repeat_apply_rel
                  << " linearity_add_rel=" << mode.linearity.linearity_add_rel
                  << " linearity_scale_rel=" << mode.linearity.linearity_scale_rel
                  << " self_consistency_passed=" << (mode.self_consistency_passed ? 1 : 0)
                  << " discrete_self_consistency_rel=" << mode.self_consistency.discrete_self_consistency_rel
                  << " rhs_unc_vs_compatible_vol=" << mode.operator_compare.rhs_unc_vs_compatible_vol
                  << " eq_res_vol_compatible=" << mode.operator_compare.eq_res_vol_compatible << std::endl;
    };

    std::cout << "test_3d_ophelie_phi_compatible_operator_mms n=" << n << " passed=" << (passed ? 1 : 0) << std::endl;
    print_mode(unc);
    print_mode(comp);
    return passed ? 0 : 1;
}
