/**
 * @file test_3d_ophelie_phi_neumann_cylinder.cpp
 * @brief Cylinder MMS: phi=alpha*z^2/2, A=-grad(phi)/omega => E=0; analytic cylinder Neumann.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline void assignQuadraticCylinderFields(BaseParticles &particles, const OphelieGlassFieldNames &names, Real alpha,
                                          Real omega)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        phi[i] = Real(0.5) * alpha * z * z;
        a_src[i] = Vecd(0.0, 0.0, -alpha * z / omega);
    }
    syncVariableToDevice<Real>(particles, names.phi_imag);
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
}

inline Real measureCylinderPhiExactResidualL2(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                            OphelieParameters &params, const OpheliePhiBoundaryGeometryContext &geom,
                                            Real dp, Real alpha, Real omega)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    assignQuadraticCylinderFields(particles, names, alpha, omega);

    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, &geom, dp, nullptr);
    applyOpheliePhiImagLhsOperator<MainExecutionPolicy>(glass_body, inner, names, params);

    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real res_l2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real d = lhs[i] - rhs[i];
        res_l2 += vol[i] * d * d;
    }
    return std::sqrt(res_l2);
}

} // namespace

int main(int, char *[])
{
    OphelieFrenchReducedCaseParams french;
    OphelieParameters params;
    french.dp = 0.06;
    french.glass_radius = 0.15;
    french.glass_half_height = 0.12;
    applyFrenchReducedDefaults(params, french);
    refreshFrenchReducedCoilStack(french);

    params.enable_power_scaling_ = false;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    params.phi_gauge_penalty_ = 0.0;
    params.phi_boundary_mode_ = OpheliePhiBoundaryMode::OneSidedNeumann;
    params.phi_boundary_grad_neumann_projection_ = true;
    params.phi_boundary_lhs_grad_neumann_ = false;
    params.phi_boundary_normal_source_ = OpheliePhiBoundaryNormalSource::AnalyticCylinder;

    const Real alpha = 1.5;
    const Real dp = french.dp;
    const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
    const Real omega = params.omega();

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, 2.0 * dp);
    SPHSystem sph_system(bounds, dp);
    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height,
                                                                              12));
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

    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
    geom.french = french;
    setupOpheliePhiBoundaryParticleFields(particles, glass_names, params, geom, dp);

    OphelieParameters params_no_neumann = params;
    params_no_neumann.phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
    params_no_neumann.phi_boundary_lhs_grad_neumann_ = false;
    const Real res_no_neumann = measureCylinderPhiExactResidualL2(glass_body, *glass_inner, glass_names,
                                                                  params_no_neumann, geom, dp, alpha, omega);

    const Real res_neumann = measureCylinderPhiExactResidualL2(glass_body, *glass_inner, glass_names, params, geom, dp,
                                                               alpha, omega);

    assignQuadraticCylinderFields(particles, glass_names, alpha, omega);
    StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, glass_names.phi_imag);
    zero_phi.exec();
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    finalizeOpheliePhiImagRhsHost(particles, glass_names, params, &geom, dp, nullptr);
    const Real phi_solver_rel_res =
        solvePhiImagWithCurrentRhs<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq(glass_body, glass_names, params);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    execOphelieScalarPhiGradient<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>,
                                 ComputeOphelieScalarPhiGradientCorrectedCK<Inner<>>,
                                 ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<>>>(
        glass_body, *glass_inner, glass_names, params);
    applyOpheliePhiBoundaryGradNeumannProjectionDynamics<MainExecutionPolicy>(glass_body, glass_names, params, false,
                                                                              &geom, dp);
    compute_ejq.exec();

    const OpheliePhiBoundaryJnMetrics jn_post = computeFrenchCylinderBoundaryJnMetrics(
        particles, glass_names, n, french, params, boundary_width);

    const bool jn_ok = jn_post.jn_boundary_rel < 1.0e-4;
    const bool solver_ok = phi_solver_rel_res < 0.15;
    const bool passed = n > 0 && jn_ok && solver_ok;

    std::cout << "test_3d_ophelie_phi_neumann_cylinder n=" << n << " res_no_neumann=" << res_no_neumann
              << " res_neumann=" << res_neumann << " phi_solver_rel_res=" << phi_solver_rel_res
              << " Jn_post_rel=" << jn_post.jn_boundary_rel << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
