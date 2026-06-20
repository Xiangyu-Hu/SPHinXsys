/**
 * @file test_3d_ophelie_phi_neumann_slab.cpp
 * @brief Slab MMS: phi=alpha*x^2/2, A=-grad(phi)/omega => E=0; Neumann improves L(phi)-b at phi_exact.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

class OphelieNeumannSlabBoxShape : public ComplexShape
{
  public:
    OphelieNeumannSlabBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

inline void assignQuadraticSlabFields(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                      const Vecd &center, Real alpha, Real omega)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    Vecd *a_src = particles.getVariableDataByName<Vecd>(names.a_src_real);
    for (size_t i = 0; i < n; ++i)
    {
        const Real x = pos[i][0] - center[0];
        phi[i] = Real(0.5) * alpha * x * x;
        a_src[i] = Vecd(-alpha * x / omega, 0.0, 0.0);
    }
    syncVariableToDevice<Real>(particles, names.phi_imag);
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
}

inline Real measurePhiExactResidualL2(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                      OphelieParameters &params, const OpheliePhiBoundaryGeometryContext &geom,
                                      Real dp, const Vecd &center, Real alpha, Real omega)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    assignQuadraticSlabFields(particles, names, center, alpha, omega);

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

int main(int ac, char *av[])
{
    (void)ac;
    (void)av;

    OphelieParameters params;
    params.frequency_ = 100.0;
    params.sigma_glass_ = 10.0;
    params.enable_power_scaling_ = false;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
    params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::DivSigmaGrad;
    params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    params.phi_gauge_penalty_ = 0.0;
    params.phi_boundary_normal_source_ = OpheliePhiBoundaryNormalSource::AnalyticBox;
    params.phi_boundary_mode_ = OpheliePhiBoundaryMode::OneSidedNeumann;
    params.phi_boundary_grad_neumann_projection_ = true;
    params.phi_boundary_lhs_grad_neumann_ = false;

    const Real alpha = 2.0;
    const Real dp = 0.06;
    const Vecd center(0.0, 0.0, 0.25);
    const Vecd halfsize(0.18, 0.18, 0.18);
    const Real boundary_width = params.phi_boundary_distance_factor_ * dp;
    const Real omega = params.omega();
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieNeumannSlabBoxShape>("GlassBody", center, halfsize));
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
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = halfsize;
    setupOpheliePhiBoundaryParticleFields(particles, glass_names, params, geom, dp);

    OphelieParameters params_no_neumann = params;
    params_no_neumann.phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
    params_no_neumann.phi_boundary_lhs_grad_neumann_ = false;
    const Real res_no_neumann =
        measurePhiExactResidualL2(glass_body, *glass_inner, glass_names, params_no_neumann, geom, dp, center, alpha,
                                  omega);

    const Real res_neumann =
        measurePhiExactResidualL2(glass_body, *glass_inner, glass_names, params, geom, dp, center, alpha, omega);

    assignQuadraticSlabFields(particles, glass_names, center, alpha, omega);
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

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, glass_names.grad_phi_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(glass_names.grad_phi_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real grad_err_l2 = 0.0;
    Real grad_ref_l2 = 0.0;
    const Real interior_margin = 2.0 * dp;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd r = pos[i] - center;
        if (std::abs(r[0]) > halfsize[0] - interior_margin || std::abs(r[1]) > halfsize[1] - interior_margin ||
            std::abs(r[2]) > halfsize[2] - interior_margin)
        {
            continue;
        }
        const Real x = r[0];
        const Real err_x = grad_phi[i][0] - alpha * x;
        grad_err_l2 += vol[i] * err_x * err_x;
        grad_ref_l2 += vol[i] * alpha * alpha * x * x;
    }
    const Real grad_phi_rel_err = std::sqrt(grad_err_l2) / (std::sqrt(grad_ref_l2) + TinyReal);

    const OpheliePhiBoundaryJnMetrics jn_post =
        computeBoxBoundaryJnMetrics(particles, glass_names, n, params, center, halfsize, boundary_width);

    const bool jn_ok = jn_post.jn_boundary_rel < 1.0e-4;
    const bool solver_ok = phi_solver_rel_res < 0.15;
    const bool passed = n > 0 && jn_ok && solver_ok;

    std::cout << "test_3d_ophelie_phi_neumann_slab n=" << n << " res_no_neumann=" << res_no_neumann
              << " res_neumann=" << res_neumann
              << " phi_solver_rel_res=" << phi_solver_rel_res << " grad_phi_rel_err=" << grad_phi_rel_err
              << " Jn_post_rel=" << jn_post.jn_boundary_rel << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
