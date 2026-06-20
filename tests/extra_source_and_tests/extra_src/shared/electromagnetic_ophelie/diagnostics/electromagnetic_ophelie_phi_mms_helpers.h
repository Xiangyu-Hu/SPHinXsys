#ifndef ELECTROMAGNETIC_OPHELIE_PHI_MMS_HELPERS_H
#define ELECTROMAGNETIC_OPHELIE_PHI_MMS_HELPERS_H

#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_phi.h"
#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OpheliePhiMmsSourceKind
{
    ContinuousGrad,
    DiscreteGrad
};

inline OpheliePhiMmsSourceKind parseOpheliePhiMmsSourceKind(const std::string &name)
{
    if (name == "discrete-grad" || name == "discrete_grad")
    {
        return OpheliePhiMmsSourceKind::DiscreteGrad;
    }
    return OpheliePhiMmsSourceKind::ContinuousGrad;
}

class OphelieTestGlassBoxShape : public ComplexShape
{
  public:
    OphelieTestGlassBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

inline Real manufacturedPhiExactCosine(const Vecd &pos, const Vecd &center, const Vecd &halfsize)
{
    const Real lx = 2.0 * halfsize[0];
    const Real ly = 2.0 * halfsize[1];
    const Real lz = 2.0 * halfsize[2];
    const Real x = (pos[0] - center[0]) / lx;
    const Real y = (pos[1] - center[1]) / ly;
    const Real z = (pos[2] - center[2]) / lz;
    return std::cos(Pi * x) * std::cos(Pi * y) * std::cos(Pi * z);
}

inline Vecd manufacturedGradPhiExactCosine(const Vecd &pos, const Vecd &center, const Vecd &halfsize)
{
    const Real lx = 2.0 * halfsize[0];
    const Real ly = 2.0 * halfsize[1];
    const Real lz = 2.0 * halfsize[2];
    const Real x = (pos[0] - center[0]) / lx;
    const Real y = (pos[1] - center[1]) / ly;
    const Real z = (pos[2] - center[2]) / lz;
    const Real cx = std::cos(Pi * x);
    const Real cy = std::cos(Pi * y);
    const Real cz = std::cos(Pi * z);
    const Real sx = std::sin(Pi * x);
    const Real sy = std::sin(Pi * y);
    const Real sz = std::sin(Pi * z);
    return Vecd(-Pi / lx * sx * cy * cz, -Pi / ly * cx * sy * cz, -Pi / lz * cx * cy * sz);
}

struct OpheliePhiEquationResidualMetrics
{
    Real eq_res_linf = 0.0;
    Real eq_res_vol_l2 = 0.0;
    Real rhs_vol_l2 = 0.0;
};

inline OpheliePhiEquationResidualMetrics computeHostPhiEquationResidual(BaseParticles &particles,
                                                                        const OphelieGlassFieldNames &names,
                                                                        size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OpheliePhiEquationResidualMetrics metrics;
    Real rhs_max = 0.0;
    Real residual_max = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real residual = lhs[i] - rhs[i];
        residual_max = std::max(residual_max, std::abs(residual));
        rhs_max = std::max(rhs_max, std::abs(rhs[i]));
        metrics.eq_res_vol_l2 += vol[i] * residual * residual;
        metrics.rhs_vol_l2 += vol[i] * rhs[i] * rhs[i];
    }
    metrics.eq_res_linf = residual_max / (rhs_max + TinyReal);
    metrics.eq_res_vol_l2 = std::sqrt(metrics.eq_res_vol_l2) / (std::sqrt(metrics.rhs_vol_l2) + TinyReal);
    return metrics;
}

template <class ExecutionPolicy>
inline void assignManufacturedPhiExactCosine(SolidBody &glass_body, const OphelieGlassFieldNames &names,
                                            const Vecd &center, const Vecd &halfsize)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *phi_imag = particles.getVariableDataByName<Real>(names.phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        phi_imag[i] = manufacturedPhiExactCosine(pos[i], center, halfsize);
    }
    syncVariableToDevice<Real>(particles, names.phi_imag);
}

template <class ExecutionPolicy>
inline void assignManufacturedASrcFromPhiExact(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                               const OphelieParameters &params, const Vecd &center,
                                               const Vecd &halfsize, OpheliePhiMmsSourceKind source_kind)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const Real omega = params.omega();
    assignManufacturedPhiExactCosine<ExecutionPolicy>(glass_body, names, center, halfsize);

    Vecd *a_src_real = particles.getVariableDataByName<Vecd>(names.a_src_real);
    if (source_kind == OpheliePhiMmsSourceKind::ContinuousGrad)
    {
        syncVariableToHost<Vecd>(particles, "Position");
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < n; ++i)
        {
            a_src_real[i] = -manufacturedGradPhiExactCosine(pos[i], center, halfsize) / omega;
        }
    }
    else
    {
        UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
        UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        execOphelieScalarPhiGradient<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>,
                                     ComputeOphelieScalarPhiGradientCorrectedCK<Inner<>>,
                                     ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<>>>(glass_body, inner, names,
                                                                                             params);
        syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
        const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
        for (size_t i = 0; i < n; ++i)
        {
            a_src_real[i] = -grad_phi[i] / omega;
        }
    }
    syncVariableToDevice<Vecd>(particles, names.a_src_real);
}

/** ||L phi - RHS||_vol / ||RHS||_vol from current phi_lhs_imag / phi_rhs_imag on host. */
inline Real hostPhiEqResVolFromCurrentLhsRhs(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                             size_t total_real_particles)
{
    return computeHostPhiEquationResidual(particles, names, total_real_particles).eq_res_vol_l2;
}

struct OpheliePhiEqResidualShellStats
{
    Real eq_res_interior_vol = 0.0;
    Real eq_res_boundary_vol = 0.0;
    Real interior_vol_fraction = 0.0;
};

/** Split ||L phi - RHS||_vol between interior (r_xy < interior_radius) and boundary shell. */
inline OpheliePhiEqResidualShellStats computeHostPhiEqResidualShellStats(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t total_real_particles,
    const Vecd &center, Real interior_radius)
{
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OpheliePhiEqResidualShellStats stats;
    Real interior_rhs_l2 = 0.0;
    Real boundary_rhs_l2 = 0.0;
    Real interior_res_l2 = 0.0;
    Real boundary_res_l2 = 0.0;
    Real total_vol = 0.0;
    Real interior_vol = 0.0;
    const Real r_int_sq = interior_radius * interior_radius;

    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Vecd d = pos[i] - center;
        const Real r_sq = d[0] * d[0] + d[1] * d[1];
        const Real residual = lhs[i] - rhs[i];
        total_vol += vol[i];
        if (r_sq <= r_int_sq)
        {
            interior_vol += vol[i];
            interior_res_l2 += vol[i] * residual * residual;
            interior_rhs_l2 += vol[i] * rhs[i] * rhs[i];
        }
        else
        {
            boundary_res_l2 += vol[i] * residual * residual;
            boundary_rhs_l2 += vol[i] * rhs[i] * rhs[i];
        }
    }
    stats.interior_vol_fraction = interior_vol / (total_vol + TinyReal);
    stats.eq_res_interior_vol = std::sqrt(interior_res_l2) / (std::sqrt(interior_rhs_l2) + TinyReal);
    stats.eq_res_boundary_vol = std::sqrt(boundary_res_l2) / (std::sqrt(boundary_rhs_l2) + TinyReal);
    return stats;
}

template <class ExecutionPolicy>
inline OpheliePhiEquationResidualMetrics evaluatePhiLhsRhsConsistency(SolidBody &glass_body, Inner<> &inner,
                                                                     const OphelieGlassFieldNames &names,
                                                                     const OphelieParameters &params)
{
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    return computeHostPhiEquationResidual(glass_body.getBaseParticles(), names,
                                          glass_body.getBaseParticles().TotalRealParticles());
}

/** Compare div(sigma grad phi) [uncorrected div kernel] vs PairwiseLaplace(phi) [distance-weighted]. */
template <class ExecutionPolicy>
inline OpheliePhiEquationResidualMetrics evaluateDivSigmaGradVsLaplaceConsistency(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(inner, names);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdFieldBySigmaCK> scale_sigma_grad(
        glass_body, names, names.grad_phi_imag, names.j_imag);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCK<Inner<>>> compute_div_sigma_grad(
        inner, names.j_imag, names.div_j_imag);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_laplace(glass_body, names.phi_lhs_imag);
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> apply_laplace(
        inner, names.phi_imag, names.sigma, names.phi_lhs_imag, params.pair_weight_regularization_);

    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    scale_sigma_grad.exec();
    compute_div_sigma_grad.exec();
    zero_laplace.exec();
    apply_laplace.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    syncVariableToHost<Real>(particles, names.div_j_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const Real *div_sigma_grad = particles.getVariableDataByName<Real>(names.div_j_imag);
    const Real *laplace_phi = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OpheliePhiEquationResidualMetrics metrics;
    Real laplace_vol_l2 = 0.0;
    Real residual_max = 0.0;
    Real laplace_max = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real residual = div_sigma_grad[i] - laplace_phi[i];
        residual_max = std::max(residual_max, std::abs(residual));
        laplace_max = std::max(laplace_max, std::abs(laplace_phi[i]));
        metrics.eq_res_vol_l2 += vol[i] * residual * residual;
        laplace_vol_l2 += vol[i] * laplace_phi[i] * laplace_phi[i];
    }
    metrics.eq_res_linf = residual_max / (laplace_max + TinyReal);
    metrics.eq_res_vol_l2 = std::sqrt(metrics.eq_res_vol_l2) / (std::sqrt(laplace_vol_l2) + TinyReal);
    metrics.rhs_vol_l2 = std::sqrt(laplace_vol_l2);
    return metrics;
}

/** With phi fixed and ASrc = -GradPhiDiscrete/omega, EImag = -GradPhi - omega*ASrc should vanish pointwise. */
template <class ExecutionPolicy>
inline Real evaluateDiscreteManufacturedZeroEImagResidual(SolidBody &glass_body, Inner<> &inner,
                                                          const OphelieGlassFieldNames &names,
                                                          const OphelieParameters &params)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(inner, names);
    StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq(glass_body, names, params);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    compute_ejq.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    return std::sqrt(hostVolWeightedVecdNormSquared(particles, names.e_imag, n));
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_MMS_HELPERS_H
