#ifndef ELECTROMAGNETIC_OPHELIE_PHI_COMPONENT_H
#define ELECTROMAGNETIC_OPHELIE_PHI_COMPONENT_H

#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_phi_gmres.h"
#include "electromagnetic_ophelie_phi_boundary_diagnostics.h"
#include "electromagnetic_ophelie_phi_rhs_diagnostics.h"
#include "electromagnetic_ophelie_progress.h"

#include <Eigen/Dense>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Generic gauge penalty: lhs += penalty * phi for any scalar phi/lhs pair. */
class OpheliePhiGaugePenaltyToLhsCK : public LocalDynamics
{
  public:
    OpheliePhiGaugePenaltyToLhsCK(SPHBody &sph_body, const std::string &phi_field, const std::string &phi_lhs_field,
                                  Real phi_gauge_penalty)
        : LocalDynamics(sph_body), phi_gauge_penalty_(phi_gauge_penalty),
          dv_phi_(particles_->template getVariableByName<Real>(phi_field)),
          dv_phi_lhs_(particles_->template getVariableByName<Real>(phi_lhs_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : phi_gauge_penalty_(encloser.phi_gauge_penalty_), phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
              phi_lhs_(encloser.dv_phi_lhs_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            phi_lhs_[index_i] += phi_gauge_penalty_ * phi_[index_i];
        }

      protected:
        Real phi_gauge_penalty_;
        Real *phi_;
        Real *phi_lhs_;
    };

  protected:
    Real phi_gauge_penalty_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Real> *dv_phi_lhs_;
};

inline Real hostPhiComponentEqResVolFromCurrentLhsRhs(BaseParticles &particles, const std::string &phi_rhs_field,
                                                      const std::string &phi_lhs_field, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, phi_rhs_field);
    syncVariableToHost<Real>(particles, phi_lhs_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(phi_rhs_field);
    const Real *lhs = particles.getVariableDataByName<Real>(phi_lhs_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real residual_sq = 0.0;
    Real rhs_sq = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real residual = lhs[i] - rhs[i];
        residual_sq += vol[i] * residual * residual;
        rhs_sq += vol[i] * rhs[i] * rhs[i];
    }
    return std::sqrt(residual_sq / (rhs_sq + TinyReal));
}

template <class ExecutionPolicy>
inline void applyOpheliePhiComponentLegacyPairwiseLhsOperator(SolidBody &glass_body, Inner<> &inner,
                                                                const OphelieGlassFieldNames &names,
                                                                const OphelieEdgeFluxComponent &component,
                                                                const OphelieParameters &params)
{
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> apply_laplace(
        inner, component.phi_field, names.sigma, component.phi_lhs_field, params.pair_weight_regularization_);
    StateDynamics<ExecutionPolicy, OpheliePhiGaugePenaltyToLhsCK> apply_gauge_penalty(
        glass_body, component.phi_field, component.phi_lhs_field, params.phi_gauge_penalty_);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_lhs(glass_body, component.phi_lhs_field);
    zero_lhs.exec();
    apply_laplace.exec();
    apply_gauge_penalty.exec();
}

template <class ExecutionPolicy>
inline void applyOpheliePhiComponentLhsOperator(SolidBody &glass_body, Inner<> &inner,
                                                const OphelieGlassFieldNames &names,
                                                const OphelieEdgeFluxComponent &component,
                                                const OphelieParameters &params)
{
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        std::cout << "[ophelie] WARNING: DivSigmaGrad LHS not yet parameterized for edge-flux component "
                  << component.chain_label << "; using legacy pairwise." << std::endl;
    }
    applyOpheliePhiComponentLegacyPairwiseLhsOperator<ExecutionPolicy>(glass_body, inner, names, component, params);
}

template <class ExecutionPolicy>
inline Real evaluateOpheliePhiComponentRelativeResidual(SolidBody &glass_body, Inner<> &inner,
                                                        const OphelieGlassFieldNames &names,
                                                        const OphelieEdgeFluxComponent &component,
                                                        const OphelieParameters &params)
{
    applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, inner, names, component, params);
    BaseParticles &particles = glass_body.getBaseParticles();
    syncVariableToHost<Real>(particles, component.phi_rhs_field);
    syncVariableToHost<Real>(particles, component.phi_lhs_field);
    const size_t n = particles.TotalRealParticles();
    const Real *rhs = particles.getVariableDataByName<Real>(component.phi_rhs_field);
    const Real *lhs = particles.getVariableDataByName<Real>(component.phi_lhs_field);
    Real rhs_max = 0.0;
    Real residual_max = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        rhs_max = std::max(rhs_max, std::abs(rhs[i]));
        residual_max = std::max(residual_max, std::abs(lhs[i] - rhs[i]));
    }
    return residual_max / (rhs_max + TinyReal);
}

template <class ExecutionPolicy>
inline void setupOphelieEdgeFluxComponentRhsFromA(SolidBody &glass_body, Inner<> &inner,
                                                  const OphelieGlassFieldNames &names,
                                                  const OphelieEdgeFluxComponent &component,
                                                  const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_rhs(glass_body, component.phi_rhs_field);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    zero_rhs.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    const Real boundary_width_m =
        params.edge_recon_boundary_width_factor_ *
        glass_body.getBaseParticles().getSPHAdaptation().ReferenceSmoothingLength();
    const OphelieEdgeFluxBoundaryClosureFlags boundary_closure =
        resolveOphelieEdgeFluxBoundaryClosure(params, boundary_width_m, glass_body.getBaseParticles());
    logOphelieEdgeFluxBoundaryClosure(boundary_closure, component.chain_label);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Inner<>>> compute_rhs(
        inner, names, component, params.omega(), params.pair_weight_regularization_, boundary_closure);
    compute_rhs.exec();
}

template <class ExecutionPolicy>
inline void finalizeOpheliePhiComponentRhsHost(BaseParticles &particles, const std::string &phi_rhs_field,
                                               OphelieParameters &params)
{
    if (params.phi_rhs_project_zero_mean_)
    {
        applyOpheliePhiRhsZeroMeanProjection(particles, phi_rhs_field, particles.TotalRealParticles());
    }
}

template <class ExecutionPolicy>
inline Real solvePhiComponentPCGWithCurrentRhs(SolidBody &glass_body, Inner<> &inner,
                                               const OphelieGlassFieldNames &names,
                                               const OphelieEdgeFluxComponent &component,
                                               const OphelieParameters &params)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, component.phi_field);
    zero_phi.exec();
    const std::string pcg_entry_stage = std::string(component.chain_label) + "_pcg_entry";
    logOpheliePhiRhsFingerprint(pcg_entry_stage.c_str(),
                                computeOpheliePhiRhsFingerprint(particles, component.phi_rhs_field, n));

    computeOpheliePhiOperatorDiagonal<ExecutionPolicy>(glass_body, inner, names, params);

    StdVec<Real> solution(n, Real(0));
    StdVec<Real> rhs(n, Real(0));
    StdVec<Real> diagonal(n, Real(0));
    StdVec<Real> residual(n, Real(0));
    StdVec<Real> preconditioned_residual(n, Real(0));
    StdVec<Real> search_direction(n, Real(0));
    StdVec<Real> operator_on_search(n, Real(0));

    hostReadScalarField(particles, component.phi_rhs_field, rhs.data(), n);
    hostReadScalarField(particles, names.phi_laplace_diag, diagonal.data(), n);

    const Real solver_local_rhs_scale = Real(1);

    auto apply_operator = [&](const Real *input, Real *output)
    {
        hostAssignScalarField(particles, component.phi_field, input, n);
        applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, inner, names, component, params);
        hostReadScalarField(particles, component.phi_lhs_field, output, n);
    };

    auto apply_preconditioner = [&](const Real *input, Real *output)
    {
        for (size_t i = 0; i < n; ++i)
        {
            output[i] = input[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
        }
    };

    for (size_t i = 0; i < n; ++i)
    {
        solution[i] = rhs[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
    }

    apply_operator(solution.data(), operator_on_search.data());
    for (size_t i = 0; i < n; ++i)
    {
        residual[i] = rhs[i] - operator_on_search[i];
    }

    const Real rhs_norm = hostVolWeightedNorm(particles, rhs.data(), n);
    const Real rhs_max = hostScalarFieldMax(particles, component.phi_rhs_field, n);
    ophelieUnscalePhiPcgSolution(solution, solver_local_rhs_scale);
    hostAssignScalarField(particles, component.phi_field, solution.data(), n);
    Real relative_residual_linf =
        evaluateOpheliePhiComponentRelativeResidual<ExecutionPolicy>(glass_body, inner, names, component, params);
    if (rhs_max > TinyReal && relative_residual_linf < params.phi_pcg_tolerance_)
    {
        return relative_residual_linf;
    }
    if (solver_local_rhs_scale < Real(1) - TinyReal)
    {
        ophelieScalePhiPcgWorkspaceRhs(solution, solver_local_rhs_scale);
    }

    apply_preconditioner(residual.data(), preconditioned_residual.data());
    search_direction = preconditioned_residual;
    Real rz_old = hostVolWeightedDot(particles, residual.data(), preconditioned_residual.data(), n);

    OphelieProgressLogger pcg_progress(std::string("phi_pcg_") + component.chain_label);
    pcg_progress.log("n=" + std::to_string(n) + " max_iter=" + std::to_string(params.phi_pcg_max_iterations_));

    for (size_t iter = 0; iter < params.phi_pcg_max_iterations_; ++iter)
    {
        apply_operator(search_direction.data(), operator_on_search.data());
        const Real search_energy =
            hostVolWeightedDot(particles, search_direction.data(), operator_on_search.data(), n);
        if (std::abs(search_energy) < TinyReal)
        {
            break;
        }

        const Real step_length = rz_old / (search_energy + TinyReal);
        hostScaledAdd(solution.data(), step_length, search_direction.data(), n);
        for (size_t i = 0; i < n; ++i)
        {
            residual[i] -= step_length * operator_on_search[i];
        }

        if ((iter + 1) % 5 == 0 || iter + 1 == params.phi_pcg_max_iterations_ || (iter + 1) % 100 == 0)
        {
            StdVec<Real> solution_physical = solution;
            ophelieUnscalePhiPcgSolution(solution_physical, solver_local_rhs_scale);
            hostAssignScalarField(particles, component.phi_field, solution_physical.data(), n);
            relative_residual_linf =
                evaluateOpheliePhiComponentRelativeResidual<ExecutionPolicy>(glass_body, inner, names, component, params);
            if ((iter + 1) % 100 == 0 || iter + 1 == params.phi_pcg_max_iterations_)
            {
                pcg_progress.log("iter " + std::to_string(iter + 1) + " rel_res_linf=" +
                                 std::to_string(relative_residual_linf));
            }
            if (relative_residual_linf < params.phi_pcg_tolerance_)
            {
                solution = solution_physical;
                break;
            }
        }

        apply_preconditioner(residual.data(), preconditioned_residual.data());
        const Real rz_new = hostVolWeightedDot(particles, residual.data(), preconditioned_residual.data(), n);
        const Real beta = rz_new / (rz_old + TinyReal);
        for (size_t i = 0; i < n; ++i)
        {
            search_direction[i] = preconditioned_residual[i] + beta * search_direction[i];
        }
        rz_old = rz_new;
    }

    ophelieUnscalePhiPcgSolution(solution, solver_local_rhs_scale);
    hostAssignScalarField(particles, component.phi_field, solution.data(), n);
    return evaluateOpheliePhiComponentRelativeResidual<ExecutionPolicy>(glass_body, inner, names, component, params);
}

template <class ExecutionPolicy>
inline Real solvePhiComponentGMRESWithCurrentRhs(SolidBody &glass_body, Inner<> &inner,
                                                 const OphelieGlassFieldNames &names,
                                                 const OphelieEdgeFluxComponent &component,
                                                 const OphelieParameters &params)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, component.phi_field);
    zero_phi.exec();
    const std::string gmres_entry_stage = std::string(component.chain_label) + "_gmres_entry";
    logOpheliePhiRhsFingerprint(gmres_entry_stage.c_str(),
                                computeOpheliePhiRhsFingerprint(particles, component.phi_rhs_field, n));

    computeOpheliePhiGMRESPreconditionerDiagonal<ExecutionPolicy>(glass_body, inner, names, params);
    const UnsignedInt restart_dimension = std::max(params.phi_gmres_restart_dimension_, UnsignedInt(1));
    const UnsignedInt max_outer_iterations = std::max(params.phi_gmres_max_outer_iterations_, UnsignedInt(1));

    OphelieProgressLogger gmres_progress(std::string("phi_gmres_") + component.chain_label);
    gmres_progress.log("n=" + std::to_string(n) + " restart=" + std::to_string(restart_dimension) +
                       " max_outer=" + std::to_string(max_outer_iterations));

    StdVec<Real> solution(n, Real(0));
    StdVec<Real> rhs(n, Real(0));
    StdVec<Real> diagonal(n, Real(0));
    StdVec<Real> operator_output(n, Real(0));
    StdVec<Real> residual(n, Real(0));
    StdVec<Real> workspace(n, Real(0));
    StdVec<Real> preconditioned(n, Real(0));
    StdVec<StdVec<Real>> krylov_basis(restart_dimension + 1, StdVec<Real>(n, Real(0)));

    hostReadScalarField(particles, component.phi_rhs_field, rhs.data(), n);
    hostReadScalarField(particles, names.phi_laplace_diag, diagonal.data(), n);

    for (size_t i = 0; i < n; ++i)
    {
        const Real preconditioner = diagonal[i] + params.phi_gauge_penalty_ + TinyReal;
        solution[i] = params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad ? Real(0)
                                                                                               : rhs[i] / preconditioner;
    }
    hostAssignScalarField(particles, component.phi_field, solution.data(), n);

    auto apply_operator = [&](const Real *input, Real *output)
    {
        hostAssignScalarField(particles, component.phi_field, input, n);
        applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, inner, names, component, params);
        hostReadScalarField(particles, component.phi_lhs_field, output, n);
    };

    auto apply_preconditioner = [&](const Real *input, Real *output)
    {
        for (size_t i = 0; i < n; ++i)
        {
            output[i] = input[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
        }
    };

    apply_operator(solution.data(), operator_output.data());
    for (size_t i = 0; i < n; ++i)
    {
        residual[i] = rhs[i] - operator_output[i];
    }

    const Real rhs_norm = hostVolWeightedNorm(particles, rhs.data(), n);
    const Real rhs_max = hostScalarFieldMax(particles, component.phi_rhs_field, n);
    Real beta = hostVolWeightedNorm(particles, residual.data(), n);
    Real relative_residual = beta / (rhs_norm + TinyReal);

    if (!std::isfinite(relative_residual) ||
        (rhs_max > TinyReal && relative_residual < params.phi_gmres_tolerance_))
    {
        hostAssignScalarField(particles, component.phi_field, solution.data(), n);
        return relative_residual;
    }

    StdVec<Real> hessenberg_cos(restart_dimension, Real(0));
    StdVec<Real> hessenberg_sin(restart_dimension, Real(0));

    for (UnsignedInt outer = 0; outer < max_outer_iterations; ++outer)
    {
        const Real inv_beta = Real(1.0) / (beta + TinyReal);
        for (size_t i = 0; i < n; ++i)
        {
            krylov_basis[0][i] = residual[i] * inv_beta;
        }

        Eigen::VectorXd givens_rhs = Eigen::VectorXd::Zero(restart_dimension + 1);
        givens_rhs(0) = beta;
        Eigen::MatrixXd hessenberg = Eigen::MatrixXd::Zero(restart_dimension + 1, restart_dimension);

        UnsignedInt krylov_dimension = 0;
        for (UnsignedInt inner_step = 0; inner_step < restart_dimension; ++inner_step)
        {
            apply_preconditioner(krylov_basis[inner_step].data(), preconditioned.data());
            apply_operator(preconditioned.data(), workspace.data());

            for (UnsignedInt i = 0; i <= inner_step; ++i)
            {
                const Real projection = hostVolWeightedDot(particles, krylov_basis[i].data(), workspace.data(), n);
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step)) = projection;
                hostSubtractScaledVector(workspace.data(), krylov_basis[i].data(), projection, n);
            }

            const Real subdiagonal_norm = hostVolWeightedNorm(particles, workspace.data(), n);
            hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step)) =
                subdiagonal_norm;

            if (!std::isfinite(subdiagonal_norm) || subdiagonal_norm < TinyReal)
            {
                krylov_dimension = inner_step + (subdiagonal_norm < TinyReal ? 1 : 0);
                break;
            }

            const Real inv_subdiagonal_norm = Real(1.0) / (subdiagonal_norm + TinyReal);
            for (size_t k = 0; k < n; ++k)
            {
                krylov_basis[inner_step + 1][k] = workspace[k] * inv_subdiagonal_norm;
            }

            for (UnsignedInt i = 0; i < inner_step; ++i)
            {
                const Real h0 = hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step));
                const Real h1 = hessenberg(static_cast<Eigen::Index>(i + 1), static_cast<Eigen::Index>(inner_step));
                const Real temp = hessenberg_cos[i] * h0 + hessenberg_sin[i] * h1;
                hessenberg(static_cast<Eigen::Index>(i + 1), static_cast<Eigen::Index>(inner_step)) =
                    -hessenberg_sin[i] * h0 + hessenberg_cos[i] * h1;
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step)) = temp;
            }

            const Real h0 = hessenberg(static_cast<Eigen::Index>(inner_step), static_cast<Eigen::Index>(inner_step));
            const Real h1 = hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step));
            const Real denom = std::sqrt(h0 * h0 + h1 * h1 + TinyReal);
            hessenberg_cos[inner_step] = h0 / denom;
            hessenberg_sin[inner_step] = h1 / denom;
            hessenberg(static_cast<Eigen::Index>(inner_step), static_cast<Eigen::Index>(inner_step)) = denom;
            hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step)) = Real(0);

            const Real temp_rhs = hessenberg_cos[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step)) +
                                  hessenberg_sin[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step + 1));
            givens_rhs(static_cast<Eigen::Index>(inner_step + 1)) =
                -hessenberg_sin[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step)) +
                hessenberg_cos[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step + 1));
            givens_rhs(static_cast<Eigen::Index>(inner_step)) = temp_rhs;

            krylov_dimension = inner_step + 1;
            relative_residual = std::abs(givens_rhs(static_cast<Eigen::Index>(inner_step + 1))) / (rhs_norm + TinyReal);
            if (relative_residual < params.phi_gmres_tolerance_)
            {
                break;
            }
        }

        if (krylov_dimension == 0)
        {
            break;
        }

        const Eigen::Index system_size = static_cast<Eigen::Index>(krylov_dimension);
        Eigen::MatrixXd hessenberg_system = hessenberg.topLeftCorner(system_size, system_size);
        Eigen::VectorXd rhs_vector = givens_rhs.head(system_size);
        Eigen::VectorXd coefficients = hessenberg_system.colPivHouseholderQr().solve(rhs_vector);

        for (Eigen::Index column = 0; column < system_size; ++column)
        {
            apply_preconditioner(krylov_basis[static_cast<size_t>(column)].data(), preconditioned.data());
            hostScaledAdd(solution.data(), coefficients(column), preconditioned.data(), n);
        }

        apply_operator(solution.data(), operator_output.data());
        for (size_t i = 0; i < n; ++i)
        {
            residual[i] = rhs[i] - operator_output[i];
        }
        beta = hostVolWeightedNorm(particles, residual.data(), n);
        relative_residual = beta / (rhs_norm + TinyReal);

        if (relative_residual < params.phi_gmres_tolerance_)
        {
            break;
        }
    }

    hostAssignScalarField(particles, component.phi_field, solution.data(), n);
    applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, inner, names, component, params);
    const Real rel_res_linf =
        evaluateOpheliePhiComponentRelativeResidual<ExecutionPolicy>(glass_body, inner, names, component, params);
    gmres_progress.finish("rel_res_linf=" + std::to_string(rel_res_linf));
    return rel_res_linf;
}

template <class ExecutionPolicy>
inline Real solvePhiComponentWithCurrentRhs(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                            const OphelieEdgeFluxComponent &component, const OphelieParameters &params)
{
    if (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES)
    {
        if (std::string(component.chain_label) == "imag")
        {
            return solvePhiImagGMRESWithCurrentRhs<ExecutionPolicy>(glass_body, inner, names, params);
        }
        return solvePhiComponentGMRESWithCurrentRhs<ExecutionPolicy>(glass_body, inner, names, component, params);
    }
    return solvePhiComponentPCGWithCurrentRhs<ExecutionPolicy>(glass_body, inner, names, component, params);
}

struct OphelieComplexEdgeFluxSolveReport
{
    Real phi_imag_solver_rel_residual = 0.0;
    Real phi_real_solver_rel_residual = 0.0;
    Real phi_eq_res_vol_imag = 0.0;
    Real phi_eq_res_vol_real = 0.0;
};

/** Solve imag (+ optional real) edge-flux scalar systems from current active A fields. */
template <class ExecutionPolicy>
inline OphelieComplexEdgeFluxSolveReport solveOphelieComplexEdgeFluxWithCurrentA(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OpheliePhiBoundaryGeometryContext *boundary_geom = nullptr, Real dp = 0.0,
    bool skip_real_phi_solve = false)
{
    OphelieComplexEdgeFluxSolveReport report;
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi_imag(glass_body, names.phi_imag);
    zero_phi_imag.exec();
    if (params.edge_flux_complex_)
    {
        StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi_real(glass_body, names.phi_real);
        zero_phi_real.exec();
    }

    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    finalizeOpheliePhiImagRhsHost(particles, names, params, boundary_geom, dp, nullptr);
    const bool solver_local_active =
        params.edge_flux_normalization_mode_ == OphelieEdgeFluxNormalizationMode::SolverLocal;
    if (solver_local_active)
    {
        params.edge_flux_solver_local_rhs_scale_ = computeOpheliePhiSolverLocalRhsScale(
            particles, names.phi_rhs_imag, n, params.edge_flux_safe_rhs_l2_, params.edge_flux_safe_rhs_max_abs_);
    }
    const Real solver_local_scale = ophelieEdgeFluxEffectiveSolverLocalRhsScale(params);
    if (solver_local_scale < Real(1) - TinyReal)
    {
        std::cout << "[ophelie] solver-local scale phi_rhs_imag on particles by " << solver_local_scale << std::endl;
        hostScaleScalarFieldInPlace(particles, names.phi_rhs_imag, solver_local_scale, n);
    }
    report.phi_imag_solver_rel_residual =
        solvePhiImagWithCurrentRhs<ExecutionPolicy>(glass_body, inner, names, params);
    if (solver_local_scale < Real(1) - TinyReal)
    {
        const Real inv_scale = Real(1) / (solver_local_scale + TinyReal);
        hostScaleScalarFieldInPlace(particles, names.phi_imag, inv_scale, n);
        hostScaleScalarFieldInPlace(particles, names.phi_rhs_imag, inv_scale, n);
    }

    if (ophelieUseEdgeFluxElectromotiveRhs(params) && params.edge_flux_complex_ && !skip_real_phi_solve)
    {
        const OphelieEdgeFluxComponent real_component = makeOphelieEdgeFluxRealComponent(names, params);
        setupOphelieEdgeFluxComponentRhsFromA<ExecutionPolicy>(glass_body, inner, names, real_component, params);
        finalizeOpheliePhiComponentRhsHost<ExecutionPolicy>(particles, real_component.phi_rhs_field, params);
        if (solver_local_scale < Real(1) - TinyReal)
        {
            std::cout << "[ophelie] solver-local scale phi_rhs_real on particles by " << solver_local_scale
                      << std::endl;
            hostScaleScalarFieldInPlace(particles, real_component.phi_rhs_field, solver_local_scale, n);
        }
        report.phi_real_solver_rel_residual = solvePhiComponentWithCurrentRhs<ExecutionPolicy>(
            glass_body, inner, names, real_component, params);
        if (solver_local_scale < Real(1) - TinyReal)
        {
            const Real inv_scale = Real(1) / (solver_local_scale + TinyReal);
            hostScaleScalarFieldInPlace(particles, real_component.phi_field, inv_scale, n);
            hostScaleScalarFieldInPlace(particles, real_component.phi_rhs_field, inv_scale, n);
        }
        applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, inner, names, real_component, params);
        report.phi_eq_res_vol_real = hostPhiComponentEqResVolFromCurrentLhsRhs(
            particles, real_component.phi_rhs_field, real_component.phi_lhs_field, n);
    }

    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    report.phi_eq_res_vol_imag = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);
    return report;
}

/** Solve complex edge-flux, reconstruct E/J/Q, sync primary fields; return P_recon (complex when enabled). */
template <class ExecutionPolicy>
inline Real execOphelieComplexEdgeFluxSolveReconAndPower(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const OpheliePhiBoundaryGeometryContext *boundary_geom = nullptr, Real dp = 0.0,
    OphelieComplexEdgeFluxSolveReport *solve_report = nullptr, bool skip_real_phi_solve = false)
{
    const OphelieComplexEdgeFluxSolveReport report = solveOphelieComplexEdgeFluxWithCurrentA<ExecutionPolicy>(
        glass_body, inner, names, params, boundary_geom, dp, skip_real_phi_solve);
    if (solve_report != nullptr)
    {
        *solve_report = report;
    }
    execOphelieEdgeFluxPostPhiPipeline<ExecutionPolicy>(glass_body, inner, names, params, skip_real_phi_solve);
    syncOphelieEdgeReconToPrimaryEJQ<ExecutionPolicy>(glass_body, names, params);
    BaseParticles &particles = glass_body.getBaseParticles();
    return hostEdgeFluxReconPower(particles, names, particles.TotalRealParticles(), params);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_COMPONENT_H
