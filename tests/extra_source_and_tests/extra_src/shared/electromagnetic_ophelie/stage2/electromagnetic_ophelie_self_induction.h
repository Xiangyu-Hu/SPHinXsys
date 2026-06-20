#ifndef ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H
#define ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H

#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_phi_gmres.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_biot_savart.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real ophelieSelfInductionPicardPhiEqResTolerance(const OphelieParameters &params)
{
    if (ophelieUseEdgeFluxElectromotiveRhs(params) && params.edge_flux_complex_)
    {
        return params.self_induction_phi_eq_res_tolerance_;
    }
    return params.phi_eq_res_vol_gate_;
}

inline bool ophelieSelfInductionPicardConverged(Real j_rel, Real phi_eq_res_vol, const OphelieParameters &params)
{
    return j_rel < params.self_induction_j_tolerance_ &&
           phi_eq_res_vol < ophelieSelfInductionPicardPhiEqResTolerance(params);
}

inline Real ophelieSelfInductionPicardPhiEqResVol(const OphelieComplexEdgeFluxSolveReport &solve_report)
{
    return std::max(solve_report.phi_eq_res_vol_imag, solve_report.phi_eq_res_vol_real);
}

template <class ExecutionPolicy>
inline Real hostRelativeVecdFieldChange(BaseParticles &particles, const std::string &field_name,
                                        const StdVec<Real> &previous_components, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    const Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real delta_x = field[i][0] - previous_components[3 * i + 0];
        const Real delta_y = field[i][1] - previous_components[3 * i + 1];
        const Real delta_z = field[i][2] - previous_components[3 * i + 2];
        numerator += delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        denominator += field[i].squaredNorm();
    }
    return std::sqrt(numerator / (denominator + TinyReal));
}

template <class ExecutionPolicy>
inline void hostStoreVecdFieldComponents(BaseParticles &particles, const std::string &field_name,
                                         StdVec<Real> &components, size_t n)
{
    syncVariableToHost<Vecd>(particles, field_name);
    const Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    components.resize(3 * n);
    for (size_t i = 0; i < n; ++i)
    {
        components[3 * i + 0] = field[i][0];
        components[3 * i + 1] = field[i][1];
        components[3 * i + 2] = field[i][2];
    }
}

template <class ExecutionPolicy>
inline void hostRelaxVecdFieldTowardPrevious(BaseParticles &particles, const std::string &field_name,
                                             const StdVec<Real> &previous_components, Real relaxation_factor, size_t n)
{
    const Real alpha = std::max(Real(0.0), std::min(relaxation_factor, Real(1.0)));
    syncVariableToHost<Vecd>(particles, field_name);
    Vecd *field = particles.getVariableDataByName<Vecd>(field_name);
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd value_new = field[i];
        const Vecd value_old(previous_components[3 * i + 0], previous_components[3 * i + 1],
                             previous_components[3 * i + 2]);
        field[i] = alpha * value_new + (Real(1.0) - alpha) * value_old;
    }
    syncVariableToDevice<Vecd>(particles, field_name);
}

template <class ExecutionPolicy>
inline Real hostRelativeJImagChange(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                    const StdVec<Real> &previous_j_imag_components, size_t n)
{
    return hostRelativeVecdFieldChange<ExecutionPolicy>(particles, names.j_imag, previous_j_imag_components, n);
}

template <class ExecutionPolicy>
inline void hostStoreJImagComponents(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                   StdVec<Real> &components, size_t n)
{
    hostStoreVecdFieldComponents<ExecutionPolicy>(particles, names.j_imag, components, n);
}

template <class ExecutionPolicy>
inline void hostRelaxJImagTowardPrevious(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                         const StdVec<Real> &previous_j_imag_components, Real relaxation_factor,
                                         size_t n)
{
    hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, names.j_imag, previous_j_imag_components,
                                                     relaxation_factor, n);
}

template <class ExecutionPolicy>
inline Real hostRelativeComplexEdgeJChange(const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                         BaseParticles &particles, const StdVec<Real> &previous_j_real_components,
                                         const StdVec<Real> &previous_j_imag_components, size_t n)
{
    const std::string &j_real_field = getOphelieAIndJRealFieldName(names, params);
    const std::string &j_imag_field = getOphelieAIndJImagFieldName(names, params);
    const Real j_real_change =
        hostRelativeVecdFieldChange<ExecutionPolicy>(particles, j_real_field, previous_j_real_components, n);
    const Real j_imag_change =
        hostRelativeVecdFieldChange<ExecutionPolicy>(particles, j_imag_field, previous_j_imag_components, n);
    return std::max(j_real_change, j_imag_change);
}

template <class ExecutionPolicy>
inline Real hostRelativeComplexAIndChange(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                          const StdVec<Real> &previous_a_ind_real_components,
                                          const StdVec<Real> &previous_a_ind_imag_components, size_t n)
{
    const Real a_real_change = hostRelativeVecdFieldChange<ExecutionPolicy>(
        particles, names.a_ind_real, previous_a_ind_real_components, n);
    const Real a_imag_change = hostRelativeVecdFieldChange<ExecutionPolicy>(
        particles, names.a_ind_imag, previous_a_ind_imag_components, n);
    return std::max(a_real_change, a_imag_change);
}

template <class ExecutionPolicy>
inline Real runOphelieComplexEdgeFluxSelfInductionWithPhiSolve(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters &params,
    Real &phi_solver_rel_residual, size_t &self_induction_iterations_used, Real &phi_eq_res_vol_out,
    bool &picard_converged_out)
{
    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        glass_body, names);
    const std::string &j_real_source = getOphelieAIndJRealFieldName(names, params);
    const std::string &j_imag_source = getOphelieAIndJImagFieldName(names, params);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> compute_self_induced_biot(
        glass_body, names, params, j_real_source, j_imag_source);

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> previous_j_real_components;
    StdVec<Real> previous_j_imag_components;
    Real relative_j_change = 0.0;
    self_induction_iterations_used = 1;
    OphelieProgressLogger self_ind_progress("self_induction_complex");

    const bool saved_use_a_total = params.use_a_total_for_edge_flux_;
    params.use_a_total_for_edge_flux_ = true;

    combine_vector_potential.exec();
    self_ind_progress.log("initial combine A_src");
    OphelieComplexEdgeFluxSolveReport solve_report;
    (void)execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(glass_body, inner, names, params, nullptr, 0.0,
                                                                        &solve_report);
    phi_solver_rel_residual = solve_report.phi_imag_solver_rel_residual;
    phi_eq_res_vol_out = ophelieSelfInductionPicardPhiEqResVol(solve_report);
    self_ind_progress.log("initial complex edge-flux solve, phi_imag_rel_res=" +
                          std::to_string(phi_solver_rel_residual) + " phi_eq_res_vol=" +
                          std::to_string(phi_eq_res_vol_out));
    hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_real_source, previous_j_real_components, n);
    hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components, n);
    picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);

    for (size_t iteration = 1; iteration < params.self_induction_max_iterations_; ++iteration)
    {
        self_induction_iterations_used = iteration + 1;
        compute_self_induced_biot.exec();
        combine_vector_potential.exec();

        (void)execOphelieComplexEdgeFluxSolveReconAndPower<ExecutionPolicy>(glass_body, inner, names, params, nullptr,
                                                                            0.0, &solve_report);
        phi_solver_rel_residual = solve_report.phi_imag_solver_rel_residual;
        phi_eq_res_vol_out = ophelieSelfInductionPicardPhiEqResVol(solve_report);

        hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, j_real_source, previous_j_real_components,
                                                          params.self_induction_relaxation_factor_, n);
        hostRelaxVecdFieldTowardPrevious<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components,
                                                          params.self_induction_relaxation_factor_, n);

        relative_j_change = hostRelativeComplexEdgeJChange<ExecutionPolicy>(
            names, params, particles, previous_j_real_components, previous_j_imag_components, n);
        hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_real_source, previous_j_real_components, n);
        hostStoreVecdFieldComponents<ExecutionPolicy>(particles, j_imag_source, previous_j_imag_components, n);
        picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);
        self_ind_progress.log("outer iter " + std::to_string(self_induction_iterations_used) + " J_rel=" +
                              std::to_string(relative_j_change) + " phi_eq_res_vol=" +
                              std::to_string(phi_eq_res_vol_out) + " phi_imag_rel_res=" +
                              std::to_string(phi_solver_rel_residual) + " phi_real_rel_res=" +
                              std::to_string(solve_report.phi_real_solver_rel_residual) +
                              " picard_converged=" + std::to_string(picard_converged_out ? 1 : 0));
        if (picard_converged_out)
        {
            break;
        }
    }

    params.use_a_total_for_edge_flux_ = saved_use_a_total;
    picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);
    self_ind_progress.finish("final_J_rel=" + std::to_string(relative_j_change) + " phi_eq_res_vol=" +
                             std::to_string(phi_eq_res_vol_out) +
                             " picard_converged=" + std::to_string(picard_converged_out ? 1 : 0));
    return relative_j_change;
}

template <class ExecutionPolicy>
inline Real runOphelieSelfInductionWithPhiSolve(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                                 OphelieParameters &params, Real &phi_solver_rel_residual,
                                                 size_t &self_induction_iterations_used, Real &phi_eq_res_vol_out,
                                                 bool &picard_converged_out)
{
    if (ophelieUseEdgeFluxElectromotiveRhs(params) && params.edge_flux_complex_)
    {
        return runOphelieComplexEdgeFluxSelfInductionWithPhiSolve<ExecutionPolicy>(
            glass_body, inner, names, params, phi_solver_rel_residual, self_induction_iterations_used, phi_eq_res_vol_out,
            picard_converged_out);
    }

    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        glass_body, names);
    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> compute_self_induced_biot(
        glass_body, names, params);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(inner, names);
    StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, names, params);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> previous_j_imag_components;
    Real relative_j_change = 0.0;
    self_induction_iterations_used = 1;
    OphelieProgressLogger self_ind_progress("self_induction");

    combine_vector_potential.exec();
    self_ind_progress.log("initial combine A_src");
    phi_solver_rel_residual = solvePhiImag<ExecutionPolicy>(glass_body, inner, names, params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    phi_eq_res_vol_out = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);
    self_ind_progress.log("initial phi solve, rel_res=" + std::to_string(phi_solver_rel_residual) +
                          " phi_eq_res_vol=" + std::to_string(phi_eq_res_vol_out));
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    compute_ejq_with_phi.exec();
    hostStoreJImagComponents<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
    syncVariableToDevice<Vecd>(particles, names.j_imag);
    picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);

    for (size_t iteration = 1; iteration < params.self_induction_max_iterations_; ++iteration)
    {
        self_induction_iterations_used = iteration + 1;
        syncVariableToDevice<Vecd>(particles, names.j_imag);
        compute_self_induced_biot.exec();
        combine_vector_potential.exec();

        phi_solver_rel_residual = solvePhiImag<ExecutionPolicy>(glass_body, inner, names, params);
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        phi_eq_res_vol_out = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        compute_grad_phi.exec();
        compute_ejq_with_phi.exec();
        hostRelaxJImagTowardPrevious<ExecutionPolicy>(particles, names, previous_j_imag_components,
                                                      params.self_induction_relaxation_factor_, n);

        relative_j_change = hostRelativeJImagChange<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
        hostStoreJImagComponents<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
        picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);
        self_ind_progress.log("outer iter " + std::to_string(self_induction_iterations_used) + " J_rel=" +
                              std::to_string(relative_j_change) + " phi_eq_res_vol=" +
                              std::to_string(phi_eq_res_vol_out) + " phi_rel_res=" +
                              std::to_string(phi_solver_rel_residual) +
                              " picard_converged=" + std::to_string(picard_converged_out ? 1 : 0));
        if (picard_converged_out)
        {
            break;
        }
    }

    picard_converged_out = ophelieSelfInductionPicardConverged(relative_j_change, phi_eq_res_vol_out, params);
    self_ind_progress.finish("final_J_rel=" + std::to_string(relative_j_change) + " phi_eq_res_vol=" +
                             std::to_string(phi_eq_res_vol_out) +
                             " picard_converged=" + std::to_string(picard_converged_out ? 1 : 0));
    return relative_j_change;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H
