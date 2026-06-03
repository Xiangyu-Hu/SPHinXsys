#ifndef ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H
#define ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H

#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_phi_gmres.h"
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

template <class ExecutionPolicy>
inline Real hostRelativeJImagChange(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                    const StdVec<Real> &previous_j_imag_components, size_t n)
{
    syncVariableToHost<Vecd>(particles, names.j_imag);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(names.j_imag);
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real delta_x = j_imag[i][0] - previous_j_imag_components[3 * i + 0];
        const Real delta_y = j_imag[i][1] - previous_j_imag_components[3 * i + 1];
        const Real delta_z = j_imag[i][2] - previous_j_imag_components[3 * i + 2];
        numerator += delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
        denominator += j_imag[i].squaredNorm();
    }
    return std::sqrt(numerator / (denominator + TinyReal));
}

template <class ExecutionPolicy>
inline void hostStoreJImagComponents(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                     StdVec<Real> &components, size_t n)
{
    syncVariableToHost<Vecd>(particles, names.j_imag);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(names.j_imag);
    components.resize(3 * n);
    for (size_t i = 0; i < n; ++i)
    {
        components[3 * i + 0] = j_imag[i][0];
        components[3 * i + 1] = j_imag[i][1];
        components[3 * i + 2] = j_imag[i][2];
    }
}

template <class ExecutionPolicy>
inline void hostRelaxJImagTowardPrevious(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                          const StdVec<Real> &previous_j_imag_components, Real relaxation_factor,
                                          size_t n)
{
    const Real alpha = std::max(Real(0.0), std::min(relaxation_factor, Real(1.0)));
    syncVariableToHost<Vecd>(particles, names.j_imag);
    Vecd *j_imag = particles.getVariableDataByName<Vecd>(names.j_imag);
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd j_new = j_imag[i];
        const Vecd j_old(previous_j_imag_components[3 * i + 0], previous_j_imag_components[3 * i + 1],
                         previous_j_imag_components[3 * i + 2]);
        j_imag[i] = alpha * j_new + (Real(1.0) - alpha) * j_old;
    }
    syncVariableToDevice<Vecd>(particles, names.j_imag);
}

template <class ExecutionPolicy>
inline Real runOphelieSelfInductionWithPhiSolve(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                                 const OphelieParameters &params, Real &phi_solver_rel_residual,
                                                 size_t &self_induction_iterations_used)
{
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
    self_ind_progress.log("initial phi solve, rel_res=" + std::to_string(phi_solver_rel_residual));
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    compute_ejq_with_phi.exec();
    hostStoreJImagComponents<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
    syncVariableToDevice<Vecd>(particles, names.j_imag);

    for (size_t iteration = 1; iteration < params.self_induction_max_iterations_; ++iteration)
    {
        self_induction_iterations_used = iteration + 1;
        syncVariableToDevice<Vecd>(particles, names.j_imag);
        compute_self_induced_biot.exec();
        combine_vector_potential.exec();
        // Do not sync a_src/b_src to device here: combine writes on device; host buffers are stale.

        phi_solver_rel_residual = solvePhiImag<ExecutionPolicy>(glass_body, inner, names, params);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        compute_grad_phi.exec();
        compute_ejq_with_phi.exec();
        hostRelaxJImagTowardPrevious<ExecutionPolicy>(particles, names, previous_j_imag_components,
                                                      params.self_induction_relaxation_factor_, n);

        relative_j_change = hostRelativeJImagChange<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
        hostStoreJImagComponents<ExecutionPolicy>(particles, names, previous_j_imag_components, n);
        self_ind_progress.log("outer iter " + std::to_string(self_induction_iterations_used) + " J_rel=" +
                              std::to_string(relative_j_change) + " phi_rel_res=" + std::to_string(phi_solver_rel_residual));
        if (relative_j_change < params.self_induction_j_tolerance_)
        {
            break;
        }
    }

    self_ind_progress.finish("final_J_rel=" + std::to_string(relative_j_change));
    return relative_j_change;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_SELF_INDUCTION_H
