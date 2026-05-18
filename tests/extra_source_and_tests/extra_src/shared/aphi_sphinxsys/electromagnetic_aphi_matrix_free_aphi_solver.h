#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_residuals.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct MatrixFreeAPhiSphNativeContext;

struct MatrixFreeAPhiSolverParameters
{
    size_t max_outer_iterations;
    ScalarComplexHelmholtzSolverParameters a_component_solver;
    ScalarComplexHelmholtzSolverParameters phi_solver;
    ScalarComplexHelmholtzSolverParameters gauge_solver;
    Real residual_tolerance;
    Real divergence_tolerance;
    bool enable_gauge_projection;
    bool remove_gauge_mean_offset;
    bool use_operator_consistent_gauge_projection;
    bool enable_gauge_penalty;
    Real gauge_penalty_coefficient;
    size_t gauge_penalty_ramp_iterations;
    Real gauge_penalty_initial_ratio;
    Real outer_relaxation_factor;
    Real update_tolerance;
    size_t residual_growth_guard_start_iteration;
    Real residual_growth_limit;
    size_t source_ramp_iterations;
    Real source_initial_ratio;

    MatrixFreeAPhiSolverParameters();
};

struct MatrixFreeAPhiGaugeStepDiagnostics
{
    bool applied;
    Real chi_l2;
    Real grad_chi_l2;
    Real chi_max_abs;
    Real phi_reference_offset_after_phi_solve_abs;
    Real phi_reference_offset_after_gauge_update_abs;
    Real phi_reference_offset_after_final_reference_abs;
    Real gauge_rhs_mean_abs_before_projection;
    Real gauge_rhs_mean_abs_after_projection;
    Real chi_mean_abs_after_solve;
    Real divergence_a_before_l2;
    Real divergence_a_after_raw_l2;
    Real divergence_a_after_final_l2;
    Real divergence_j_before_l2;
    Real divergence_j_after_raw_l2;
    Real divergence_j_after_final_l2;
    Real electric_field_before_l2;
    Real electric_field_after_raw_l2;
    Real electric_field_after_final_l2;
    Real electric_field_change_raw_l2;
    Real electric_field_change_final_l2;
    Real current_density_before_l2;
    Real current_density_after_raw_l2;
    Real current_density_after_final_l2;
    Real current_density_change_raw_l2;
    Real current_density_change_final_l2;

    MatrixFreeAPhiGaugeStepDiagnostics();
};

struct MatrixFreeAPhiGaugeProjectionResult
{
    ScalarComplexHelmholtzSolverState solver_state;
    MatrixFreeAPhiGaugeStepDiagnostics diagnostics;

    MatrixFreeAPhiGaugeProjectionResult();
};

enum class MatrixFreeAPhiScalarHelmholtzBackendKind
{
    GraphHostControlled,
    GraphFusedSycl,
    SphNativeHostControlled,
    ReservedSphNativeSpecialized
};

struct MatrixFreeAPhiSolverOperatorSemantics
{
    MatrixFreeAPhiDiscreteOperatorSemanticsKind particle_kernel;
    MatrixFreeAPhiDiscreteOperatorSemanticsKind legacy_solver;
    MatrixFreeAPhiDiscreteOperatorSemanticsKind gauge;
    MatrixFreeAPhiDiscreteOperatorSemanticsKind gauge_diagnostic_field;
    MatrixFreeAPhiDiscreteOperatorSemanticsKind gauge_current_divergence;

    bool usesParticleKernel() const
    {
        return particle_kernel == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel;
    }
};

struct MatrixFreeAPhiWorkspace
{
    MatrixFreeAPhiFields previous_fields;
    MatrixFreeAPhiSources effective_sources;
    StdVec<Vec3c> sigma_grad_phi;
    StdVec<Vec3c> gauge_penalty_gradient;
    StdVec<Complex> rhs_ax;
    StdVec<Complex> rhs_ay;
    StdVec<Complex> rhs_az;
    StdVec<Complex> divergence_sigma_a;
    StdVec<Complex> rhs_phi;

    void resize(size_t size);
};

struct MatrixFreeAPhiGaugeEvaluationBundle
{
    StdVec<Complex> divergence_a;
    StdVec<Vec3c> electric_field;
    StdVec<Vec3c> current_density;
    StdVec<Complex> divergence_j;
};

struct MatrixFreeAPhiStaggeredOperatorBundle
{
    StdVec<Vec3c> sigma_grad_phi;
    StdVec<Vec3c> gauge_penalty_gradient;

    void resize(size_t size);
};

struct MatrixFreeAPhiSolverState
{
    size_t outer_iterations;
    bool converged;
    ScalarComplexHelmholtzSolverState ax_state;
    ScalarComplexHelmholtzSolverState ay_state;
    ScalarComplexHelmholtzSolverState az_state;
    ScalarComplexHelmholtzSolverState phi_state;
    ScalarComplexHelmholtzSolverState gauge_state;
    MatrixFreeAPhiGaugeStepDiagnostics gauge_diagnostics;
    MatrixFreeAPhiResiduals residuals;
    Real effective_gauge_penalty_coefficient;
    Real effective_source_scale;
    Real field_update_l2;
    Real relative_field_update_l2;
    MatrixFreeAPhiIterationHistory iteration_history;

    MatrixFreeAPhiSolverState();
};

void enforcePhiReference(MatrixFreeAPhiFields &fields, const MatrixFreeAPhiParameters &parameters);

MatrixFreeAPhiGaugeProjectionResult applyMatrixFreeAPhiGaugeProjectionStep(
    const MatrixFreePairwiseGraph &graph, MatrixFreeAPhiFields &fields, const StdVec<Real> &electrical_conductivity,
    const MatrixFreeAPhiParameters &parameters, const ScalarComplexHelmholtzSolverParameters &gauge_solver_parameters,
    bool remove_gauge_mean_offset, bool enforce_phi_reference_after_gauge, bool use_operator_consistent_gauge_projection,
    MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr);

MatrixFreeAPhiScalarHelmholtzBackendKind selectMatrixFreeAPhiScalarHelmholtzBackend(
    const MatrixFreePairwiseGraph &graph, size_t field_size, bool prefer_sph_native, bool require_flat_graph_compatibility,
    MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr);

/** Extension point for future device-resident Krylov / multigrid; current production path remains host-controlled iteration
 *  with SYCL-accelerated operator applies (Jacobi / matrix-free graph kernels). */
enum class MatrixFreeAPhiLinearAlgebraBackendKind
{
    HostControlledWithSyclPrimitives,
    ReservedDeviceKrylov
};

#if SPHINXSYS_USE_SYCL
MatrixFreeAPhiLinearAlgebraBackendKind matrixFreeAPhiLinearAlgebraBackend();
void setMatrixFreeAPhiLinearAlgebraBackend(MatrixFreeAPhiLinearAlgebraBackendKind kind);
/** Dynamics / drivers call this once per substep (or after graph rebuild): SYCL graph + scalar pools + optional USM topology mirror. */
void matrixFreeAPhiSyclBindTimestepResources(const MatrixFreePairwiseGraph &graph, size_t particle_count);
#else
inline MatrixFreeAPhiLinearAlgebraBackendKind matrixFreeAPhiLinearAlgebraBackend()
{
    return MatrixFreeAPhiLinearAlgebraBackendKind::HostControlledWithSyclPrimitives;
}
inline void setMatrixFreeAPhiLinearAlgebraBackend(MatrixFreeAPhiLinearAlgebraBackendKind) {}
inline void matrixFreeAPhiSyclBindTimestepResources(const MatrixFreePairwiseGraph &, size_t) {}
#endif

#if SPHINXSYS_USE_SYCL
struct MatrixFreeSyclPolicySnapshot
{
    bool gradient = false;
    bool harmonic_gradient = false;
    bool laplace_residual = false;
    bool jacobi = false;
};

MatrixFreeSyclPolicySnapshot captureMatrixFreeSyclPolicySnapshot();
void applyMatrixFreeSyclPolicySnapshot(const MatrixFreeSyclPolicySnapshot &snapshot);
#else
struct MatrixFreeSyclPolicySnapshot
{
    bool gradient = false;
    bool harmonic_gradient = false;
    bool laplace_residual = false;
    bool jacobi = false;
};

inline MatrixFreeSyclPolicySnapshot captureMatrixFreeSyclPolicySnapshot()
{
    return MatrixFreeSyclPolicySnapshot{};
}

inline void applyMatrixFreeSyclPolicySnapshot(const MatrixFreeSyclPolicySnapshot &) {}
#endif

MatrixFreeAPhiSolverState solveMatrixFreeAPhiStaggered(const MatrixFreePairwiseGraph &graph,
                                                       MatrixFreeAPhiFields &fields,
                                                       const StdVec<Real> &electrical_conductivity,
                                                       const StdVec<Real> &magnetic_reluctivity,
                                                       const MatrixFreeAPhiParameters &parameters,
                                                       const MatrixFreeAPhiSources &sources,
                                                       const MatrixFreeAPhiSolverParameters &solver_parameters,
                                                       MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.hpp"

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_H
