#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct ScalarComplexHelmholtzSolverParameters
{
    size_t max_iterations_;
    Real relaxation_factor_;
    Real absolute_tolerance_;
    Real diagonal_regularization_;
    /** If positive, abort iteration when residual L2 exceeds this factor times the last stable pre-update L2 (gauge-style guard). */
    Real residual_growth_limit_factor_ = Real(0.0);
    /** After each Jacobi update, pull the field to host, subtract its mean, and re-upload (gauge null-space step). */
    bool remove_field_mean_after_jacobi_ = false;

    ScalarComplexHelmholtzSolverParameters();
};

struct ScalarComplexHelmholtzSolverState
{
    size_t iterations_;
    bool converged_;
    Real initial_residual_l2_;
    Real current_residual_l2_;
    Real current_mean_abs_;
    Real current_max_abs_;
    Real current_min_diagonal_abs_;
    Real current_max_diagonal_abs_;
    size_t current_nonfinite_diagonal_count_;

    ScalarComplexHelmholtzSolverState();
};

void applyScalarComplexHelmholtzJacobiUpdate(StdVec<Complex> &field,
                                             const StdVec<Complex> &reaction_coefficient,
                                             const ScalarComplexHelmholtzResiduals &residuals,
                                             Real relaxation_factor, Real diagonal_regularization);

#if SPHINXSYS_USE_SYCL
void setMatrixFreeJacobiUseSycl(bool enable);
bool matrixFreeJacobiUseSycl();
size_t matrixFreeJacobiReactionUploadCount();
size_t matrixFreeFusedGraphHelmholtzSolveCount();
size_t matrixFreeFusedGraphHelmholtzMetricDownloadCount();
size_t matrixFreeFusedGraphHelmholtzMetricScalarDownloadCount();
size_t matrixFreeFusedGraphHelmholtzRhsUploadCount();
size_t matrixFreeFusedGraphHelmholtzReactionUploadCount();
#else
inline void setMatrixFreeJacobiUseSycl(bool) {}
inline bool matrixFreeJacobiUseSycl()
{
    return false;
}
inline size_t matrixFreeJacobiReactionUploadCount()
{
    return 0;
}
inline size_t matrixFreeFusedGraphHelmholtzSolveCount()
{
    return 0;
}
inline size_t matrixFreeFusedGraphHelmholtzMetricDownloadCount()
{
    return 0;
}
inline size_t matrixFreeFusedGraphHelmholtzMetricScalarDownloadCount()
{
    return 0;
}
inline size_t matrixFreeFusedGraphHelmholtzRhsUploadCount()
{
    return 0;
}
inline size_t matrixFreeFusedGraphHelmholtzReactionUploadCount()
{
    return 0;
}
#endif

template <class ResidualBuilder>
ScalarComplexHelmholtzSolverState solveScalarComplexHelmholtz(StdVec<Complex> &field, const StdVec<Complex> &rhs,
                                                              const StdVec<Complex> &reaction_coefficient,
                                                              ScalarComplexHelmholtzResiduals &residuals,
                                                              const ScalarComplexHelmholtzSolverParameters &parameters,
                                                              ResidualBuilder &&residual_builder);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.hpp"

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_H
