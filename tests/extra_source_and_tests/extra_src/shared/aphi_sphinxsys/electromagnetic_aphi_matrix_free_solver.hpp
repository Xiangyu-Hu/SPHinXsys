#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.h"
#if SPHINXSYS_USE_SYCL
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_sycl_kernels.hpp"
#endif
#include <limits>

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

namespace
{
inline Complex computeComplexFieldMeanForSolver(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return Complex(0.0, 0.0);
    }

    Complex sum(0.0, 0.0);
    for (const Complex &value : field)
    {
        sum += value;
    }
    return sum / static_cast<Real>(field.size());
}

inline void removeComplexFieldMeanForSolver(StdVec<Complex> &field)
{
    if (field.empty())
    {
        return;
    }

    const Complex mean_value = computeComplexFieldMeanForSolver(field);
    for (Complex &value : field)
    {
        value -= mean_value;
    }
}
} // namespace

inline ScalarComplexHelmholtzSolverParameters::ScalarComplexHelmholtzSolverParameters()
    : max_iterations_(200), relaxation_factor_(0.5), absolute_tolerance_(1.0e-8), diagonal_regularization_(1.0e-12)
{
}

inline ScalarComplexHelmholtzSolverState::ScalarComplexHelmholtzSolverState()
    : iterations_(0), converged_(false), initial_residual_l2_(0.0), current_residual_l2_(0.0), current_mean_abs_(0.0),
      current_max_abs_(0.0), current_min_diagonal_abs_(0.0), current_max_diagonal_abs_(0.0), current_nonfinite_diagonal_count_(0)
{
}

inline void applyScalarComplexHelmholtzJacobiUpdate(StdVec<Complex> &field,
                                                    const StdVec<Complex> &reaction_coefficient,
                                                    const ScalarComplexHelmholtzResiduals &residuals,
                                                    Real relaxation_factor, Real diagonal_regularization)
{
#if SPHINXSYS_USE_SYCL
    if (matrixFreeJacobiUseSycl())
    {
        applyScalarComplexHelmholtzJacobiUpdateSycl(field, reaction_coefficient, residuals, relaxation_factor,
                                                    diagonal_regularization);
        return;
    }
#endif
    const size_t number_of_particles = field.size();
    for (size_t index_i = 0; index_i != number_of_particles; ++index_i)
    {
        const Complex local_diagonal =
            Complex(residuals.diagonal_scale_[index_i] + diagonal_regularization, 0.0) + reaction_coefficient[index_i];

        if (std::abs(local_diagonal) <= diagonal_regularization)
        {
            continue;
        }

        field[index_i] -= relaxation_factor * residuals.residual_[index_i] / local_diagonal;
    }
}

template <class ResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarComplexHelmholtz(StdVec<Complex> &field, const StdVec<Complex> &rhs,
                                                                     const StdVec<Complex> &reaction_coefficient,
                                                                     ScalarComplexHelmholtzResiduals &residuals,
                                                                     const ScalarComplexHelmholtzSolverParameters &parameters,
                                                                     ResidualBuilder &&residual_builder)
{
    ScalarComplexHelmholtzSolverState state;

    if (field.size() != rhs.size() || field.size() != reaction_coefficient.size())
    {
        return state;
    }

    if (residuals.residual_.size() != field.size())
    {
        residuals.resize(field.size());
    }

    const bool track_last_stable = (parameters.residual_growth_limit_factor_ > Real(0.0)) ||
                                   parameters.remove_field_mean_after_jacobi_;
    StdVec<Complex> last_stable_field;
    Real last_stable_residual_l2 = MaxReal;
    if (track_last_stable)
    {
        last_stable_field = field;
    }

    for (size_t iteration = 0; iteration != parameters.max_iterations_; ++iteration)
    {
        residuals.clear();
        residual_builder(field, residuals);
        finalizeScalarHelmholtzResiduals(field, rhs, reaction_coefficient, residuals);

        state.iterations_ = iteration + 1;
        state.current_residual_l2_ = residuals.l2_norm_;
        state.current_mean_abs_ = residuals.mean_abs_;
        state.current_max_abs_ = residuals.max_abs_;
        state.current_min_diagonal_abs_ = std::numeric_limits<Real>::max();
        state.current_max_diagonal_abs_ = 0.0;
        state.current_nonfinite_diagonal_count_ = 0;
        for (size_t index_i = 0; index_i != field.size(); ++index_i)
        {
            const Complex local_diagonal =
                Complex(residuals.diagonal_scale_[index_i] + parameters.diagonal_regularization_, 0.0) + reaction_coefficient[index_i];
            const Real diagonal_abs = std::abs(local_diagonal);
            if (!std::isfinite(diagonal_abs))
            {
                ++state.current_nonfinite_diagonal_count_;
                continue;
            }
            state.current_min_diagonal_abs_ = std::min(state.current_min_diagonal_abs_, diagonal_abs);
            state.current_max_diagonal_abs_ = std::max(state.current_max_diagonal_abs_, diagonal_abs);
        }
        if (state.current_min_diagonal_abs_ == std::numeric_limits<Real>::max())
        {
            state.current_min_diagonal_abs_ = 0.0;
        }

        if (iteration == 0)
        {
            state.initial_residual_l2_ = residuals.l2_norm_;
        }

        if (!std::isfinite(residuals.l2_norm_))
        {
            if (track_last_stable)
            {
                field = last_stable_field;
                state.current_residual_l2_ = last_stable_residual_l2;
            }
            break;
        }

        if (parameters.residual_growth_limit_factor_ > Real(0.0) && last_stable_residual_l2 < MaxReal &&
            residuals.l2_norm_ >
                parameters.residual_growth_limit_factor_ * (last_stable_residual_l2 + TinyReal))
        {
            if (track_last_stable)
            {
                field = last_stable_field;
                state.current_residual_l2_ = last_stable_residual_l2;
            }
            break;
        }

        if (residuals.l2_norm_ <= parameters.absolute_tolerance_)
        {
            state.converged_ = true;
            break;
        }

        applyScalarComplexHelmholtzJacobiUpdate(field, reaction_coefficient, residuals, parameters.relaxation_factor_,
                                                parameters.diagonal_regularization_);

        if (parameters.remove_field_mean_after_jacobi_)
        {
            removeComplexFieldMeanForSolver(field);
        }

        bool field_is_finite = true;
        for (const Complex &value : field)
        {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                field_is_finite = false;
                break;
            }
        }
        if (!field_is_finite)
        {
            if (track_last_stable)
            {
                field = last_stable_field;
                state.current_residual_l2_ = last_stable_residual_l2;
            }
            break;
        }

        if (track_last_stable)
        {
            last_stable_field = field;
            last_stable_residual_l2 = residuals.l2_norm_;
        }
    }

    return state;
}

#if SPHINXSYS_USE_SYCL
namespace
{
bool g_matrix_free_jacobi_use_sycl = false;
} // namespace

inline void setMatrixFreeJacobiUseSycl(bool enable)
{
    g_matrix_free_jacobi_use_sycl = enable;
}

inline bool matrixFreeJacobiUseSycl()
{
    return g_matrix_free_jacobi_use_sycl;
}
#endif

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_SOLVER_HPP
