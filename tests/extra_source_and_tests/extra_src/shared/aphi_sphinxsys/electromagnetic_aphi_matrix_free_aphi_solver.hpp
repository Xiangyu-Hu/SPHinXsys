#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"
#include <algorithm>
#if SPHINXSYS_USE_SYCL
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_sycl_queue.hpp"
#endif

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

#if SPHINXSYS_USE_SYCL
namespace
{
inline MatrixFreeAPhiLinearAlgebraBackendKind &matrixFreeAPhiLinearAlgebraBackendStorage()
{
    static MatrixFreeAPhiLinearAlgebraBackendKind storage =
        MatrixFreeAPhiLinearAlgebraBackendKind::HostControlledWithSyclPrimitives;
    return storage;
}
} // namespace

inline MatrixFreeSyclPolicySnapshot captureMatrixFreeSyclPolicySnapshot()
{
    MatrixFreeSyclPolicySnapshot snapshot;
    snapshot.gradient = matrixFreeGradientUseSycl();
    snapshot.harmonic_gradient = matrixFreeHarmonicGradientUseSycl();
    snapshot.laplace_residual = matrixFreeLaplaceResidualUseSycl();
    snapshot.jacobi = matrixFreeJacobiUseSycl();
    return snapshot;
}

inline void applyMatrixFreeSyclPolicySnapshot(const MatrixFreeSyclPolicySnapshot &snapshot)
{
    setMatrixFreeGradientUseSycl(snapshot.gradient);
    setMatrixFreeHarmonicGradientUseSycl(snapshot.harmonic_gradient);
    setMatrixFreeLaplaceResidualUseSycl(snapshot.laplace_residual);
    setMatrixFreeJacobiUseSycl(snapshot.jacobi);
}

inline MatrixFreeAPhiLinearAlgebraBackendKind matrixFreeAPhiLinearAlgebraBackend()
{
    return matrixFreeAPhiLinearAlgebraBackendStorage();
}
inline void setMatrixFreeAPhiLinearAlgebraBackend(MatrixFreeAPhiLinearAlgebraBackendKind kind)
{
    matrixFreeAPhiLinearAlgebraBackendStorage() = kind;
}

inline void matrixFreeAPhiSyclBindTimestepResources(const MatrixFreePairwiseGraph &graph, size_t particle_count)
{
    const bool any_operator = matrixFreeGradientUseSycl() || matrixFreeHarmonicGradientUseSycl() ||
                              matrixFreeLaplaceResidualUseSycl() || matrixFreeJacobiUseSycl();
    const bool reserve_device_la =
        matrixFreeAPhiLinearAlgebraBackend() == MatrixFreeAPhiLinearAlgebraBackendKind::ReservedDeviceKrylov;
    if (!any_operator && !reserve_device_la)
    {
        return;
    }
    prepareMatrixFreeAPhiSyclDeviceResources(graph, particle_count);
}
#endif

namespace
{
inline StdVec<Complex> zeroComplexField(size_t size)
{
    return StdVec<Complex>(size, Complex(0.0, 0.0));
}

Complex removeMeanOffset(StdVec<Complex> &field);

#if SPHINXSYS_USE_SYCL
inline ScalarComplexHelmholtzSolverState solveScalarComplexHelmholtzFromGraphSycl(
    const MatrixFreePairwiseGraph &graph, StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient,
    const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, bool reaction_coefficient_is_zero);
#endif

template <class ResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarComplexHelmholtzHostControlled(
    StdVec<Complex> &field, const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, ResidualBuilder &&residual_builder)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(field.size());
    return solveScalarComplexHelmholtz(field, rhs, reaction_coefficient, residuals, parameters,
                                       std::forward<ResidualBuilder>(residual_builder));
}

template <class NativeResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarLaplaceHelmholtzSphNativeHostControlled(
    StdVec<Complex> &field, const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, NativeResidualBuilder &&native_residual_builder)
{
    return solveScalarComplexHelmholtzHostControlled(field, rhs, reaction_coefficient, parameters,
                                                     std::forward<NativeResidualBuilder>(native_residual_builder));
}

template <class NativeResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarLaplaceHelmholtzSphNativeSpecialized(
    StdVec<Complex> &field, const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, NativeResidualBuilder &&native_residual_builder)
{
    // First specialized SPH-native stage: keep the validated native residual
    // build and host-controlled iteration loop, but route the Jacobi update
    // through the device primitive for this backend only.
#if SPHINXSYS_USE_SYCL
    const MatrixFreeSyclPolicySnapshot snapshot = captureMatrixFreeSyclPolicySnapshot();
    MatrixFreeSyclPolicySnapshot specialized_snapshot = snapshot;
    specialized_snapshot.jacobi = true;
    applyMatrixFreeSyclPolicySnapshot(specialized_snapshot);

    const ScalarComplexHelmholtzSolverState state = solveScalarLaplaceHelmholtzSphNativeHostControlled(
        field, rhs, reaction_coefficient, parameters, std::forward<NativeResidualBuilder>(native_residual_builder));

    applyMatrixFreeSyclPolicySnapshot(snapshot);
    return state;
#else
    return solveScalarLaplaceHelmholtzSphNativeHostControlled(
        field, rhs, reaction_coefficient, parameters, std::forward<NativeResidualBuilder>(native_residual_builder));
#endif
}

template <class GraphResidualBuilder, class NativeResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarLaplaceHelmholtzWithSelectedBackend(
    MatrixFreeAPhiScalarHelmholtzBackendKind backend, const MatrixFreePairwiseGraph &graph, StdVec<Complex> &field,
    const StdVec<Real> &diffusion_coefficient, const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, bool fused_graph_remove_mean,
    GraphResidualBuilder &&graph_residual_builder, NativeResidualBuilder &&native_residual_builder)
{
    if (backend == MatrixFreeAPhiScalarHelmholtzBackendKind::SphNativeHostControlled)
    {
        return solveScalarLaplaceHelmholtzSphNativeHostControlled(
            field, rhs, reaction_coefficient, parameters, std::forward<NativeResidualBuilder>(native_residual_builder));
    }
    if (backend == MatrixFreeAPhiScalarHelmholtzBackendKind::ReservedSphNativeSpecialized)
    {
        return solveScalarLaplaceHelmholtzSphNativeSpecialized(
            field, rhs, reaction_coefficient, parameters, std::forward<NativeResidualBuilder>(native_residual_builder));
    }
#if SPHINXSYS_USE_SYCL
    if (backend == MatrixFreeAPhiScalarHelmholtzBackendKind::GraphFusedSycl)
    {
        return solveScalarComplexHelmholtzFromGraphSycl(graph, field, diffusion_coefficient, rhs, reaction_coefficient,
                                                        parameters, fused_graph_remove_mean);
    }
#endif
    return solveScalarComplexHelmholtzHostControlled(field, rhs, reaction_coefficient, parameters,
                                                     std::forward<GraphResidualBuilder>(graph_residual_builder));
}

template <class ResidualBuilder>
inline ScalarComplexHelmholtzSolverState solveScalarHelmholtzWithNormalizedResidualStep(
    StdVec<Complex> &field, const StdVec<Complex> &rhs, const ScalarComplexHelmholtzSolverParameters &parameters,
    ResidualBuilder &&residual_builder)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(field.size());
    const StdVec<Complex> zero_reaction = zeroComplexField(field.size());

    ScalarComplexHelmholtzSolverState state;
    StdVec<Complex> last_stable_field = field;
    Real last_stable_residual_l2 = MaxReal;
    const Real residual_growth_limit = Real(5.0);

    for (size_t iteration = 0; iteration != parameters.max_iterations_; ++iteration)
    {
        residuals.clear();
        residual_builder(field, residuals);
        finalizeScalarHelmholtzResiduals(field, rhs, zero_reaction, residuals);

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
                Complex(residuals.diagonal_scale_[index_i] + parameters.diagonal_regularization_, 0.0);
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
            field = last_stable_field;
            state.current_residual_l2_ = last_stable_residual_l2;
            break;
        }

        if (last_stable_residual_l2 < MaxReal &&
            residuals.l2_norm_ > residual_growth_limit * (last_stable_residual_l2 + TinyReal))
        {
            field = last_stable_field;
            state.current_residual_l2_ = last_stable_residual_l2;
            break;
        }

        if (residuals.l2_norm_ <= parameters.absolute_tolerance_)
        {
            state.converged_ = true;
            break;
        }

        Real maximum_residual_abs = TinyReal;
        for (const Complex &value : residuals.residual_)
        {
            maximum_residual_abs = std::max(maximum_residual_abs, std::abs(value));
        }
        const Real normalized_step =
            parameters.relaxation_factor_ / std::max(maximum_residual_abs, TinyReal);

        StdVec<Complex> candidate_field = field;
        for (size_t index_i = 0; index_i != field.size(); ++index_i)
        {
            candidate_field[index_i] -= normalized_step * residuals.residual_[index_i];
        }
        removeMeanOffset(candidate_field);

        bool candidate_is_finite = true;
        for (const Complex &value : candidate_field)
        {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                candidate_is_finite = false;
                break;
            }
        }
        if (!candidate_is_finite)
        {
            field = last_stable_field;
            state.current_residual_l2_ = last_stable_residual_l2;
            break;
        }

        field = candidate_field;
        last_stable_field = field;
        last_stable_residual_l2 = residuals.l2_norm_;
    }

    return state;
}

#if SPHINXSYS_USE_SYCL
struct MatrixFreeSyclGraphHelmholtzBuffers
{
    size_t size_ = 0;
    size_t solve_count_ = 0;
    size_t metric_download_count_ = 0;
    size_t metric_scalar_download_count_ = 0;
    size_t rhs_upload_count_ = 0;
    size_t reaction_upload_count_ = 0;
    const Complex *rhs_source_ = nullptr;
    size_t rhs_source_size_ = 0;
    const Complex *reaction_source_ = nullptr;
    size_t reaction_source_size_ = 0;
    bool reaction_cache_is_zero_ = false;
    std::unique_ptr<sycl::buffer<Real, 1>> rhs_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> rhs_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> reaction_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> reaction_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> residual_abs_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> diagonal_abs_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> nonfinite_diagonal_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> metric_scalar_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> partial_metric_buffer_;
    size_t partial_group_count_ = 0;

    void resize(size_t size)
    {
        if (size_ == size && rhs_real_buffer_ && rhs_imag_buffer_ && reaction_real_buffer_ && reaction_imag_buffer_ &&
            residual_abs_buffer_ && diagonal_abs_buffer_ && nonfinite_diagonal_buffer_ && metric_scalar_buffer_)
        {
            return;
        }

        size_ = size;
        rhs_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        rhs_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        reaction_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        reaction_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        residual_abs_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        diagonal_abs_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        nonfinite_diagonal_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        metric_scalar_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(6));
        partial_metric_buffer_.reset();
        partial_group_count_ = 0;
        rhs_source_ = nullptr;
        rhs_source_size_ = 0;
        reaction_source_ = nullptr;
        reaction_source_size_ = 0;
        reaction_cache_is_zero_ = false;
    }

    void loadInputsIfNeeded(const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
                            bool reaction_coefficient_is_zero)
    {
        if (rhs_source_ != rhs.data() || rhs_source_size_ != rhs.size())
        {
            sycl::host_accessor<Real, 1, sycl::access::mode::write> rhs_real(*rhs_real_buffer_);
            sycl::host_accessor<Real, 1, sycl::access::mode::write> rhs_imag(*rhs_imag_buffer_);
            for (size_t i = 0; i != rhs.size(); ++i)
            {
                rhs_real[i] = rhs[i].real();
                rhs_imag[i] = rhs[i].imag();
            }
            rhs_source_ = rhs.data();
            rhs_source_size_ = rhs.size();
            ++rhs_upload_count_;
        }

        if (reaction_coefficient_is_zero)
        {
            if (!reaction_cache_is_zero_ || reaction_source_size_ != reaction_coefficient.size())
            {
                sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_real(*reaction_real_buffer_);
                sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_imag(*reaction_imag_buffer_);
                for (size_t i = 0; i != reaction_coefficient.size(); ++i)
                {
                    reaction_real[i] = 0.0;
                    reaction_imag[i] = 0.0;
                }
                reaction_source_ = nullptr;
                reaction_source_size_ = reaction_coefficient.size();
                reaction_cache_is_zero_ = true;
                ++reaction_upload_count_;
            }
            return;
        }

        if (reaction_cache_is_zero_ || reaction_source_ != reaction_coefficient.data() ||
            reaction_source_size_ != reaction_coefficient.size())
        {
            sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_real(*reaction_real_buffer_);
            sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_imag(*reaction_imag_buffer_);
            for (size_t i = 0; i != reaction_coefficient.size(); ++i)
            {
                reaction_real[i] = reaction_coefficient[i].real();
                reaction_imag[i] = reaction_coefficient[i].imag();
            }
            reaction_source_ = reaction_coefficient.data();
            reaction_source_size_ = reaction_coefficient.size();
            reaction_cache_is_zero_ = false;
            ++reaction_upload_count_;
        }
    }

    void aggregateMetricsOnDevice()
    {
        const size_t size = size_;
        const size_t local_size = 256;
        const size_t group_count = size == 0 ? size_t(1) : (size + local_size - 1) / local_size;
        const size_t global_size = group_count * local_size;
        if (!partial_metric_buffer_ || partial_group_count_ != group_count)
        {
            partial_metric_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(group_count * 6));
            partial_group_count_ = group_count;
        }

        matrixFreeSyclSubmit([&](sycl::handler &cgh) {
            auto residual_abs = residual_abs_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto diagonal_abs = diagonal_abs_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto nonfinite_diagonal = nonfinite_diagonal_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto partial_metrics = partial_metric_buffer_->get_access<sycl::access::mode::write>(cgh);
            sycl::local_accessor<Real, 1> local_residual_square_sum(sycl::range<1>(local_size), cgh);
            sycl::local_accessor<Real, 1> local_residual_abs_sum(sycl::range<1>(local_size), cgh);
            sycl::local_accessor<Real, 1> local_max_residual_abs(sycl::range<1>(local_size), cgh);
            sycl::local_accessor<Real, 1> local_min_diagonal_abs(sycl::range<1>(local_size), cgh);
            sycl::local_accessor<Real, 1> local_max_diagonal_abs(sycl::range<1>(local_size), cgh);
            sycl::local_accessor<Real, 1> local_nonfinite_diagonal_count(sycl::range<1>(local_size), cgh);
            const Real max_real = MaxReal;

            cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>(global_size), sycl::range<1>(local_size)),
                             [=](sycl::nd_item<1> item) {
                const size_t global_id = item.get_global_linear_id();
                const size_t local_id = item.get_local_linear_id();
                const size_t group_id = item.get_group_linear_id();

                Real residual_square = 0.0;
                Real residual_sum = 0.0;
                Real residual_max = 0.0;
                Real diagonal_min = max_real;
                Real diagonal_max = 0.0;
                Real nonfinite_count = 0.0;
                if (global_id < size)
                {
                    const Real local_residual_abs = residual_abs[global_id];
                    residual_square = local_residual_abs * local_residual_abs;
                    residual_sum = local_residual_abs;
                    residual_max = local_residual_abs;
                    if (nonfinite_diagonal[global_id] != 0.0)
                    {
                        nonfinite_count = 1.0;
                    }
                    else
                    {
                        diagonal_min = diagonal_abs[global_id];
                        diagonal_max = diagonal_abs[global_id];
                    }
                }

                local_residual_square_sum[local_id] = residual_square;
                local_residual_abs_sum[local_id] = residual_sum;
                local_max_residual_abs[local_id] = residual_max;
                local_min_diagonal_abs[local_id] = diagonal_min;
                local_max_diagonal_abs[local_id] = diagonal_max;
                local_nonfinite_diagonal_count[local_id] = nonfinite_count;
                item.barrier(sycl::access::fence_space::local_space);

                for (size_t stride = local_size / 2; stride != 0; stride /= 2)
                {
                    if (local_id < stride)
                    {
                        local_residual_square_sum[local_id] += local_residual_square_sum[local_id + stride];
                        local_residual_abs_sum[local_id] += local_residual_abs_sum[local_id + stride];
                        local_max_residual_abs[local_id] =
                            local_max_residual_abs[local_id] > local_max_residual_abs[local_id + stride]
                                ? local_max_residual_abs[local_id]
                                : local_max_residual_abs[local_id + stride];
                        local_min_diagonal_abs[local_id] =
                            local_min_diagonal_abs[local_id] < local_min_diagonal_abs[local_id + stride]
                                ? local_min_diagonal_abs[local_id]
                                : local_min_diagonal_abs[local_id + stride];
                        local_max_diagonal_abs[local_id] =
                            local_max_diagonal_abs[local_id] > local_max_diagonal_abs[local_id + stride]
                                ? local_max_diagonal_abs[local_id]
                                : local_max_diagonal_abs[local_id + stride];
                        local_nonfinite_diagonal_count[local_id] += local_nonfinite_diagonal_count[local_id + stride];
                    }
                    item.barrier(sycl::access::fence_space::local_space);
                }

                if (local_id == 0)
                {
                    const size_t offset = group_id * 6;
                    partial_metrics[offset + 0] = local_residual_square_sum[0];
                    partial_metrics[offset + 1] = local_residual_abs_sum[0];
                    partial_metrics[offset + 2] = local_max_residual_abs[0];
                    partial_metrics[offset + 3] = local_min_diagonal_abs[0];
                    partial_metrics[offset + 4] = local_max_diagonal_abs[0];
                    partial_metrics[offset + 5] = local_nonfinite_diagonal_count[0];
                }
            });
        });

        matrixFreeSyclSubmit([&](sycl::handler &cgh) {
            auto partial_metrics = partial_metric_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto metric_scalars = metric_scalar_buffer_->get_access<sycl::access::mode::write>(cgh);
            const size_t groups = group_count;
            const Real max_real = MaxReal;

            cgh.single_task([=]() {
                Real residual_square_sum = 0.0;
                Real residual_abs_sum = 0.0;
                Real max_residual_abs = 0.0;
                Real min_diagonal_abs = max_real;
                Real max_diagonal_abs = 0.0;
                Real nonfinite_diagonal_count = 0.0;
                for (size_t group = 0; group != groups; ++group)
                {
                    const size_t offset = group * 6;
                    residual_square_sum += partial_metrics[offset + 0];
                    residual_abs_sum += partial_metrics[offset + 1];
                    max_residual_abs =
                        partial_metrics[offset + 2] > max_residual_abs ? partial_metrics[offset + 2] : max_residual_abs;
                    min_diagonal_abs =
                        partial_metrics[offset + 3] < min_diagonal_abs ? partial_metrics[offset + 3] : min_diagonal_abs;
                    max_diagonal_abs =
                        partial_metrics[offset + 4] > max_diagonal_abs ? partial_metrics[offset + 4] : max_diagonal_abs;
                    nonfinite_diagonal_count += partial_metrics[offset + 5];
                }

                metric_scalars[0] = residual_square_sum;
                metric_scalars[1] = residual_abs_sum;
                metric_scalars[2] = max_residual_abs;
                metric_scalars[3] = min_diagonal_abs == max_real ? Real(0.0) : min_diagonal_abs;
                metric_scalars[4] = max_diagonal_abs;
                metric_scalars[5] = nonfinite_diagonal_count;
            });
        });
        matrixFreeSyclQueueWait();
    }

    void readMetrics(ScalarComplexHelmholtzSolverState &state)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::read> metric_scalars(*metric_scalar_buffer_);

        if (size_ != 0)
        {
            state.current_residual_l2_ = std::sqrt(metric_scalars[0] / static_cast<Real>(size_));
            state.current_mean_abs_ = metric_scalars[1] / static_cast<Real>(size_);
        }
        else
        {
            state.current_residual_l2_ = 0.0;
            state.current_mean_abs_ = 0.0;
        }
        state.current_max_abs_ = metric_scalars[2];
        state.current_min_diagonal_abs_ = metric_scalars[3];
        state.current_max_diagonal_abs_ = metric_scalars[4];
        state.current_nonfinite_diagonal_count_ = static_cast<size_t>(metric_scalars[5]);
        ++metric_download_count_;
        metric_scalar_download_count_ += 6;
    }
};

inline MatrixFreeSyclGraphHelmholtzBuffers &matrixFreeSyclGraphHelmholtzBuffers()
{
    static MatrixFreeSyclGraphHelmholtzBuffers buffers;
    return buffers;
}

inline void readMatrixFreeSyclScalarField(StdVec<Complex> &field, MatrixFreeSyclScalarValueBuffers &value_buffers)
{
    value_buffers.syncHostFieldBuffersFromDeviceUsm();
    sycl::host_accessor<Real, 1, sycl::access::mode::read> field_real(*value_buffers.field_real_buffer_);
    sycl::host_accessor<Real, 1, sycl::access::mode::read> field_imag(*value_buffers.field_imag_buffer_);
    for (size_t i = 0; i != field.size(); ++i)
    {
        field[i] = Complex(field_real[i], field_imag[i]);
    }
    value_buffers.markFieldCurrentForSource(field);
}

inline ScalarComplexHelmholtzSolverState solveScalarComplexHelmholtzFromGraphSycl(
    const MatrixFreePairwiseGraph &graph, StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient,
    const StdVec<Complex> &rhs, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzSolverParameters &parameters, bool reaction_coefficient_is_zero = false)
{
    ScalarComplexHelmholtzSolverState state;
    const size_t number_of_rows = graph.rows_.size();
    if (!matrixFreeLaplaceResidualUseSycl() || !matrixFreeJacobiUseSycl() ||
        !isFlatGraphFieldCompatible(graph, field.size()) || diffusion_coefficient.size() != number_of_rows ||
        rhs.size() != number_of_rows || reaction_coefficient.size() != number_of_rows)
    {
        return state;
    }

    MatrixFreeAPhiSyclWorkspace &sycl_workspace = matrixFreeAPhiSyclWorkspace();
    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    MatrixFreeSyclGraphHelmholtzBuffers &helmholtz_buffers = matrixFreeSyclGraphHelmholtzBuffers();

    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    value_buffers.resize(number_of_rows);
    value_buffers.loadFieldIfNeeded(field);
    value_buffers.loadCoefficientIfNeeded(diffusion_coefficient);
    helmholtz_buffers.resize(number_of_rows);
    helmholtz_buffers.loadInputsIfNeeded(rhs, reaction_coefficient, reaction_coefficient_is_zero);
    ++helmholtz_buffers.solve_count_;

    const bool track_last_stable = (parameters.residual_growth_limit_factor_ > Real(0.0)) ||
                                   parameters.remove_field_mean_after_jacobi_;
    StdVec<Complex> last_stable_field;
    Real last_stable_residual_l2 = MaxReal;
    if (track_last_stable)
    {
        last_stable_field = field;
    }

    const Real *helm_field_usm_r = value_buffers.deviceFieldRealUsmRead();
    const Real *helm_field_usm_i = value_buffers.deviceFieldImagUsmRead();
    const bool helm_field_reads_from_usm = helm_field_usm_r != nullptr && helm_field_usm_i != nullptr;

    for (size_t iteration = 0; iteration != parameters.max_iterations_; ++iteration)
    {
        if (graph_buffers.deviceTopologyUsmReady())
        {
            const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
            const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
            const Real *edge_w_dev = graph_buffers.deviceTopologyDiffusionWeights();
            if (helm_field_reads_from_usm)
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto rhs_real = helmholtz_buffers.rhs_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto rhs_imag = helmholtz_buffers.rhs_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_abs = helmholtz_buffers.residual_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal_abs = helmholtz_buffers.diagonal_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto nonfinite_diagonal =
                        helmholtz_buffers.nonfinite_diagonal_buffer_->get_access<sycl::access::mode::write>(cgh);
                    const Real reg = parameters.diagonal_regularization_;
                    const Real *const f_fr = helm_field_usm_r;
                    const Real *const f_fi = helm_field_usm_i;

                    cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        const Real value_i_real = f_fr[index_i];
                        const Real value_i_imag = f_fi[index_i];
                        Real laplace_real = 0.0;
                        Real laplace_imag = 0.0;
                        Real diagonal_scale = 0.0;
                        for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                        {
                            const size_t index_j = col_dev[edge];
                            const Real coefficient_ij =
                                2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                                (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                            const Real pair_weight = coefficient_ij * edge_w_dev[edge];
                            laplace_real += pair_weight * (value_i_real - f_fr[index_j]);
                            laplace_imag += pair_weight * (value_i_imag - f_fi[index_j]);
                            diagonal_scale += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                        }

                        const Real reaction_term_real =
                            reaction_real[index_i] * value_i_real - reaction_imag[index_i] * value_i_imag;
                        const Real reaction_term_imag =
                            reaction_real[index_i] * value_i_imag + reaction_imag[index_i] * value_i_real;
                        const Real local_residual_real = laplace_real + reaction_term_real - rhs_real[index_i];
                        const Real local_residual_imag = laplace_imag + reaction_term_imag - rhs_imag[index_i];
                        residual_real[index_i] = local_residual_real;
                        residual_imag[index_i] = local_residual_imag;
                        diagonal[index_i] = diagonal_scale;

                        const Real local_diagonal_real = diagonal_scale + reg + reaction_real[index_i];
                        const Real local_diagonal_imag = reaction_imag[index_i];
                        const Real local_diagonal_abs = sycl::sqrt(local_diagonal_real * local_diagonal_real +
                                                                   local_diagonal_imag * local_diagonal_imag);
                        residual_abs[index_i] = sycl::sqrt(local_residual_real * local_residual_real +
                                                           local_residual_imag * local_residual_imag);
                        diagonal_abs[index_i] = local_diagonal_abs;
                        nonfinite_diagonal[index_i] = sycl::isfinite(local_diagonal_abs) ? Real(0.0) : Real(1.0);
                    });
                });
            }
            else
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto rhs_real = helmholtz_buffers.rhs_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto rhs_imag = helmholtz_buffers.rhs_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_abs = helmholtz_buffers.residual_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal_abs = helmholtz_buffers.diagonal_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto nonfinite_diagonal =
                        helmholtz_buffers.nonfinite_diagonal_buffer_->get_access<sycl::access::mode::write>(cgh);
                    const Real reg = parameters.diagonal_regularization_;

                    cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        const Real value_i_real = f_real[index_i];
                        const Real value_i_imag = f_imag[index_i];
                        Real laplace_real = 0.0;
                        Real laplace_imag = 0.0;
                        Real diagonal_scale = 0.0;
                        for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                        {
                            const size_t index_j = col_dev[edge];
                            const Real coefficient_ij =
                                2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                                (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                            const Real pair_weight = coefficient_ij * edge_w_dev[edge];
                            laplace_real += pair_weight * (value_i_real - f_real[index_j]);
                            laplace_imag += pair_weight * (value_i_imag - f_imag[index_j]);
                            diagonal_scale += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                        }

                        const Real reaction_term_real =
                            reaction_real[index_i] * value_i_real - reaction_imag[index_i] * value_i_imag;
                        const Real reaction_term_imag =
                            reaction_real[index_i] * value_i_imag + reaction_imag[index_i] * value_i_real;
                        const Real local_residual_real = laplace_real + reaction_term_real - rhs_real[index_i];
                        const Real local_residual_imag = laplace_imag + reaction_term_imag - rhs_imag[index_i];
                        residual_real[index_i] = local_residual_real;
                        residual_imag[index_i] = local_residual_imag;
                        diagonal[index_i] = diagonal_scale;

                        const Real local_diagonal_real = diagonal_scale + reg + reaction_real[index_i];
                        const Real local_diagonal_imag = reaction_imag[index_i];
                        const Real local_diagonal_abs = sycl::sqrt(local_diagonal_real * local_diagonal_real +
                                                                   local_diagonal_imag * local_diagonal_imag);
                        residual_abs[index_i] =
                            sycl::sqrt(local_residual_real * local_residual_real + local_residual_imag * local_residual_imag);
                        diagonal_abs[index_i] = local_diagonal_abs;
                        nonfinite_diagonal[index_i] = sycl::isfinite(local_diagonal_abs) ? Real(0.0) : Real(1.0);
                    });
                });
            }
        }
        else
        {
            if (helm_field_reads_from_usm)
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto edge_weights = graph_buffers.diffusion_weights_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto rhs_real = helmholtz_buffers.rhs_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto rhs_imag = helmholtz_buffers.rhs_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_abs = helmholtz_buffers.residual_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal_abs = helmholtz_buffers.diagonal_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto nonfinite_diagonal =
                        helmholtz_buffers.nonfinite_diagonal_buffer_->get_access<sycl::access::mode::write>(cgh);
                    const Real reg = parameters.diagonal_regularization_;
                    const Real *const f_fr = helm_field_usm_r;
                    const Real *const f_fi = helm_field_usm_i;

                    cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        const Real value_i_real = f_fr[index_i];
                        const Real value_i_imag = f_fi[index_i];
                        Real laplace_real = 0.0;
                        Real laplace_imag = 0.0;
                        Real diagonal_scale = 0.0;
                        for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                        {
                            const size_t index_j = column_indices[edge];
                            const Real coefficient_ij =
                                2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                                (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                            const Real pair_weight = coefficient_ij * edge_weights[edge];
                            laplace_real += pair_weight * (value_i_real - f_fr[index_j]);
                            laplace_imag += pair_weight * (value_i_imag - f_fi[index_j]);
                            diagonal_scale += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                        }

                        const Real reaction_term_real =
                            reaction_real[index_i] * value_i_real - reaction_imag[index_i] * value_i_imag;
                        const Real reaction_term_imag =
                            reaction_real[index_i] * value_i_imag + reaction_imag[index_i] * value_i_real;
                        const Real local_residual_real = laplace_real + reaction_term_real - rhs_real[index_i];
                        const Real local_residual_imag = laplace_imag + reaction_term_imag - rhs_imag[index_i];
                        residual_real[index_i] = local_residual_real;
                        residual_imag[index_i] = local_residual_imag;
                        diagonal[index_i] = diagonal_scale;

                        const Real local_diagonal_real = diagonal_scale + reg + reaction_real[index_i];
                        const Real local_diagonal_imag = reaction_imag[index_i];
                        const Real local_diagonal_abs = sycl::sqrt(local_diagonal_real * local_diagonal_real +
                                                                   local_diagonal_imag * local_diagonal_imag);
                        residual_abs[index_i] = sycl::sqrt(local_residual_real * local_residual_real +
                                                           local_residual_imag * local_residual_imag);
                        diagonal_abs[index_i] = local_diagonal_abs;
                        nonfinite_diagonal[index_i] = sycl::isfinite(local_diagonal_abs) ? Real(0.0) : Real(1.0);
                    });
                });
            }
            else
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto edge_weights = graph_buffers.diffusion_weights_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto rhs_real = helmholtz_buffers.rhs_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto rhs_imag = helmholtz_buffers.rhs_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto residual_abs = helmholtz_buffers.residual_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto diagonal_abs = helmholtz_buffers.diagonal_abs_buffer_->get_access<sycl::access::mode::write>(cgh);
                    auto nonfinite_diagonal =
                        helmholtz_buffers.nonfinite_diagonal_buffer_->get_access<sycl::access::mode::write>(cgh);
                    const Real reg = parameters.diagonal_regularization_;

                    cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        const Real value_i_real = f_real[index_i];
                        const Real value_i_imag = f_imag[index_i];
                        Real laplace_real = 0.0;
                        Real laplace_imag = 0.0;
                        Real diagonal_scale = 0.0;
                        for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                        {
                            const size_t index_j = column_indices[edge];
                            const Real coefficient_ij =
                                2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                                (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                            const Real pair_weight = coefficient_ij * edge_weights[edge];
                            laplace_real += pair_weight * (value_i_real - f_real[index_j]);
                            laplace_imag += pair_weight * (value_i_imag - f_imag[index_j]);
                            diagonal_scale += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                        }

                        const Real reaction_term_real =
                            reaction_real[index_i] * value_i_real - reaction_imag[index_i] * value_i_imag;
                        const Real reaction_term_imag =
                            reaction_real[index_i] * value_i_imag + reaction_imag[index_i] * value_i_real;
                        const Real local_residual_real = laplace_real + reaction_term_real - rhs_real[index_i];
                        const Real local_residual_imag = laplace_imag + reaction_term_imag - rhs_imag[index_i];
                        residual_real[index_i] = local_residual_real;
                        residual_imag[index_i] = local_residual_imag;
                        diagonal[index_i] = diagonal_scale;

                        const Real local_diagonal_real = diagonal_scale + reg + reaction_real[index_i];
                        const Real local_diagonal_imag = reaction_imag[index_i];
                        const Real local_diagonal_abs = sycl::sqrt(local_diagonal_real * local_diagonal_real +
                                                                   local_diagonal_imag * local_diagonal_imag);
                        residual_abs[index_i] =
                            sycl::sqrt(local_residual_real * local_residual_real + local_residual_imag * local_residual_imag);
                        diagonal_abs[index_i] = local_diagonal_abs;
                        nonfinite_diagonal[index_i] = sycl::isfinite(local_diagonal_abs) ? Real(0.0) : Real(1.0);
                    });
                });
            }
        }

        state.iterations_ = iteration + 1;
        helmholtz_buffers.aggregateMetricsOnDevice();
        helmholtz_buffers.readMetrics(state);
        if (iteration == 0)
        {
            state.initial_residual_l2_ = state.current_residual_l2_;
        }
        const Real residual_l2_before_jacobi = state.current_residual_l2_;

        if (!std::isfinite(residual_l2_before_jacobi))
        {
            if (track_last_stable)
            {
                field = last_stable_field;
                value_buffers.invalidateFieldCache();
                value_buffers.loadFieldIfNeeded(field);
                state.current_residual_l2_ = last_stable_residual_l2;
            }
            break;
        }

        if (parameters.residual_growth_limit_factor_ > Real(0.0) && last_stable_residual_l2 < MaxReal &&
            residual_l2_before_jacobi >
                parameters.residual_growth_limit_factor_ * (last_stable_residual_l2 + TinyReal))
        {
            field = last_stable_field;
            value_buffers.invalidateFieldCache();
            value_buffers.loadFieldIfNeeded(field);
            state.current_residual_l2_ = last_stable_residual_l2;
            break;
        }

        if (residual_l2_before_jacobi <= parameters.absolute_tolerance_)
        {
            state.converged_ = true;
            break;
        }

        Real *field_re_usm = value_buffers.mutableDeviceFieldRealUsmForJacobi();
        Real *field_im_usm = value_buffers.mutableDeviceFieldImagUsmForJacobi();
        if (field_re_usm != nullptr && field_im_usm != nullptr)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                const Real relax = parameters.relaxation_factor_;
                const Real reg = parameters.diagonal_regularization_;
                Real *const f_re = field_re_usm;
                Real *const f_im = field_im_usm;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t i = id[0];
                    const Real ld_re = diagonal[i] + reg + reaction_real[i];
                    const Real ld_im = reaction_imag[i];
                    const Real denom = ld_re * ld_re + ld_im * ld_im;
                    if (sycl::sqrt(denom) <= reg || denom <= 0.0)
                    {
                        return;
                    }
                    const Real update_real = (residual_real[i] * ld_re + residual_imag[i] * ld_im) / denom;
                    const Real update_imag = (residual_imag[i] * ld_re - residual_real[i] * ld_im) / denom;
                    f_re[i] -= relax * update_real;
                    f_im[i] -= relax * update_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto residual_real = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto residual_imag = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto reaction_real = helmholtz_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto reaction_imag = helmholtz_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                const Real relax = parameters.relaxation_factor_;
                const Real reg = parameters.diagonal_regularization_;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t i = id[0];
                    const Real ld_re = diagonal[i] + reg + reaction_real[i];
                    const Real ld_im = reaction_imag[i];
                    const Real denom = ld_re * ld_re + ld_im * ld_im;
                    if (sycl::sqrt(denom) <= reg || denom <= 0.0)
                    {
                        return;
                    }
                    const Real update_real = (residual_real[i] * ld_re + residual_imag[i] * ld_im) / denom;
                    const Real update_imag = (residual_imag[i] * ld_re - residual_real[i] * ld_im) / denom;
                    f_real[i] -= relax * update_real;
                    f_imag[i] -= relax * update_imag;
                });
            });
        }

        if (track_last_stable)
        {
            readMatrixFreeSyclScalarField(field, value_buffers);
            if (parameters.remove_field_mean_after_jacobi_)
            {
                removeMeanOffset(field);
                value_buffers.invalidateFieldCache();
                value_buffers.loadFieldIfNeeded(field);
            }
            last_stable_field = field;
            last_stable_residual_l2 = residual_l2_before_jacobi;
        }
    }

    readMatrixFreeSyclScalarField(field, value_buffers);
    return state;
}
#endif

inline StdVec<Complex> buildReactionCoefficient(const StdVec<Real> &electrical_conductivity, Real angular_frequency)
{
    const Complex imaginary_unit(0.0, 1.0);
    StdVec<Complex> coefficient(electrical_conductivity.size(), Complex(0.0, 0.0));
    for (size_t i = 0; i != electrical_conductivity.size(); ++i)
    {
        coefficient[i] = imaginary_unit * angular_frequency * electrical_conductivity[i];
    }
    return coefficient;
}

inline StdVec<Complex> computeDivergenceOfSigmaA(const MatrixFreePairwiseGraph &graph,
                                                 const MatrixFreeAPhiFields &fields,
                                                 const StdVec<Real> &electrical_conductivity,
                                                 MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                                 MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                                     MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeHarmonicDivergenceOfVectorField(
            fields.ax, fields.ay, fields.az, electrical_conductivity);
    }
    return computeHarmonicDivergenceOfVectorField(graph, fields.ax, fields.ay, fields.az, electrical_conductivity);
}

inline MatrixFreeAPhiSolverOperatorSemantics selectMatrixFreeAPhiSolverOperatorSemantics(
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const MatrixFreeAPhiDiscreteOperatorSemanticsKind particle_kernel =
        useMatrixFreeAPhiSphNativePath(sph_native_context)
            ? MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel
            : MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility;
    return {particle_kernel, MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility,
            MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility, particle_kernel, particle_kernel};
}

inline StdVec<Vec3c> applyMatrixFreeScalarGradientWithSemantics(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
    MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics, MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeStandardGradient(field);
    }
    return applyMatrixFreeGradient(graph, field);
}

inline StdVec<Complex> computeGradientDivergenceOfVectorFieldWithSemantics(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
    const StdVec<Complex> &field_z, MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeGradientDivergenceOfVectorField(field_x, field_y, field_z);
    }
    return computeGradientDivergenceOfVectorField(graph, field_x, field_y, field_z);
}

inline void accumulateScalarDivergenceOfGradientHelmholtzResidualsWithSemantics(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field, ScalarComplexHelmholtzResiduals &residuals,
    MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics, MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        sph_native_context->operatorAssembly().accumulateScalarDivergenceOfGradientHelmholtzResiduals(field, residuals);
        return;
    }
    accumulateScalarDivergenceOfGradientResidualsFromGraph(graph, field, residuals);
}

inline StdVec<Complex> computeDivergenceOfA(const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
                                            MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                            MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                                MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    return computeGradientDivergenceOfVectorFieldWithSemantics(
        graph, fields.ax, fields.ay, fields.az, semantics, sph_native_context);
}

inline StdVec<Vec3c> computeGaugePenaltyGradient(const MatrixFreePairwiseGraph &graph,
                                                 const MatrixFreeAPhiFields &fields,
                                                 MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                                 MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                                     MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    const size_t n = fields.ax.size();
    StdVec<Vec3c> gradient(n, Vec3c::Zero());
#if SPHINXSYS_USE_SYCL
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility &&
        matrixFreeGradientUseSycl() && graph.flat_.isValidFor(graph.rows_.size()) && graph.rows_.size() == n &&
        tryMatrixFreeDivergenceAndGradientOfDivergenceOfVectorPotentialSycl(graph, fields.ax, fields.ay, fields.az, nullptr,
                                                                           gradient))
    {
        return gradient;
    }
#endif
    const StdVec<Complex> divergence_a = computeDivergenceOfA(graph, fields, sph_native_context, semantics);
    return applyMatrixFreeScalarGradientWithSemantics(graph, divergence_a, semantics, sph_native_context);
}

inline Complex computePhiReferenceOffset(const MatrixFreeAPhiFields &fields, const MatrixFreeAPhiParameters &parameters)
{
    if (!parameters.fix_phi_reference || fields.phi.empty() || parameters.phi_reference_index >= fields.phi.size())
    {
        return Complex(0.0, 0.0);
    }
    return fields.phi[parameters.phi_reference_index] - parameters.phi_reference_value;
}

inline Complex computeComplexFieldMean(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return Complex(0.0, 0.0);
    }
    Complex mean_value(0.0, 0.0);
    for (const Complex &value : field)
    {
        mean_value += value;
    }
    return mean_value / static_cast<Real>(field.size());
}

inline Complex removeMeanOffset(StdVec<Complex> &field)
{
    if (field.empty())
    {
        return Complex(0.0, 0.0);
    }
    const Complex mean_value = computeComplexFieldMean(field);
    for (Complex &value : field)
    {
        value -= mean_value;
    }
    return mean_value;
}

inline Real computeComplexFieldL2Norm(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (const Complex &value : field)
    {
        squared_sum += std::norm(value);
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

inline Real computeComplexFieldMaxAbs(const StdVec<Complex> &field)
{
    Real max_abs = 0.0;
    for (const Complex &value : field)
    {
        max_abs = SMAX(max_abs, std::abs(value));
    }
    return max_abs;
}

inline Real computeVectorComplexFieldL2Norm(const StdVec<Vec3c> &field)
{
    if (field.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (const Vec3c &value : field)
    {
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            squared_sum += std::norm(value[axis]);
        }
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

inline Real computeVectorComplexFieldL2Difference(const StdVec<Vec3c> &lhs, const StdVec<Vec3c> &rhs)
{
    if (lhs.size() != rhs.size() || lhs.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != lhs.size(); ++i)
    {
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            squared_sum += std::norm(lhs[i][axis] - rhs[i][axis]);
        }
    }
    return std::sqrt(squared_sum / static_cast<Real>(lhs.size()));
}


inline Real computeAPhiFieldL2Norm(const MatrixFreeAPhiFields &fields)
{
    const size_t number_of_particles = fields.ax.size();
    if (number_of_particles == 0)
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        squared_sum += std::norm(fields.ax[i]);
        squared_sum += std::norm(fields.ay[i]);
        squared_sum += std::norm(fields.az[i]);
        squared_sum += std::norm(fields.phi[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(number_of_particles));
}

inline Real computeAPhiFieldUpdateL2(const MatrixFreeAPhiFields &previous_fields,
                                     const MatrixFreeAPhiFields &current_fields)
{
    const size_t number_of_particles = previous_fields.ax.size();
    if (number_of_particles == 0)
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        squared_sum += std::norm(current_fields.ax[i] - previous_fields.ax[i]);
        squared_sum += std::norm(current_fields.ay[i] - previous_fields.ay[i]);
        squared_sum += std::norm(current_fields.az[i] - previous_fields.az[i]);
        squared_sum += std::norm(current_fields.phi[i] - previous_fields.phi[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(number_of_particles));
}

inline void blendAPhiFields(const MatrixFreeAPhiFields &previous_fields,
                            MatrixFreeAPhiFields &current_fields,
                            Real relaxation_factor)
{
    const Real alpha = SMAX(Real(0.0), SMIN(Real(1.0), relaxation_factor));
    if (alpha >= Real(1.0))
    {
        return;
    }
    for (size_t i = 0; i != current_fields.ax.size(); ++i)
    {
        current_fields.ax[i] = previous_fields.ax[i] + alpha * (current_fields.ax[i] - previous_fields.ax[i]);
        current_fields.ay[i] = previous_fields.ay[i] + alpha * (current_fields.ay[i] - previous_fields.ay[i]);
        current_fields.az[i] = previous_fields.az[i] + alpha * (current_fields.az[i] - previous_fields.az[i]);
        current_fields.phi[i] = previous_fields.phi[i] + alpha * (current_fields.phi[i] - previous_fields.phi[i]);
    }
#if SPHINXSYS_USE_SYCL
    invalidateMatrixFreeAPhiSyclValueFieldCache();
#endif
}

inline StdVec<Vec3c> applyMatrixFreeScalarGradient(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                   MatrixFreeAPhiSphNativeContext *sph_native_context,
                                                   MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                                       MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    return applyMatrixFreeScalarGradientWithSemantics(graph, field, semantics, sph_native_context);
}

inline void accumulateScalarDivergenceOfGradientHelmholtzResiduals(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field, ScalarComplexHelmholtzResiduals &residuals,
    MatrixFreeAPhiSphNativeContext *sph_native_context,
    MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
        MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    accumulateScalarDivergenceOfGradientHelmholtzResidualsWithSemantics(graph, field, residuals, semantics,
                                                                        sph_native_context);
}

inline StdVec<Vec3c> computeElectricField(const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
                                          const MatrixFreeAPhiParameters &parameters,
                                          MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                          MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                              MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    const size_t number_of_particles = fields.ax.size();
    const Complex imaginary_unit(0.0, 1.0);
    const StdVec<Vec3c> grad_phi = applyMatrixFreeScalarGradient(graph, fields.phi, sph_native_context, semantics);
    StdVec<Vec3c> electric_field(number_of_particles, Vec3c::Zero());
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        electric_field[i][0] = -imaginary_unit * parameters.angular_frequency * fields.ax[i] - grad_phi[i][0];
        electric_field[i][1] = -imaginary_unit * parameters.angular_frequency * fields.ay[i] - grad_phi[i][1];
        electric_field[i][2] = -imaginary_unit * parameters.angular_frequency * fields.az[i] - grad_phi[i][2];
    }
    return electric_field;
}

inline StdVec<Vec3c> computeCurrentDensity(const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
                                           const StdVec<Real> &electrical_conductivity,
                                           const MatrixFreeAPhiParameters &parameters,
                                           MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                           MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                               MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    const StdVec<Vec3c> electric_field = computeElectricField(graph, fields, parameters, sph_native_context, semantics);
    StdVec<Vec3c> current_density(electric_field.size(), Vec3c::Zero());
    for (size_t i = 0; i != electric_field.size(); ++i)
    {
        current_density[i] = electrical_conductivity[i] * electric_field[i];
    }
    return current_density;
}

inline StdVec<Complex> computeDivergenceOfCurrentDensity(const MatrixFreePairwiseGraph &graph,
                                                         const StdVec<Vec3c> &current_density,
                                                         MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr,
                                                         MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics =
                                                             MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility)
{
    const size_t number_of_particles = current_density.size();
    StdVec<Complex> current_x(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> current_y(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> current_z(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        current_x[i] = current_density[i][0];
        current_y[i] = current_density[i][1];
        current_z[i] = current_density[i][2];
    }
    return computeGradientDivergenceOfVectorFieldWithSemantics(
        graph, current_x, current_y, current_z, semantics, sph_native_context);
}

inline MatrixFreeAPhiGaugeEvaluationBundle evaluateMatrixFreeAPhiGaugeFields(
    const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
    const StdVec<Real> &electrical_conductivity, const MatrixFreeAPhiParameters &parameters,
    MatrixFreeAPhiSphNativeContext *sph_native_context,
    const MatrixFreeAPhiSolverOperatorSemantics &semantics)
{
    MatrixFreeAPhiGaugeEvaluationBundle bundle;
    bundle.divergence_a = computeDivergenceOfA(graph, fields, sph_native_context, semantics.gauge);
    bundle.electric_field =
        computeElectricField(graph, fields, parameters, sph_native_context, semantics.gauge_diagnostic_field);
    bundle.current_density = computeCurrentDensity(graph, fields, electrical_conductivity, parameters,
                                                   sph_native_context, semantics.gauge_diagnostic_field);
    bundle.divergence_j =
        computeDivergenceOfCurrentDensity(graph, bundle.current_density, sph_native_context,
                                          semantics.gauge_current_divergence);
    return bundle;
}

inline MatrixFreeAPhiStaggeredOperatorBundle assembleMatrixFreeAPhiStaggeredOperators(
    const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
    const StdVec<Real> &electrical_conductivity, bool enable_gauge_penalty,
    const MatrixFreeAPhiSolverOperatorSemantics &semantics, MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    MatrixFreeAPhiStaggeredOperatorBundle bundle;
    bundle.resize(fields.ax.size());

    if (semantics.usesParticleKernel())
    {
        bundle.sigma_grad_phi =
            sph_native_context->operatorAssembly().computeSigmaGradPhi(fields.phi, electrical_conductivity);
    }
    else
    {
        bundle.sigma_grad_phi =
            applyMatrixFreeHarmonicWeightedGradient(graph, fields.phi, electrical_conductivity);
    }

    if (enable_gauge_penalty)
    {
        bundle.gauge_penalty_gradient =
            computeGaugePenaltyGradient(graph, fields, sph_native_context, semantics.legacy_solver);
    }
    else
    {
        std::fill(bundle.gauge_penalty_gradient.begin(), bundle.gauge_penalty_gradient.end(), Vec3c::Zero());
    }
    return bundle;
}

inline void scaleMatrixFreeAPhiSources(const MatrixFreeAPhiSources &sources, Real effective_source_scale,
                                       MatrixFreeAPhiSources &effective_sources)
{
    for (size_t i = 0; i != effective_sources.source_ax.size(); ++i)
    {
        effective_sources.source_ax[i] = effective_source_scale * sources.source_ax[i];
        effective_sources.source_ay[i] = effective_source_scale * sources.source_ay[i];
        effective_sources.source_az[i] = effective_source_scale * sources.source_az[i];
        effective_sources.source_phi[i] = effective_source_scale * sources.source_phi[i];
    }
}

inline Real computeMatrixFreeAPhiEffectiveGaugePenalty(const MatrixFreeAPhiSolverParameters &solver_parameters,
                                                       size_t outer_iteration)
{
    if (!solver_parameters.enable_gauge_penalty)
    {
        return 0.0;
    }

    const Real initial_ratio =
        SMAX(Real(0.0), SMIN(Real(1.0), solver_parameters.gauge_penalty_initial_ratio));
    if (solver_parameters.gauge_penalty_ramp_iterations == 0)
    {
        return solver_parameters.gauge_penalty_coefficient;
    }

    const Real ramp_fraction = SMIN(Real(1.0), static_cast<Real>(outer_iteration + 1) /
                                                 static_cast<Real>(solver_parameters.gauge_penalty_ramp_iterations));
    const Real blended_ratio = initial_ratio + (1.0 - initial_ratio) * ramp_fraction;
    return solver_parameters.gauge_penalty_coefficient * blended_ratio;
}

inline Real computeMatrixFreeAPhiEffectiveSourceScale(const MatrixFreeAPhiSolverParameters &solver_parameters,
                                                      size_t outer_iteration)
{
    const Real source_initial_ratio =
        SMAX(Real(0.0), SMIN(Real(1.0), solver_parameters.source_initial_ratio));
    if (solver_parameters.source_ramp_iterations == 0)
    {
        return 1.0;
    }

    const Real source_ramp_fraction =
        SMIN(Real(1.0), static_cast<Real>(outer_iteration + 1) /
                             static_cast<Real>(solver_parameters.source_ramp_iterations));
    return source_initial_ratio + (1.0 - source_initial_ratio) * source_ramp_fraction;
}

inline void assembleMatrixFreeAPhiARhs(const MatrixFreeAPhiSources &effective_sources,
                                       const MatrixFreeAPhiStaggeredOperatorBundle &operators,
                                       Real effective_penalty, MatrixFreeAPhiWorkspace &workspace)
{
    for (size_t i = 0; i != workspace.rhs_ax.size(); ++i)
    {
        workspace.rhs_ax[i] = effective_sources.source_ax[i] - operators.sigma_grad_phi[i][0] -
                              effective_penalty * operators.gauge_penalty_gradient[i][0];
        workspace.rhs_ay[i] = effective_sources.source_ay[i] - operators.sigma_grad_phi[i][1] -
                              effective_penalty * operators.gauge_penalty_gradient[i][1];
        workspace.rhs_az[i] = effective_sources.source_az[i] - operators.sigma_grad_phi[i][2] -
                              effective_penalty * operators.gauge_penalty_gradient[i][2];
    }
}

inline void assembleMatrixFreeAPhiPhiRhs(const MatrixFreeAPhiSources &effective_sources,
                                         const StdVec<Complex> &divergence_sigma_a,
                                         const MatrixFreeAPhiParameters &parameters,
                                         StdVec<Complex> &rhs_phi)
{
    const Complex imaginary_unit(0.0, 1.0);
    for (size_t i = 0; i != rhs_phi.size(); ++i)
    {
        rhs_phi[i] = imaginary_unit * parameters.angular_frequency * divergence_sigma_a[i] +
                     effective_sources.source_phi[i];
    }
}

inline void finalizeMatrixFreeAPhiFieldUpdateMetrics(const MatrixFreeAPhiFields &previous_fields,
                                                     const MatrixFreeAPhiFields &current_fields,
                                                     MatrixFreeAPhiSolverState &state)
{
    state.field_update_l2 = computeAPhiFieldUpdateL2(previous_fields, current_fields);
    state.relative_field_update_l2 = state.field_update_l2 / SMAX(computeAPhiFieldL2Norm(current_fields), TinyReal);
}

inline void appendMatrixFreeAPhiIterationRecord(MatrixFreeAPhiSolverState &state)
{
    MatrixFreeAPhiIterationRecord iteration_record;
    iteration_record.outer_iteration_ = state.outer_iterations;
    iteration_record.effective_source_scale_ = state.effective_source_scale;
    iteration_record.field_update_l2_ = state.field_update_l2;
    iteration_record.relative_field_update_l2_ = state.relative_field_update_l2;
    iteration_record.residual_ax_l2_ = state.residuals.residual_ax_l2;
    iteration_record.residual_ay_l2_ = state.residuals.residual_ay_l2;
    iteration_record.residual_az_l2_ = state.residuals.residual_az_l2;
    iteration_record.residual_phi_l2_ = state.residuals.residual_phi_l2;
    iteration_record.divergence_a_l2_ = state.residuals.divergence_a_l2;
    iteration_record.divergence_j_l2_ = state.residuals.divergence_j_l2;
    state.iteration_history.push_back(iteration_record);
}

inline Real computeMatrixFreeAPhiMaxFieldResidual(const MatrixFreeAPhiResiduals &residuals)
{
    return SMAX(SMAX(residuals.residual_ax_l2, residuals.residual_ay_l2),
                SMAX(residuals.residual_az_l2, residuals.residual_phi_l2));
}

inline ScalarComplexHelmholtzSolverState solveOneAComponent(const MatrixFreePairwiseGraph &graph,
                                                            StdVec<Complex> &component_field,
                                                            const StdVec<Real> &magnetic_reluctivity,
                                                            const StdVec<Complex> &reaction_coefficient,
                                                            const StdVec<Complex> &rhs,
                                                            const ScalarComplexHelmholtzSolverParameters &parameters,
                                                            MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const MatrixFreeAPhiScalarHelmholtzBackendKind backend =
        selectMatrixFreeAPhiScalarHelmholtzBackend(graph, component_field.size(), true, false, sph_native_context);
    return solveScalarLaplaceHelmholtzWithSelectedBackend(
        backend, graph, component_field, magnetic_reluctivity, rhs, reaction_coefficient, parameters, false,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromClearedGraph(graph, current_field, magnetic_reluctivity, current_residuals);
        },
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            sph_native_context->operatorAssembly().accumulateScalarLaplaceHelmholtzResiduals(current_field, magnetic_reluctivity,
                                                                                             current_residuals);
        });
}

inline ScalarComplexHelmholtzSolverState solvePhiComponent(const MatrixFreePairwiseGraph &graph,
                                                           StdVec<Complex> &phi,
                                                           const StdVec<Real> &electrical_conductivity,
                                                           const StdVec<Complex> &rhs,
                                                           const ScalarComplexHelmholtzSolverParameters &parameters,
                                                           MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const StdVec<Complex> zero_reaction = zeroComplexField(phi.size());
    const MatrixFreeAPhiScalarHelmholtzBackendKind backend =
        selectMatrixFreeAPhiScalarHelmholtzBackend(graph, phi.size(), true, false, sph_native_context);
    return solveScalarLaplaceHelmholtzWithSelectedBackend(
        backend, graph, phi, electrical_conductivity, rhs, zero_reaction, parameters, true,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromClearedGraph(graph, current_field, electrical_conductivity, current_residuals);
        },
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            sph_native_context->operatorAssembly().accumulateScalarLaplaceHelmholtzResiduals(current_field, electrical_conductivity,
                                                                                             current_residuals);
        });
}

inline ScalarComplexHelmholtzSolverState solveGaugeComponent(const MatrixFreePairwiseGraph &graph,
                                                             StdVec<Complex> &chi, const StdVec<Complex> &rhs,
                                                             const ScalarComplexHelmholtzSolverParameters &parameters,
                                                             bool use_operator_consistent_projection,
                                                             MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const StdVec<Real> unit_diffusion(chi.size(), Real(1.0));
    const MatrixFreeAPhiSolverOperatorSemantics semantics =
        selectMatrixFreeAPhiSolverOperatorSemantics(sph_native_context);
    const MatrixFreeAPhiScalarHelmholtzBackendKind backend = selectMatrixFreeAPhiScalarHelmholtzBackend(
        graph, chi.size(), !use_operator_consistent_projection, true, sph_native_context);

    if (!use_operator_consistent_projection)
    {
        const StdVec<Complex> zero_reaction = zeroComplexField(chi.size());
        ScalarComplexHelmholtzSolverParameters gauge_parameters = parameters;
        gauge_parameters.residual_growth_limit_factor_ = Real(5.0);
        gauge_parameters.remove_field_mean_after_jacobi_ = true;
        return solveScalarLaplaceHelmholtzWithSelectedBackend(
            backend, graph, chi, unit_diffusion, rhs, zero_reaction, gauge_parameters, true,
            [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
            {
                accumulateScalarLaplaceResidualsFromClearedGraph(graph, current_field, unit_diffusion, current_residuals);
            },
            [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
            {
                sph_native_context->operatorAssembly().accumulateScalarLaplaceHelmholtzResiduals(current_field, unit_diffusion,
                                                                                                 current_residuals);
            });
    }
    return solveScalarHelmholtzWithNormalizedResidualStep(
        chi, rhs, parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarDivergenceOfGradientHelmholtzResiduals(graph, current_field, current_residuals,
                                                                  sph_native_context, semantics.gauge);
        });
}
} // namespace

inline MatrixFreeAPhiScalarHelmholtzBackendKind selectMatrixFreeAPhiScalarHelmholtzBackend(
    const MatrixFreePairwiseGraph &graph, size_t field_size, bool prefer_sph_native, bool require_flat_graph_compatibility,
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (prefer_sph_native && useMatrixFreeAPhiSphNativeHelmholtz(sph_native_context))
    {
        if (matrixFreeAPhiLinearAlgebraBackend() == MatrixFreeAPhiLinearAlgebraBackendKind::ReservedDeviceKrylov)
        {
            return MatrixFreeAPhiScalarHelmholtzBackendKind::ReservedSphNativeSpecialized;
        }
        return MatrixFreeAPhiScalarHelmholtzBackendKind::SphNativeHostControlled;
    }
#if SPHINXSYS_USE_SYCL
    const bool graph_sycl_ready = matrixFreeLaplaceResidualUseSycl() && matrixFreeJacobiUseSycl();
    const bool flat_graph_ready =
        !require_flat_graph_compatibility || isFlatGraphFieldCompatible(graph, field_size);
    if (graph_sycl_ready && flat_graph_ready)
    {
        return MatrixFreeAPhiScalarHelmholtzBackendKind::GraphFusedSycl;
    }
#endif
    return MatrixFreeAPhiScalarHelmholtzBackendKind::GraphHostControlled;
}

#if SPHINXSYS_USE_SYCL
inline size_t matrixFreeFusedGraphHelmholtzSolveCount()
{
    return matrixFreeSyclGraphHelmholtzBuffers().solve_count_;
}

inline size_t matrixFreeFusedGraphHelmholtzMetricDownloadCount()
{
    return matrixFreeSyclGraphHelmholtzBuffers().metric_download_count_;
}

inline size_t matrixFreeFusedGraphHelmholtzMetricScalarDownloadCount()
{
    return matrixFreeSyclGraphHelmholtzBuffers().metric_scalar_download_count_;
}

inline size_t matrixFreeFusedGraphHelmholtzRhsUploadCount()
{
    return matrixFreeSyclGraphHelmholtzBuffers().rhs_upload_count_;
}

inline size_t matrixFreeFusedGraphHelmholtzReactionUploadCount()
{
    return matrixFreeSyclGraphHelmholtzBuffers().reaction_upload_count_;
}
#endif

inline MatrixFreeAPhiSolverParameters::MatrixFreeAPhiSolverParameters()
    : max_outer_iterations(12), residual_tolerance(1.0e-5), divergence_tolerance(1.0e-5),
      enable_gauge_projection(true), remove_gauge_mean_offset(true), use_operator_consistent_gauge_projection(true),
      enable_gauge_penalty(false), gauge_penalty_coefficient(0.0), gauge_penalty_ramp_iterations(0),
      gauge_penalty_initial_ratio(1.0), outer_relaxation_factor(1.0), update_tolerance(1.0e-7), residual_growth_guard_start_iteration(1), residual_growth_limit(5.0), source_ramp_iterations(0), source_initial_ratio(1.0)
{
    a_component_solver.max_iterations_ = 300;
    a_component_solver.relaxation_factor_ = 1.0;
    a_component_solver.absolute_tolerance_ = 1.0e-5;

    phi_solver.max_iterations_ = 300;
    phi_solver.relaxation_factor_ = 1.0;
    phi_solver.absolute_tolerance_ = 1.0e-5;

    gauge_solver.max_iterations_ = 500;
    gauge_solver.relaxation_factor_ = 1.0;
    gauge_solver.absolute_tolerance_ = 1.0e-5;
}

inline MatrixFreeAPhiGaugeStepDiagnostics::MatrixFreeAPhiGaugeStepDiagnostics()
    : applied(false), chi_l2(0.0), grad_chi_l2(0.0), chi_max_abs(0.0), phi_reference_offset_after_phi_solve_abs(0.0),
      phi_reference_offset_after_gauge_update_abs(0.0), phi_reference_offset_after_final_reference_abs(0.0),
      gauge_rhs_mean_abs_before_projection(0.0), gauge_rhs_mean_abs_after_projection(0.0), chi_mean_abs_after_solve(0.0),
      divergence_a_before_l2(0.0), divergence_a_after_raw_l2(0.0), divergence_a_after_final_l2(0.0),
      divergence_j_before_l2(0.0), divergence_j_after_raw_l2(0.0), divergence_j_after_final_l2(0.0),
      electric_field_before_l2(0.0), electric_field_after_raw_l2(0.0), electric_field_after_final_l2(0.0),
      electric_field_change_raw_l2(0.0), electric_field_change_final_l2(0.0), current_density_before_l2(0.0),
      current_density_after_raw_l2(0.0), current_density_after_final_l2(0.0), current_density_change_raw_l2(0.0),
      current_density_change_final_l2(0.0)
{
}

inline MatrixFreeAPhiGaugeProjectionResult::MatrixFreeAPhiGaugeProjectionResult() {}

inline void MatrixFreeAPhiWorkspace::resize(size_t size)
{
    previous_fields.ax.resize(size);
    previous_fields.ay.resize(size);
    previous_fields.az.resize(size);
    previous_fields.phi.resize(size);
    effective_sources.source_ax.resize(size);
    effective_sources.source_ay.resize(size);
    effective_sources.source_az.resize(size);
    effective_sources.source_phi.resize(size);
    sigma_grad_phi.resize(size, Vec3c::Zero());
    gauge_penalty_gradient.resize(size, Vec3c::Zero());
    rhs_ax.resize(size, Complex(0.0, 0.0));
    rhs_ay.resize(size, Complex(0.0, 0.0));
    rhs_az.resize(size, Complex(0.0, 0.0));
    divergence_sigma_a.resize(size, Complex(0.0, 0.0));
    rhs_phi.resize(size, Complex(0.0, 0.0));
}

inline void MatrixFreeAPhiStaggeredOperatorBundle::resize(size_t size)
{
    sigma_grad_phi.resize(size, Vec3c::Zero());
    gauge_penalty_gradient.resize(size, Vec3c::Zero());
}

inline MatrixFreeAPhiSolverState::MatrixFreeAPhiSolverState()
    : outer_iterations(0), converged(false), effective_gauge_penalty_coefficient(0.0),
      effective_source_scale(1.0), field_update_l2(0.0), relative_field_update_l2(0.0) {}

inline void enforcePhiReference(MatrixFreeAPhiFields &fields, const MatrixFreeAPhiParameters &parameters)
{
    if (!parameters.fix_phi_reference || fields.phi.empty() || parameters.phi_reference_index >= fields.phi.size())
    {
        return;
    }
    const Complex offset = fields.phi[parameters.phi_reference_index] - parameters.phi_reference_value;
    for (Complex &value : fields.phi)
    {
        value -= offset;
    }
#if SPHINXSYS_USE_SYCL
    invalidateMatrixFreeAPhiSyclValueFieldCache();
#endif
}

inline MatrixFreeAPhiGaugeProjectionResult applyMatrixFreeAPhiGaugeProjectionStep(
    const MatrixFreePairwiseGraph &graph, MatrixFreeAPhiFields &fields, const StdVec<Real> &electrical_conductivity,
    const MatrixFreeAPhiParameters &parameters, const ScalarComplexHelmholtzSolverParameters &gauge_solver_parameters,
    bool remove_gauge_mean_offset, bool enforce_phi_reference_after_gauge, bool use_operator_consistent_gauge_projection,
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    MatrixFreeAPhiGaugeProjectionResult result;
    result.diagnostics.applied = true;
    const MatrixFreeAPhiSolverOperatorSemantics semantics =
        selectMatrixFreeAPhiSolverOperatorSemantics(sph_native_context);

    const Complex imaginary_unit(0.0, 1.0);
    const MatrixFreeAPhiGaugeEvaluationBundle before_fields = evaluateMatrixFreeAPhiGaugeFields(
        graph, fields, electrical_conductivity, parameters, sph_native_context, semantics);

    result.diagnostics.divergence_a_before_l2 = computeComplexFieldL2Norm(before_fields.divergence_a);
    result.diagnostics.divergence_j_before_l2 = computeComplexFieldL2Norm(before_fields.divergence_j);
    result.diagnostics.electric_field_before_l2 = computeVectorComplexFieldL2Norm(before_fields.electric_field);
    result.diagnostics.current_density_before_l2 = computeVectorComplexFieldL2Norm(before_fields.current_density);

    StdVec<Complex> chi(fields.ax.size(), Complex(0.0, 0.0));
    StdVec<Complex> gauge_rhs = before_fields.divergence_a;
    result.diagnostics.gauge_rhs_mean_abs_before_projection = std::abs(computeComplexFieldMean(gauge_rhs));
    removeMeanOffset(gauge_rhs);
    result.diagnostics.gauge_rhs_mean_abs_after_projection = std::abs(computeComplexFieldMean(gauge_rhs));
    result.solver_state = solveGaugeComponent(graph, chi, gauge_rhs, gauge_solver_parameters,
                                              use_operator_consistent_gauge_projection, sph_native_context);
    if (remove_gauge_mean_offset)
    {
        removeMeanOffset(chi);
    }
    result.diagnostics.chi_mean_abs_after_solve = std::abs(computeComplexFieldMean(chi));
    result.diagnostics.chi_l2 = computeComplexFieldL2Norm(chi);
    result.diagnostics.chi_max_abs = computeComplexFieldMaxAbs(chi);

    const StdVec<Vec3c> grad_chi = applyMatrixFreeScalarGradient(graph, chi, sph_native_context, semantics.gauge);
    result.diagnostics.grad_chi_l2 = computeVectorComplexFieldL2Norm(grad_chi);
    for (size_t i = 0; i != fields.ax.size(); ++i)
    {
        fields.ax[i] -= grad_chi[i][0];
        fields.ay[i] -= grad_chi[i][1];
        fields.az[i] -= grad_chi[i][2];
        fields.phi[i] += imaginary_unit * parameters.angular_frequency * chi[i];
    }
#if SPHINXSYS_USE_SYCL
    invalidateMatrixFreeAPhiSyclValueFieldCache();
#endif
    result.diagnostics.phi_reference_offset_after_gauge_update_abs =
        std::abs(computePhiReferenceOffset(fields, parameters));

    const MatrixFreeAPhiGaugeEvaluationBundle raw_after_fields = evaluateMatrixFreeAPhiGaugeFields(
        graph, fields, electrical_conductivity, parameters, sph_native_context, semantics);

    result.diagnostics.divergence_a_after_raw_l2 = computeComplexFieldL2Norm(raw_after_fields.divergence_a);
    result.diagnostics.divergence_j_after_raw_l2 = computeComplexFieldL2Norm(raw_after_fields.divergence_j);
    result.diagnostics.electric_field_after_raw_l2 = computeVectorComplexFieldL2Norm(raw_after_fields.electric_field);
    result.diagnostics.current_density_after_raw_l2 = computeVectorComplexFieldL2Norm(raw_after_fields.current_density);
    result.diagnostics.electric_field_change_raw_l2 =
        computeVectorComplexFieldL2Difference(raw_after_fields.electric_field, before_fields.electric_field);
    result.diagnostics.current_density_change_raw_l2 =
        computeVectorComplexFieldL2Difference(raw_after_fields.current_density, before_fields.current_density);

    if (enforce_phi_reference_after_gauge)
    {
        enforcePhiReference(fields, parameters);
    }
    result.diagnostics.phi_reference_offset_after_final_reference_abs =
        std::abs(computePhiReferenceOffset(fields, parameters));

    const MatrixFreeAPhiGaugeEvaluationBundle final_fields = evaluateMatrixFreeAPhiGaugeFields(
        graph, fields, electrical_conductivity, parameters, sph_native_context, semantics);

    result.diagnostics.divergence_a_after_final_l2 = computeComplexFieldL2Norm(final_fields.divergence_a);
    result.diagnostics.divergence_j_after_final_l2 = computeComplexFieldL2Norm(final_fields.divergence_j);
    result.diagnostics.electric_field_after_final_l2 = computeVectorComplexFieldL2Norm(final_fields.electric_field);
    result.diagnostics.current_density_after_final_l2 = computeVectorComplexFieldL2Norm(final_fields.current_density);
    result.diagnostics.electric_field_change_final_l2 =
        computeVectorComplexFieldL2Difference(final_fields.electric_field, before_fields.electric_field);
    result.diagnostics.current_density_change_final_l2 =
        computeVectorComplexFieldL2Difference(final_fields.current_density, before_fields.current_density);

    return result;
}

inline MatrixFreeAPhiSolverState solveMatrixFreeAPhiStaggered(const MatrixFreePairwiseGraph &graph,
                                                              MatrixFreeAPhiFields &fields,
                                                              const StdVec<Real> &electrical_conductivity,
                                                              const StdVec<Real> &magnetic_reluctivity,
                                                              const MatrixFreeAPhiParameters &parameters,
                                                              const MatrixFreeAPhiSources &sources,
                                                              const MatrixFreeAPhiSolverParameters &solver_parameters,
                                                              MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    MatrixFreeAPhiSolverState state;
    const MatrixFreeAPhiSolverOperatorSemantics semantics =
        selectMatrixFreeAPhiSolverOperatorSemantics(sph_native_context);
    const size_t number_of_particles = fields.ax.size();
    const StdVec<Complex> reaction_coefficient =
        buildReactionCoefficient(electrical_conductivity, parameters.angular_frequency);
    MatrixFreeAPhiWorkspace workspace;
    workspace.resize(number_of_particles);
#if SPHINXSYS_USE_SYCL
    matrixFreeAPhiSyclBindTimestepResources(graph, number_of_particles);
#endif

    MatrixFreeAPhiFields last_stable_fields = fields;
    MatrixFreeAPhiResiduals last_stable_residuals;
    last_stable_residuals.resize(number_of_particles);
    last_stable_residuals.clear();
    Real last_stable_max_field_residual = MaxReal;
    const Real residual_growth_limit = solver_parameters.residual_growth_limit;

    for (size_t outer_iteration = 0; outer_iteration != solver_parameters.max_outer_iterations; ++outer_iteration)
    {
        workspace.previous_fields = fields;
        const Real effective_penalty = computeMatrixFreeAPhiEffectiveGaugePenalty(solver_parameters, outer_iteration);
        state.effective_gauge_penalty_coefficient = effective_penalty;

        const Real effective_source_scale = computeMatrixFreeAPhiEffectiveSourceScale(solver_parameters, outer_iteration);
        state.effective_source_scale = effective_source_scale;

        MatrixFreeAPhiSources &effective_sources = workspace.effective_sources;
        scaleMatrixFreeAPhiSources(sources, effective_source_scale, effective_sources);

        const MatrixFreeAPhiStaggeredOperatorBundle operators = assembleMatrixFreeAPhiStaggeredOperators(
            graph, fields, electrical_conductivity, solver_parameters.enable_gauge_penalty, semantics,
            sph_native_context);
        workspace.sigma_grad_phi = operators.sigma_grad_phi;
        workspace.gauge_penalty_gradient = operators.gauge_penalty_gradient;
        assembleMatrixFreeAPhiARhs(effective_sources, operators, effective_penalty, workspace);

        state.ax_state = solveOneAComponent(graph, fields.ax, magnetic_reluctivity, reaction_coefficient, workspace.rhs_ax,
                                            solver_parameters.a_component_solver, sph_native_context);
        state.ay_state = solveOneAComponent(graph, fields.ay, magnetic_reluctivity, reaction_coefficient, workspace.rhs_ay,
                                            solver_parameters.a_component_solver, sph_native_context);
        state.az_state = solveOneAComponent(graph, fields.az, magnetic_reluctivity, reaction_coefficient, workspace.rhs_az,
                                            solver_parameters.a_component_solver, sph_native_context);

        workspace.divergence_sigma_a =
            computeDivergenceOfSigmaA(graph, fields, electrical_conductivity, sph_native_context,
                                      semantics.particle_kernel);
        assembleMatrixFreeAPhiPhiRhs(effective_sources, workspace.divergence_sigma_a, parameters, workspace.rhs_phi);
        state.phi_state =
            solvePhiComponent(graph, fields.phi, electrical_conductivity, workspace.rhs_phi, solver_parameters.phi_solver,
                              sph_native_context);
        const Real phi_reference_offset_after_phi_solve_abs =
            std::abs(computePhiReferenceOffset(fields, parameters));
        state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs = phi_reference_offset_after_phi_solve_abs;
        enforcePhiReference(fields, parameters);

        if (solver_parameters.enable_gauge_projection)
        {
            const MatrixFreeAPhiGaugeProjectionResult gauge_result = applyMatrixFreeAPhiGaugeProjectionStep(
                graph, fields, electrical_conductivity, parameters, solver_parameters.gauge_solver,
                solver_parameters.remove_gauge_mean_offset, true, solver_parameters.use_operator_consistent_gauge_projection,
                sph_native_context);
            state.gauge_state = gauge_result.solver_state;
            state.gauge_diagnostics = gauge_result.diagnostics;
            state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs = phi_reference_offset_after_phi_solve_abs;
        }

        blendAPhiFields(workspace.previous_fields, fields, solver_parameters.outer_relaxation_factor);
        enforcePhiReference(fields, parameters);
        finalizeMatrixFreeAPhiFieldUpdateMetrics(workspace.previous_fields, fields, state);

        state.residuals = evaluateMatrixFreeAPhiResiduals(graph, fields, electrical_conductivity, magnetic_reluctivity,
                                                          parameters, effective_sources,
                                                          solver_parameters.enable_gauge_penalty,
                                                          state.effective_gauge_penalty_coefficient, sph_native_context);
        state.outer_iterations = outer_iteration + 1;
        appendMatrixFreeAPhiIterationRecord(state);

        const Real max_field_residual = computeMatrixFreeAPhiMaxFieldResidual(state.residuals);

        if (!std::isfinite(max_field_residual) || !std::isfinite(state.residuals.divergence_a_l2) ||
            !std::isfinite(state.residuals.divergence_j_l2))
        {
            fields = last_stable_fields;
            state.residuals = last_stable_residuals;
            break;
        }

        const bool growth_guard_active =
            state.outer_iterations > solver_parameters.residual_growth_guard_start_iteration;
        if (growth_guard_active && last_stable_max_field_residual < MaxReal &&
            max_field_residual > residual_growth_limit * (last_stable_max_field_residual + TinyReal))
        {
            fields = last_stable_fields;
            state.residuals = last_stable_residuals;
            break;
        }

        last_stable_fields = fields;
        last_stable_residuals = state.residuals;
        last_stable_max_field_residual = max_field_residual;

        const bool divergence_ok = !solver_parameters.enable_gauge_projection ||
                                   state.residuals.divergence_a_l2 <= solver_parameters.divergence_tolerance;
        const bool residual_ok = max_field_residual <= solver_parameters.residual_tolerance;
        const bool update_ok = state.relative_field_update_l2 <= solver_parameters.update_tolerance;
        if (divergence_ok && (residual_ok || update_ok))
        {
            state.converged = true;
            break;
        }
    }

    return state;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_SOLVER_HPP
