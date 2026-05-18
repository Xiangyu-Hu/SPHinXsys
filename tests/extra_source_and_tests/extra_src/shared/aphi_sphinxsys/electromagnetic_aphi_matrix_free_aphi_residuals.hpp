#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline MatrixFreeAPhiResiduals::MatrixFreeAPhiResiduals()
    : residual_ax_l2(0.0), residual_ay_l2(0.0), residual_az_l2(0.0), residual_phi_l2(0.0), divergence_a_l2(0.0),
      divergence_j_l2(0.0), interface_edge_count(0), high_contrast_edge_count(0), high_contrast_particle_count(0),
      max_interface_sigma_ratio(1.0), high_contrast_sigma_grad_phi_l2(0.0), high_contrast_div_sigma_a_l2(0.0),
      high_contrast_residual_a_l2(0.0), high_contrast_residual_phi_l2(0.0)
{
}

inline void MatrixFreeAPhiResiduals::resize(size_t size)
{
    residual_ax.resize(size, Complex(0.0, 0.0));
    residual_ay.resize(size, Complex(0.0, 0.0));
    residual_az.resize(size, Complex(0.0, 0.0));
    residual_phi.resize(size, Complex(0.0, 0.0));
    divergence_a.resize(size, Complex(0.0, 0.0));
    divergence_j.resize(size, Complex(0.0, 0.0));
    clear();
}

inline void MatrixFreeAPhiResiduals::clear()
{
    std::fill(residual_ax.begin(), residual_ax.end(), Complex(0.0, 0.0));
    std::fill(residual_ay.begin(), residual_ay.end(), Complex(0.0, 0.0));
    std::fill(residual_az.begin(), residual_az.end(), Complex(0.0, 0.0));
    std::fill(residual_phi.begin(), residual_phi.end(), Complex(0.0, 0.0));
    std::fill(divergence_a.begin(), divergence_a.end(), Complex(0.0, 0.0));
    std::fill(divergence_j.begin(), divergence_j.end(), Complex(0.0, 0.0));
    residual_ax_l2 = 0.0;
    residual_ay_l2 = 0.0;
    residual_az_l2 = 0.0;
    residual_phi_l2 = 0.0;
    divergence_a_l2 = 0.0;
    divergence_j_l2 = 0.0;
    interface_edge_count = 0;
    high_contrast_edge_count = 0;
    high_contrast_particle_count = 0;
    max_interface_sigma_ratio = 1.0;
    high_contrast_sigma_grad_phi_l2 = 0.0;
    high_contrast_div_sigma_a_l2 = 0.0;
    high_contrast_residual_a_l2 = 0.0;
    high_contrast_residual_phi_l2 = 0.0;
}

inline Real computeScalarFieldL2(const StdVec<Complex> &field)
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

inline StdVec<Complex> computeHarmonicDivergenceOfVectorField(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
    const StdVec<Complex> &field_z, const StdVec<Real> &edge_weight_coefficient)
{
    const size_t number_of_particles = field_x.size();
    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    if (field_y.size() != number_of_particles || field_z.size() != number_of_particles ||
        edge_weight_coefficient.size() != number_of_particles)
    {
        return divergence;
    }

#if SPHINXSYS_USE_SYCL
    if (matrixFreeHarmonicGradientUseSycl() && graph.flat_.isValidFor(graph.rows_.size()) &&
        graph.rows_.size() == number_of_particles)
    {
        StdVec<Real> field_x_real(number_of_particles, 0.0);
        StdVec<Real> field_x_imag(number_of_particles, 0.0);
        StdVec<Real> field_y_real(number_of_particles, 0.0);
        StdVec<Real> field_y_imag(number_of_particles, 0.0);
        StdVec<Real> field_z_real(number_of_particles, 0.0);
        StdVec<Real> field_z_imag(number_of_particles, 0.0);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            field_x_real[i] = field_x[i].real();
            field_x_imag[i] = field_x[i].imag();
            field_y_real[i] = field_y[i].real();
            field_y_imag[i] = field_y[i].imag();
            field_z_real[i] = field_z[i].real();
            field_z_imag[i] = field_z[i].imag();
        }

        MatrixFreeAPhiSyclWorkspace &sycl_workspace = matrixFreeAPhiSyclWorkspace();
        if (tryMatrixFreeHarmonicScalarDivergenceThreeComponentsWithStagingSycl(
                graph, sycl_workspace, field_x_real, field_x_imag, field_y_real, field_y_imag, field_z_real, field_z_imag,
                edge_weight_coefficient, divergence))
        {
            return divergence;
        }

        MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
        graph_buffers.loadIfNeeded(graph.flat_, number_of_particles);
        MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
        value_buffers.resize(number_of_particles);
        value_buffers.loadCoefficientIfNeeded(edge_weight_coefficient);

        sycl::buffer<Real, 1> field_x_real_buffer(field_x_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_x_imag_buffer(field_x_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_y_real_buffer(field_y_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_y_imag_buffer(field_y_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_z_real_buffer(field_z_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_z_imag_buffer(field_z_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> divergence_real_buffer{sycl::range<1>(number_of_particles)};
        sycl::buffer<Real, 1> divergence_imag_buffer{sycl::range<1>(number_of_particles)};

        {
            if (graph_buffers.deviceTopologyUsmReady())
            {
                const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
                const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
                const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
                const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
                const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto fxr = field_x_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fxi = field_x_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyr = field_y_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyi = field_y_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzr = field_z_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzi = field_z_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto divr = divergence_real_buffer.get_access<sycl::access::mode::write>(cgh);
                    auto divi = divergence_imag_buffer.get_access<sycl::access::mode::write>(cgh);

                    cgh.parallel_for(sycl::range<1>(number_of_particles), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        Real div_real = 0.0;
                        Real div_imag = 0.0;
                        for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                        {
                            const size_t index_j = col_dev[edge];
                            const Real coefficient_ij =
                                2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                            const Real weighted_x = coefficient_ij * gwx_dev[edge];
                            const Real weighted_y = coefficient_ij * gwy_dev[edge];
                            const Real weighted_z = coefficient_ij * gwz_dev[edge];
                            div_real += weighted_x * (fxr[index_j] - fxr[index_i]) +
                                        weighted_y * (fyr[index_j] - fyr[index_i]) +
                                        weighted_z * (fzr[index_j] - fzr[index_i]);
                            div_imag += weighted_x * (fxi[index_j] - fxi[index_i]) +
                                        weighted_y * (fyi[index_j] - fyi[index_i]) +
                                        weighted_z * (fzi[index_j] - fzi[index_i]);
                        }
                        divr[index_i] = div_real;
                        divi[index_i] = div_imag;
                    });
                });
            }
            else
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto fxr = field_x_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fxi = field_x_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyr = field_y_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyi = field_y_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzr = field_z_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzi = field_z_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto divr = divergence_real_buffer.get_access<sycl::access::mode::write>(cgh);
                    auto divi = divergence_imag_buffer.get_access<sycl::access::mode::write>(cgh);

                    cgh.parallel_for(sycl::range<1>(number_of_particles), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        Real div_real = 0.0;
                        Real div_imag = 0.0;
                        for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                        {
                            const size_t index_j = column_indices[edge];
                            const Real coefficient_ij =
                                2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                            const Real weighted_x = coefficient_ij * gwx[edge];
                            const Real weighted_y = coefficient_ij * gwy[edge];
                            const Real weighted_z = coefficient_ij * gwz[edge];
                            div_real += weighted_x * (fxr[index_j] - fxr[index_i]) +
                                        weighted_y * (fyr[index_j] - fyr[index_i]) +
                                        weighted_z * (fzr[index_j] - fzr[index_i]);
                            div_imag += weighted_x * (fxi[index_j] - fxi[index_i]) +
                                        weighted_y * (fyi[index_j] - fyi[index_i]) +
                                        weighted_z * (fzi[index_j] - fzi[index_i]);
                        }
                        divr[index_i] = div_real;
                        divi[index_i] = div_imag;
                    });
                });
            }
        }

        sycl::host_accessor<Real, 1, sycl::access::mode::read> div_real(divergence_real_buffer);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> div_imag(divergence_imag_buffer);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            divergence[i] = Complex(div_real[i], div_imag[i]);
        }
        return divergence;
    }
#endif

    const StdVec<Vec3c> grad_x = applyMatrixFreeHarmonicWeightedGradient(graph, field_x, edge_weight_coefficient);
    const StdVec<Vec3c> grad_y = applyMatrixFreeHarmonicWeightedGradient(graph, field_y, edge_weight_coefficient);
    const StdVec<Vec3c> grad_z = applyMatrixFreeHarmonicWeightedGradient(graph, field_z, edge_weight_coefficient);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = grad_x[i][0] + grad_y[i][1] + grad_z[i][2];
    }
    return divergence;
}

/** Scalar divergence div(F) = d(Fx)/dx + d(Fy)/dy + d(Fz)/dz using the same pairwise gradient weights as
 *  `applyMatrixFreeGradient` (no harmonic edge weighting). SYCL path avoids three full Vec3 gradient downloads. */
inline StdVec<Complex> computeGradientDivergenceOfVectorField(const MatrixFreePairwiseGraph &graph,
                                                             const StdVec<Complex> &field_x,
                                                             const StdVec<Complex> &field_y,
                                                             const StdVec<Complex> &field_z)
{
    const size_t number_of_particles = field_x.size();
    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    if (field_y.size() != number_of_particles || field_z.size() != number_of_particles)
    {
        return divergence;
    }

#if SPHINXSYS_USE_SYCL
    if (matrixFreeGradientUseSycl() && graph.flat_.isValidFor(graph.rows_.size()) &&
        graph.rows_.size() == number_of_particles)
    {
        StdVec<Real> field_x_real(number_of_particles, 0.0);
        StdVec<Real> field_x_imag(number_of_particles, 0.0);
        StdVec<Real> field_y_real(number_of_particles, 0.0);
        StdVec<Real> field_y_imag(number_of_particles, 0.0);
        StdVec<Real> field_z_real(number_of_particles, 0.0);
        StdVec<Real> field_z_imag(number_of_particles, 0.0);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            field_x_real[i] = field_x[i].real();
            field_x_imag[i] = field_x[i].imag();
            field_y_real[i] = field_y[i].real();
            field_y_imag[i] = field_y[i].imag();
            field_z_real[i] = field_z[i].real();
            field_z_imag[i] = field_z[i].imag();
        }

        MatrixFreeAPhiSyclWorkspace &sycl_workspace = matrixFreeAPhiSyclWorkspace();
        if (tryMatrixFreeGradientScalarDivergenceThreeComponentsWithStagingSycl(
                graph, sycl_workspace, field_x_real, field_x_imag, field_y_real, field_y_imag, field_z_real, field_z_imag,
                divergence))
        {
            return divergence;
        }

        MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
        graph_buffers.loadIfNeeded(graph.flat_, number_of_particles);

        sycl::buffer<Real, 1> field_x_real_buffer(field_x_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_x_imag_buffer(field_x_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_y_real_buffer(field_y_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_y_imag_buffer(field_y_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_z_real_buffer(field_z_real.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> field_z_imag_buffer(field_z_imag.data(), sycl::range<1>(number_of_particles));
        sycl::buffer<Real, 1> divergence_real_buffer{sycl::range<1>(number_of_particles)};
        sycl::buffer<Real, 1> divergence_imag_buffer{sycl::range<1>(number_of_particles)};

        {
            if (graph_buffers.deviceTopologyUsmReady())
            {
                const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
                const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
                const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
                const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
                const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto fxr = field_x_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fxi = field_x_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyr = field_y_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyi = field_y_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzr = field_z_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzi = field_z_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto divr = divergence_real_buffer.get_access<sycl::access::mode::write>(cgh);
                    auto divi = divergence_imag_buffer.get_access<sycl::access::mode::write>(cgh);

                    cgh.parallel_for(sycl::range<1>(number_of_particles), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        Real div_real = 0.0;
                        Real div_imag = 0.0;
                        for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                        {
                            const size_t index_j = col_dev[edge];
                            div_real += gwx_dev[edge] * (fxr[index_j] - fxr[index_i]) +
                                        gwy_dev[edge] * (fyr[index_j] - fyr[index_i]) +
                                        gwz_dev[edge] * (fzr[index_j] - fzr[index_i]);
                            div_imag += gwx_dev[edge] * (fxi[index_j] - fxi[index_i]) +
                                        gwy_dev[edge] * (fyi[index_j] - fyi[index_i]) +
                                        gwz_dev[edge] * (fzi[index_j] - fzi[index_i]);
                        }
                        divr[index_i] = div_real;
                        divi[index_i] = div_imag;
                    });
                });
            }
            else
            {
                matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                    auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                    auto fxr = field_x_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fxi = field_x_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyr = field_y_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fyi = field_y_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzr = field_z_real_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto fzi = field_z_imag_buffer.get_access<sycl::access::mode::read>(cgh);
                    auto divr = divergence_real_buffer.get_access<sycl::access::mode::write>(cgh);
                    auto divi = divergence_imag_buffer.get_access<sycl::access::mode::write>(cgh);

                    cgh.parallel_for(sycl::range<1>(number_of_particles), [=](sycl::id<1> id) {
                        const size_t index_i = id[0];
                        Real div_real = 0.0;
                        Real div_imag = 0.0;
                        for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                        {
                            const size_t index_j = column_indices[edge];
                            div_real += gwx[edge] * (fxr[index_j] - fxr[index_i]) + gwy[edge] * (fyr[index_j] - fyr[index_i]) +
                                        gwz[edge] * (fzr[index_j] - fzr[index_i]);
                            div_imag += gwx[edge] * (fxi[index_j] - fxi[index_i]) + gwy[edge] * (fyi[index_j] - fyi[index_i]) +
                                        gwz[edge] * (fzi[index_j] - fzi[index_i]);
                        }
                        divr[index_i] = div_real;
                        divi[index_i] = div_imag;
                    });
                });
            }
        }

        sycl::host_accessor<Real, 1, sycl::access::mode::read> div_real(divergence_real_buffer);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> div_imag(divergence_imag_buffer);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            divergence[i] = Complex(div_real[i], div_imag[i]);
        }
        return divergence;
    }
#endif

    const StdVec<Vec3c> grad_x = applyMatrixFreeGradient(graph, field_x);
    const StdVec<Vec3c> grad_y = applyMatrixFreeGradient(graph, field_y);
    const StdVec<Vec3c> grad_z = applyMatrixFreeGradient(graph, field_z);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = grad_x[i][0] + grad_y[i][1] + grad_z[i][2];
    }
    return divergence;
}

inline StdVec<Complex> computeResidualDivergenceOfCurrentDensity(
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

inline StdVec<Vec3c> computeResidualScalarGradient(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
    MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeStandardGradient(field);
    }
    return applyMatrixFreeGradient(graph, field);
}

inline StdVec<Complex> computeResidualDivergenceOfA(const MatrixFreePairwiseGraph &graph,
                                                    const MatrixFreeAPhiFields &fields,
                                                    MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
                                                    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeGradientDivergenceOfVectorField(fields.ax, fields.ay, fields.az);
    }
    return computeGradientDivergenceOfVectorField(graph, fields.ax, fields.ay, fields.az);
}

inline StdVec<Vec3c> computeResidualGaugePenaltyGradient(const MatrixFreePairwiseGraph &graph,
                                                         const StdVec<Complex> &divergence_a_field,
                                                         MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
                                                         MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    return computeResidualScalarGradient(graph, divergence_a_field, semantics, sph_native_context);
}

inline StdVec<Complex> computeResidualNegativeLaplace(const MatrixFreePairwiseGraph &graph,
                                                      const StdVec<Complex> &field,
                                                      const StdVec<Real> &diffusion_coefficient,
                                                      MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
                                                      MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().laplaceOperator().apply(field, diffusion_coefficient);
    }
    return applyScalarNegativeLaplaceFromGraph(graph, field, diffusion_coefficient);
}

inline StdVec<Vec3c> computeResidualSigmaGradPhi(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &phi,
                                                 const StdVec<Real> &electrical_conductivity,
                                                 MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
                                                 MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeSigmaGradPhi(phi, electrical_conductivity);
    }
    return applyMatrixFreeHarmonicWeightedGradient(graph, phi, electrical_conductivity);
}

inline StdVec<Complex> computeResidualDivergenceOfSigmaA(const MatrixFreePairwiseGraph &graph,
                                                         const MatrixFreeAPhiFields &fields,
                                                         const StdVec<Real> &electrical_conductivity,
                                                         MatrixFreeAPhiDiscreteOperatorSemanticsKind semantics,
                                                         MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    if (semantics == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel)
    {
        return sph_native_context->operatorAssembly().computeHarmonicDivergenceOfVectorField(
            fields.ax, fields.ay, fields.az, electrical_conductivity);
    }
    return computeHarmonicDivergenceOfVectorField(graph, fields.ax, fields.ay, fields.az, electrical_conductivity);
}

inline MatrixFreeAPhiResidualOperatorSemantics selectMatrixFreeAPhiResidualOperatorSemantics(
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const MatrixFreeAPhiDiscreteOperatorSemanticsKind particle_kernel =
        useMatrixFreeAPhiSphNativePath(sph_native_context)
            ? MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel
            : MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility;
    return {particle_kernel, MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility};
}

struct MatrixFreeAPhiResidualOperatorBundle
{
    StdVec<Complex> laplace_ax;
    StdVec<Complex> laplace_ay;
    StdVec<Complex> laplace_az;
    StdVec<Complex> laplace_phi;
    StdVec<Vec3c> grad_phi;
    StdVec<Vec3c> sigma_grad_phi;
    StdVec<Complex> divergence_a_field;
    StdVec<Vec3c> gauge_penalty_gradient;
    StdVec<Complex> divergence_sigma_a;

    void resize(size_t size)
    {
        grad_phi.assign(size, Vec3c::Zero());
        sigma_grad_phi.assign(size, Vec3c::Zero());
        divergence_a_field.assign(size, Complex(0.0, 0.0));
        gauge_penalty_gradient.assign(size, Vec3c::Zero());
    }
};

inline MatrixFreeAPhiResidualOperatorBundle assembleMatrixFreeAPhiResidualOperators(
    const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiFields &fields,
    const StdVec<Real> &electrical_conductivity, const StdVec<Real> &magnetic_reluctivity,
    bool enable_gauge_penalty, const MatrixFreeAPhiResidualOperatorSemantics &semantics,
    MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    MatrixFreeAPhiResidualOperatorBundle operators;
    operators.resize(fields.ax.size());

    operators.laplace_ax = computeResidualNegativeLaplace(graph, fields.ax, magnetic_reluctivity,
                                                          semantics.particle_kernel, sph_native_context);
    operators.laplace_ay = computeResidualNegativeLaplace(graph, fields.ay, magnetic_reluctivity,
                                                          semantics.particle_kernel, sph_native_context);
    operators.laplace_az = computeResidualNegativeLaplace(graph, fields.az, magnetic_reluctivity,
                                                          semantics.particle_kernel, sph_native_context);

    if (semantics.usesParticleKernel())
    {
        operators.grad_phi = computeResidualScalarGradient(graph, fields.phi, semantics.legacy_solver_gradient,
                                                           sph_native_context);
        operators.sigma_grad_phi = computeResidualSigmaGradPhi(graph, fields.phi, electrical_conductivity,
                                                               semantics.particle_kernel, sph_native_context);
        operators.divergence_a_field = computeResidualDivergenceOfA(graph, fields, semantics.legacy_solver_gradient,
                                                                    sph_native_context);
        if (enable_gauge_penalty)
        {
            operators.gauge_penalty_gradient = computeResidualGaugePenaltyGradient(
                graph, operators.divergence_a_field, semantics.legacy_solver_gradient, sph_native_context);
        }
    }
    else
    {
#if SPHINXSYS_USE_SYCL
        const bool fused_phi_gradients_ok = tryApplyPhiStandardAndHarmonicWeightedGradientsFusedSycl(
            graph, fields.phi, electrical_conductivity, operators.grad_phi, operators.sigma_grad_phi);
#else
        const bool fused_phi_gradients_ok = false;
#endif
        if (!fused_phi_gradients_ok)
        {
            operators.grad_phi = computeResidualScalarGradient(graph, fields.phi, semantics.legacy_solver_gradient,
                                                               sph_native_context);
            operators.sigma_grad_phi = computeResidualSigmaGradPhi(graph, fields.phi, electrical_conductivity,
                                                                   semantics.particle_kernel, sph_native_context);
        }
#if SPHINXSYS_USE_SYCL
        const bool fused_div_a_and_gauge_grad =
            enable_gauge_penalty &&
            semantics.legacy_solver_gradient ==
                MatrixFreeAPhiDiscreteOperatorSemanticsKind::LegacyGraphCompatibility &&
            matrixFreeGradientUseSycl() &&
            tryMatrixFreeDivergenceAndGradientOfDivergenceOfVectorPotentialSycl(
                graph, fields.ax, fields.ay, fields.az, &operators.divergence_a_field, operators.gauge_penalty_gradient);
#else
        const bool fused_div_a_and_gauge_grad = false;
#endif
        if (!fused_div_a_and_gauge_grad)
        {
            operators.divergence_a_field = computeResidualDivergenceOfA(graph, fields, semantics.legacy_solver_gradient,
                                                                        sph_native_context);
            if (enable_gauge_penalty)
            {
                operators.gauge_penalty_gradient = computeResidualGaugePenaltyGradient(
                    graph, operators.divergence_a_field, semantics.legacy_solver_gradient, sph_native_context);
            }
        }
    }

    operators.laplace_phi = computeResidualNegativeLaplace(graph, fields.phi, electrical_conductivity,
                                                           semantics.particle_kernel, sph_native_context);
    operators.divergence_sigma_a = computeResidualDivergenceOfSigmaA(graph, fields, electrical_conductivity,
                                                                     semantics.particle_kernel, sph_native_context);
    return operators;
}

struct MatrixFreeAPhiCurrentComponentFields
{
    StdVec<Complex> x;
    StdVec<Complex> y;
    StdVec<Complex> z;

    void resize(size_t size)
    {
        x.assign(size, Complex(0.0, 0.0));
        y.assign(size, Complex(0.0, 0.0));
        z.assign(size, Complex(0.0, 0.0));
    }
};

struct MatrixFreeAPhiInterfaceContrastRegion
{
    StdVec<uint8_t> high_contrast_particle_mask;
    size_t interface_edge_count = 0;
    size_t high_contrast_edge_count = 0;
    size_t high_contrast_particle_count = 0;
    Real max_interface_sigma_ratio = 1.0;
};

inline MatrixFreeAPhiInterfaceContrastRegion detectMatrixFreeAPhiInterfaceContrastRegion(
    const MatrixFreePairwiseGraph &graph, const StdVec<Real> &electrical_conductivity,
    Real interface_contrast_threshold)
{
    MatrixFreeAPhiInterfaceContrastRegion region;
    region.high_contrast_particle_mask.assign(electrical_conductivity.size(), uint8_t(0));

    for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
    {
        for (const MatrixFreePairwiseNeighborEntry &neighbor : graph.rows_[index_i].neighbors_)
        {
            const size_t index_j = neighbor.index_j_;
            if (index_i >= index_j)
            {
                continue;
            }

            const Real sigma_i = std::abs(electrical_conductivity[index_i]);
            const Real sigma_j = std::abs(electrical_conductivity[index_j]);
            const Real sigma_min = SMAX(SMIN(sigma_i, sigma_j), TinyReal);
            const Real sigma_max = SMAX(sigma_i, sigma_j);
            const Real sigma_ratio = sigma_max / sigma_min;
            region.max_interface_sigma_ratio = SMAX(region.max_interface_sigma_ratio, sigma_ratio);

            if (sigma_ratio > Real(1.0) + TinyReal)
            {
                ++region.interface_edge_count;
            }
            if (sigma_ratio + TinyReal >= interface_contrast_threshold)
            {
                ++region.high_contrast_edge_count;
                region.high_contrast_particle_mask[index_i] = uint8_t(1);
                region.high_contrast_particle_mask[index_j] = uint8_t(1);
            }
        }
    }

    for (size_t i = 0; i != region.high_contrast_particle_mask.size(); ++i)
    {
        region.high_contrast_particle_count +=
            region.high_contrast_particle_mask[i] != 0 ? size_t(1) : size_t(0);
    }

    return region;
}

inline void assembleMatrixFreeAPhiResidualFieldContributions(
    const MatrixFreeAPhiFields &fields, const StdVec<Real> &electrical_conductivity,
    const MatrixFreeAPhiParameters &parameters, const MatrixFreeAPhiSources &sources,
    const MatrixFreeAPhiResidualOperatorBundle &operators, MatrixFreeAPhiResiduals &residuals,
    MatrixFreeAPhiCurrentComponentFields &current_components)
{
    const Complex imaginary_unit(0.0, 1.0);
    const size_t number_of_particles = fields.ax.size();
    current_components.resize(number_of_particles);

    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.divergence_a[i] = operators.divergence_a_field[i];

        residuals.residual_ax[i] = operators.laplace_ax[i] + imaginary_unit * parameters.angular_frequency *
                                                      electrical_conductivity[i] * fields.ax[i] +
                                   operators.sigma_grad_phi[i][0] - sources.source_ax[i];
        residuals.residual_ay[i] = operators.laplace_ay[i] + imaginary_unit * parameters.angular_frequency *
                                                      electrical_conductivity[i] * fields.ay[i] +
                                   operators.sigma_grad_phi[i][1] - sources.source_ay[i];
        residuals.residual_az[i] = operators.laplace_az[i] + imaginary_unit * parameters.angular_frequency *
                                                      electrical_conductivity[i] * fields.az[i] +
                                   operators.sigma_grad_phi[i][2] - sources.source_az[i];

        residuals.residual_phi[i] =
            operators.laplace_phi[i] - imaginary_unit * parameters.angular_frequency * operators.divergence_sigma_a[i] -
            sources.source_phi[i];

        current_components.x[i] = -electrical_conductivity[i] *
                                  (imaginary_unit * parameters.angular_frequency * fields.ax[i] + operators.grad_phi[i][0]);
        current_components.y[i] = -electrical_conductivity[i] *
                                  (imaginary_unit * parameters.angular_frequency * fields.ay[i] + operators.grad_phi[i][1]);
        current_components.z[i] = -electrical_conductivity[i] *
                                  (imaginary_unit * parameters.angular_frequency * fields.az[i] + operators.grad_phi[i][2]);
    }
}

inline void applyMatrixFreeAPhiGaugePenaltyToResiduals(const MatrixFreeAPhiResidualOperatorBundle &operators,
                                                       bool enable_gauge_penalty, Real gauge_penalty_coefficient,
                                                       MatrixFreeAPhiResiduals &residuals)
{
    if (!enable_gauge_penalty)
    {
        return;
    }
    for (size_t i = 0; i != residuals.residual_ax.size(); ++i)
    {
        residuals.residual_ax[i] += gauge_penalty_coefficient * operators.gauge_penalty_gradient[i][0];
        residuals.residual_ay[i] += gauge_penalty_coefficient * operators.gauge_penalty_gradient[i][1];
        residuals.residual_az[i] += gauge_penalty_coefficient * operators.gauge_penalty_gradient[i][2];
    }
}

inline void accumulateMatrixFreeAPhiHighContrastResidualMetrics(
    const StdVec<uint8_t> &high_contrast_particle_mask, const MatrixFreeAPhiResidualOperatorBundle &operators,
    const MatrixFreeAPhiResiduals &residuals, MatrixFreeAPhiResiduals &metrics_target)
{
    if (metrics_target.high_contrast_particle_count == 0)
    {
        return;
    }

    Real sigma_grad_phi_squared_sum = 0.0;
    Real divergence_sigma_a_squared_sum = 0.0;
    Real residual_a_squared_sum = 0.0;
    Real residual_phi_squared_sum = 0.0;
    for (size_t i = 0; i != high_contrast_particle_mask.size(); ++i)
    {
        if (high_contrast_particle_mask[i] == 0)
        {
            continue;
        }
        sigma_grad_phi_squared_sum += std::norm(operators.sigma_grad_phi[i][0]) +
                                      std::norm(operators.sigma_grad_phi[i][1]) +
                                      std::norm(operators.sigma_grad_phi[i][2]);
        divergence_sigma_a_squared_sum += std::norm(operators.divergence_sigma_a[i]);
        residual_a_squared_sum +=
            std::norm(residuals.residual_ax[i]) + std::norm(residuals.residual_ay[i]) + std::norm(residuals.residual_az[i]);
        residual_phi_squared_sum += std::norm(residuals.residual_phi[i]);
    }
    const Real mask_count = static_cast<Real>(metrics_target.high_contrast_particle_count);
    metrics_target.high_contrast_sigma_grad_phi_l2 = std::sqrt(sigma_grad_phi_squared_sum / mask_count);
    metrics_target.high_contrast_div_sigma_a_l2 = std::sqrt(divergence_sigma_a_squared_sum / mask_count);
    metrics_target.high_contrast_residual_a_l2 = std::sqrt(residual_a_squared_sum / mask_count);
    metrics_target.high_contrast_residual_phi_l2 = std::sqrt(residual_phi_squared_sum / mask_count);
}

inline void finalizeMatrixFreeAPhiResidualDiagnostics(
    const MatrixFreePairwiseGraph &graph, const MatrixFreeAPhiCurrentComponentFields &current_components,
    const MatrixFreeAPhiResidualOperatorSemantics &semantics,
    const MatrixFreeAPhiInterfaceContrastRegion &contrast_region, const MatrixFreeAPhiResidualOperatorBundle &operators,
    MatrixFreeAPhiResiduals &residuals, MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const StdVec<Complex> divergence_j_field = computeResidualDivergenceOfCurrentDensity(
        graph, current_components.x, current_components.y, current_components.z, semantics.particle_kernel,
        sph_native_context);
    for (size_t i = 0; i != residuals.divergence_j.size(); ++i)
    {
        residuals.divergence_j[i] = divergence_j_field[i];
    }

    residuals.residual_ax_l2 = computeScalarFieldL2(residuals.residual_ax);
    residuals.residual_ay_l2 = computeScalarFieldL2(residuals.residual_ay);
    residuals.residual_az_l2 = computeScalarFieldL2(residuals.residual_az);
    residuals.residual_phi_l2 = computeScalarFieldL2(residuals.residual_phi);
    residuals.divergence_a_l2 = computeScalarFieldL2(residuals.divergence_a);
    residuals.divergence_j_l2 = computeScalarFieldL2(residuals.divergence_j);

    accumulateMatrixFreeAPhiHighContrastResidualMetrics(contrast_region.high_contrast_particle_mask, operators,
                                                        residuals, residuals);
}

inline MatrixFreeAPhiResiduals evaluateMatrixFreeAPhiResiduals(const MatrixFreePairwiseGraph &graph,
                                                               const MatrixFreeAPhiFields &fields,
                                                               const StdVec<Real> &electrical_conductivity,
                                                               const StdVec<Real> &magnetic_reluctivity,
                                                               const MatrixFreeAPhiParameters &parameters,
                                                               const MatrixFreeAPhiSources &sources,
                                                               bool enable_gauge_penalty,
                                                               Real gauge_penalty_coefficient,
                                                               MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    MatrixFreeAPhiResiduals residuals;
    const size_t number_of_particles = fields.ax.size();
    residuals.resize(number_of_particles);

    const MatrixFreeAPhiResidualOperatorSemantics semantics =
        selectMatrixFreeAPhiResidualOperatorSemantics(sph_native_context);
    const MatrixFreeAPhiResidualOperatorBundle operators = assembleMatrixFreeAPhiResidualOperators(
        graph, fields, electrical_conductivity, magnetic_reluctivity, enable_gauge_penalty, semantics,
        sph_native_context);

    const Real interface_contrast_threshold = SMAX(parameters.interface_contrast_threshold, Real(1.0));
    const MatrixFreeAPhiInterfaceContrastRegion contrast_region =
        detectMatrixFreeAPhiInterfaceContrastRegion(graph, electrical_conductivity, interface_contrast_threshold);
    residuals.interface_edge_count = contrast_region.interface_edge_count;
    residuals.high_contrast_edge_count = contrast_region.high_contrast_edge_count;
    residuals.high_contrast_particle_count = contrast_region.high_contrast_particle_count;
    residuals.max_interface_sigma_ratio = contrast_region.max_interface_sigma_ratio;

    MatrixFreeAPhiCurrentComponentFields current_components;
    assembleMatrixFreeAPhiResidualFieldContributions(fields, electrical_conductivity, parameters, sources, operators,
                                                     residuals, current_components);
    applyMatrixFreeAPhiGaugePenaltyToResiduals(operators, enable_gauge_penalty, gauge_penalty_coefficient, residuals);
    finalizeMatrixFreeAPhiResidualDiagnostics(graph, current_components, semantics, contrast_region, operators,
                                              residuals, sph_native_context);

    return residuals;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_HPP
