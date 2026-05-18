#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include <algorithm>
#include <cstdlib>
#include <memory>
#if SPHINXSYS_USE_SYCL
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_sycl_queue.hpp"
#include "implementation_sycl.h"
#include <sycl/sycl.hpp>
#endif

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline MatrixFreePairwiseNeighborEntry::MatrixFreePairwiseNeighborEntry(size_t index_j, const Vecd &gradient_weight,
                                                                        Real diffusion_weight, bool is_contact)
    : index_j_(index_j), gradient_weight_(gradient_weight), diffusion_weight_(diffusion_weight), is_contact_(is_contact)
{
}

inline void MatrixFreeFlatPairwiseGraph::clear()
{
    row_offsets_.clear();
    column_indices_.clear();
    gradient_weights_.clear();
    diffusion_weights_.clear();
    contact_flags_.clear();
}

inline bool MatrixFreeFlatPairwiseGraph::isValidFor(size_t number_of_rows) const
{
    return row_offsets_.size() == number_of_rows + 1 &&
           column_indices_.size() == gradient_weights_.size() &&
           column_indices_.size() == diffusion_weights_.size() &&
           column_indices_.size() == contact_flags_.size();
}

inline void refreshMatrixFreeFlatPairwiseGraph(MatrixFreePairwiseGraph &graph)
{
    graph.flat_.clear();
    graph.flat_.row_offsets_.resize(graph.rows_.size() + 1, 0);

    size_t edge_count = 0;
    for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
    {
        graph.flat_.row_offsets_[index_i] = edge_count;
        edge_count += graph.rows_[index_i].neighbors_.size();
    }
    graph.flat_.row_offsets_[graph.rows_.size()] = edge_count;

    graph.flat_.column_indices_.reserve(edge_count);
    graph.flat_.gradient_weights_.reserve(edge_count);
    graph.flat_.diffusion_weights_.reserve(edge_count);
    graph.flat_.contact_flags_.reserve(edge_count);

    for (const MatrixFreePairwiseRow &row : graph.rows_)
    {
        for (const MatrixFreePairwiseNeighborEntry &neighbor : row.neighbors_)
        {
            graph.flat_.column_indices_.push_back(neighbor.index_j_);
            graph.flat_.gradient_weights_.push_back(neighbor.gradient_weight_);
            graph.flat_.diffusion_weights_.push_back(neighbor.diffusion_weight_);
            graph.flat_.contact_flags_.push_back(neighbor.is_contact_ ? uint8_t(1) : uint8_t(0));
        }
    }
}

inline MatrixFreePairwiseGraph buildMatrixFreePairwiseGraph(const MatrixFreeAPhiDiscreteView &discrete_view,
                                                            const MatrixFreeAPhiParameters &parameters)
{
    MatrixFreePairwiseGraph graph;
    graph.rows_.resize(discrete_view.number_of_particles);

    if (discrete_view.particle_configuration == nullptr || discrete_view.positions == nullptr ||
        discrete_view.volumetric_measure == nullptr)
    {
        return graph;
    }

    for (size_t index_i = 0; index_i != discrete_view.number_of_particles; ++index_i)
    {
        const Neighborhood &neighborhood = (*discrete_view.particle_configuration)[index_i];
        MatrixFreePairwiseRow &row = graph.rows_[index_i];
        row.neighbors_.reserve(neighborhood.current_size_);

        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            const size_t index_j = neighborhood.j_[n];
            if (neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }

            const bool is_contact_neighbor =
                discrete_view.neighbor_is_contact != nullptr && index_i < discrete_view.neighbor_is_contact->size() &&
                n < (*discrete_view.neighbor_is_contact)[index_i].size() &&
                (*discrete_view.neighbor_is_contact)[index_i][n] != 0;

            const Real gradient_scale = is_contact_neighbor ? parameters.contact_gradient_scale : Real(1.0);
            const Real diffusion_scale = is_contact_neighbor ? parameters.contact_diffusion_scale : Real(1.0);

            Vecd gradient_weight =
                neighborhood.dW_ij_[n] * discrete_view.volumetric_measure[index_j] * neighborhood.e_ij_[n];
            gradient_weight *= gradient_scale;
            if (!std::isfinite(gradient_weight.squaredNorm()))
            {
                continue;
            }

            Real smoothing_length_i = discrete_view.reference_smoothing_length;
            Real smoothing_length_j = discrete_view.reference_smoothing_length;
            if (discrete_view.smoothing_length_ratio != nullptr)
            {
                smoothing_length_i *= discrete_view.smoothing_length_ratio[index_i];
                smoothing_length_j *= discrete_view.smoothing_length_ratio[index_j];
            }
            const Real pair_smoothing_length = SMAX(smoothing_length_i, smoothing_length_j);
            const Real dW_ijV_j = neighborhood.dW_ij_[n] * discrete_view.volumetric_measure[index_j];
            const Real diffusion_weight =
                diffusion_scale * (-2.0 * dW_ijV_j /
                                   (neighborhood.r_ij_[n] + parameters.pair_weight_regularization * pair_smoothing_length +
                                    TinyReal));

            row.neighbors_.emplace_back(index_j, gradient_weight, diffusion_weight, is_contact_neighbor);
        }
    }

    refreshMatrixFreeFlatPairwiseGraph(graph);
    return graph;
}

#if SPHINXSYS_USE_SYCL
namespace
{
bool g_matrix_free_laplace_residual_use_sycl = false;

struct MatrixFreeSyclFlatGraphBuffers
{
    size_t row_count_ = 0;
    size_t edge_count_ = 0;
    size_t upload_count_ = 0;
    const size_t *row_offsets_source_ = nullptr;
    const size_t *column_indices_source_ = nullptr;
    const Vecd *gradient_weights_source_ = nullptr;
    const Real *diffusion_weights_source_ = nullptr;
    std::unique_ptr<sycl::buffer<size_t, 1>> row_offsets_buffer_;
    std::unique_ptr<sycl::buffer<size_t, 1>> column_indices_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_weight_x_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_weight_y_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_weight_z_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> diffusion_weights_buffer_;

    size_t *dev_row_offsets_ = nullptr;
    size_t *dev_column_indices_ = nullptr;
    Real *dev_gradient_weight_x_ = nullptr;
    Real *dev_gradient_weight_y_ = nullptr;
    Real *dev_gradient_weight_z_ = nullptr;
    Real *dev_diffusion_weights_ = nullptr;
    size_t dev_topology_rows_ = 0;
    size_t dev_topology_edges_ = 0;

    void releaseDeviceTopologyUsm()
    {
        if (dev_row_offsets_)
        {
            freeDeviceData(dev_row_offsets_);
            dev_row_offsets_ = nullptr;
        }
        if (dev_column_indices_)
        {
            freeDeviceData(dev_column_indices_);
            dev_column_indices_ = nullptr;
        }
        if (dev_gradient_weight_x_)
        {
            freeDeviceData(dev_gradient_weight_x_);
            dev_gradient_weight_x_ = nullptr;
        }
        if (dev_gradient_weight_y_)
        {
            freeDeviceData(dev_gradient_weight_y_);
            dev_gradient_weight_y_ = nullptr;
        }
        if (dev_gradient_weight_z_)
        {
            freeDeviceData(dev_gradient_weight_z_);
            dev_gradient_weight_z_ = nullptr;
        }
        if (dev_diffusion_weights_)
        {
            freeDeviceData(dev_diffusion_weights_);
            dev_diffusion_weights_ = nullptr;
        }
        dev_topology_rows_ = 0;
        dev_topology_edges_ = 0;
    }

    void ensureDeviceTopologyUsmMirror(const MatrixFreeFlatPairwiseGraph &flat_graph, size_t row_count)
    {
        if (!matrixFreeSyclGraphTopologyUsesDeviceUsm())
        {
            releaseDeviceTopologyUsm();
            return;
        }
        const size_t edge_count = flat_graph.column_indices_.size();
        const bool layout_ok = dev_row_offsets_ && dev_column_indices_ && dev_gradient_weight_x_ &&
                               dev_gradient_weight_y_ && dev_gradient_weight_z_ && dev_diffusion_weights_ &&
                               dev_topology_rows_ == row_count && dev_topology_edges_ == edge_count;
        if (!layout_ok)
        {
            releaseDeviceTopologyUsm();
            dev_row_offsets_ = allocateDeviceOnly<size_t>(row_count + 1);
            dev_column_indices_ = allocateDeviceOnly<size_t>(edge_count);
            dev_gradient_weight_x_ = allocateDeviceOnly<Real>(edge_count);
            dev_gradient_weight_y_ = allocateDeviceOnly<Real>(edge_count);
            dev_gradient_weight_z_ = allocateDeviceOnly<Real>(edge_count);
            dev_diffusion_weights_ = allocateDeviceOnly<Real>(edge_count);
            dev_topology_rows_ = row_count;
            dev_topology_edges_ = edge_count;
        }
        copyToDevice(flat_graph.row_offsets_.data(), dev_row_offsets_, row_count + 1);
        copyToDevice(flat_graph.column_indices_.data(), dev_column_indices_, edge_count);
        copyToDevice(flat_graph.diffusion_weights_.data(), dev_diffusion_weights_, edge_count);
        StdVec<Real> gx(edge_count), gy(edge_count), gz(edge_count);
        for (size_t e = 0; e != edge_count; ++e)
        {
            gx[e] = flat_graph.gradient_weights_[e][0];
            gy[e] = flat_graph.gradient_weights_[e][1];
            gz[e] = flat_graph.gradient_weights_[e][2];
        }
        copyToDevice(gx.data(), dev_gradient_weight_x_, edge_count);
        copyToDevice(gy.data(), dev_gradient_weight_y_, edge_count);
        copyToDevice(gz.data(), dev_gradient_weight_z_, edge_count);
    }

    bool deviceTopologyUsmReady() const
    {
        return dev_row_offsets_ != nullptr && matrixFreeSyclGraphTopologyUsesDeviceUsm();
    }

    const size_t *deviceTopologyRowOffsets() const { return dev_row_offsets_; }
    const size_t *deviceTopologyColumnIndices() const { return dev_column_indices_; }
    const Real *deviceTopologyDiffusionWeights() const { return dev_diffusion_weights_; }
    const Real *deviceTopologyGradientWeightX() const { return dev_gradient_weight_x_; }
    const Real *deviceTopologyGradientWeightY() const { return dev_gradient_weight_y_; }
    const Real *deviceTopologyGradientWeightZ() const { return dev_gradient_weight_z_; }

    void resize(size_t row_count, size_t edge_count)
    {
        if (row_count_ == row_count && edge_count_ == edge_count && row_offsets_buffer_ && column_indices_buffer_ &&
            gradient_weight_x_buffer_ && gradient_weight_y_buffer_ && gradient_weight_z_buffer_ &&
            diffusion_weights_buffer_)
        {
            return;
        }

        releaseDeviceTopologyUsm();

        row_count_ = row_count;
        edge_count_ = edge_count;
        row_offsets_source_ = nullptr;
        column_indices_source_ = nullptr;
        gradient_weights_source_ = nullptr;
        diffusion_weights_source_ = nullptr;
        row_offsets_buffer_ = std::make_unique<sycl::buffer<size_t, 1>>(sycl::range<1>(row_count + 1));
        column_indices_buffer_ = std::make_unique<sycl::buffer<size_t, 1>>(sycl::range<1>(edge_count));
        gradient_weight_x_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(edge_count));
        gradient_weight_y_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(edge_count));
        gradient_weight_z_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(edge_count));
        diffusion_weights_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(edge_count));
    }

    bool matches(const MatrixFreeFlatPairwiseGraph &flat_graph, size_t row_count) const
    {
        return row_count_ == row_count && edge_count_ == flat_graph.column_indices_.size() &&
               row_offsets_source_ == flat_graph.row_offsets_.data() &&
               column_indices_source_ == flat_graph.column_indices_.data() &&
               gradient_weights_source_ == flat_graph.gradient_weights_.data() &&
               diffusion_weights_source_ == flat_graph.diffusion_weights_.data() &&
               row_offsets_buffer_ && column_indices_buffer_ && gradient_weight_x_buffer_ &&
               gradient_weight_y_buffer_ && gradient_weight_z_buffer_ && diffusion_weights_buffer_;
    }

    void loadIfNeeded(const MatrixFreeFlatPairwiseGraph &flat_graph, size_t row_count)
    {
        const size_t edge_count = flat_graph.column_indices_.size();
        if (matches(flat_graph, row_count))
        {
            ensureDeviceTopologyUsmMirror(flat_graph, row_count);
            return;
        }
        resize(row_count, edge_count);

        sycl::host_accessor<size_t, 1, sycl::access::mode::write> row_offsets(*row_offsets_buffer_);
        sycl::host_accessor<size_t, 1, sycl::access::mode::write> column_indices(*column_indices_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> gradient_weight_x(*gradient_weight_x_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> gradient_weight_y(*gradient_weight_y_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> gradient_weight_z(*gradient_weight_z_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> diffusion_weights(*diffusion_weights_buffer_);

        for (size_t i = 0; i != row_count + 1; ++i)
        {
            row_offsets[i] = flat_graph.row_offsets_[i];
        }
        for (size_t edge = 0; edge != edge_count; ++edge)
        {
            column_indices[edge] = flat_graph.column_indices_[edge];
            gradient_weight_x[edge] = flat_graph.gradient_weights_[edge][0];
            gradient_weight_y[edge] = flat_graph.gradient_weights_[edge][1];
            gradient_weight_z[edge] = flat_graph.gradient_weights_[edge][2];
            diffusion_weights[edge] = flat_graph.diffusion_weights_[edge];
        }

        row_offsets_source_ = flat_graph.row_offsets_.data();
        column_indices_source_ = flat_graph.column_indices_.data();
        gradient_weights_source_ = flat_graph.gradient_weights_.data();
        diffusion_weights_source_ = flat_graph.diffusion_weights_.data();
        ++upload_count_;
        ensureDeviceTopologyUsmMirror(flat_graph, row_count);
    }
};

struct MatrixFreeSyclScalarValueBuffers
{
    size_t row_count_ = 0;
    std::unique_ptr<sycl::buffer<Real, 1>> field_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> field_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> coefficient_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_x_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_x_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_y_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_y_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_z_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> gradient_z_imag_buffer_;
    /** Second gradient slot for fused phi standard + harmonic weighted gradient. */
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_x_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_x_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_y_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_y_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_z_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> harmonic_gradient_z_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> laplace_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> laplace_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> diagonal_scale_buffer_;
    size_t field_upload_count_ = 0;
    size_t field_upload_skip_count_ = 0;
    size_t coefficient_upload_count_ = 0;
    size_t laplace_upload_count_ = 0;
    size_t gradient_download_count_ = 0;
    size_t laplace_download_count_ = 0;
    const Complex *field_source_ = nullptr;
    size_t field_source_size_ = 0;
    const Real *coefficient_source_ = nullptr;
    size_t coefficient_source_size_ = 0;

    Real *dev_field_real_ = nullptr;
    Real *dev_field_imag_ = nullptr;

    void releaseFieldUsm()
    {
        if (dev_field_real_)
        {
            freeDeviceData(dev_field_real_);
            dev_field_real_ = nullptr;
        }
        if (dev_field_imag_)
        {
            freeDeviceData(dev_field_imag_);
            dev_field_imag_ = nullptr;
        }
    }

    bool deviceFieldUsmReady() const
    {
        return dev_field_real_ != nullptr && matrixFreeSyclFieldValuesUseDeviceUsm();
    }

    void syncHostFieldBuffersFromDeviceUsm()
    {
        if (!deviceFieldUsmReady())
        {
            return;
        }
        StdVec<Real> re(row_count_), im(row_count_);
        copyFromDevice(re.data(), dev_field_real_, row_count_);
        copyFromDevice(im.data(), dev_field_imag_, row_count_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> fr(*field_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> fi(*field_imag_buffer_);
        for (size_t i = 0; i != row_count_; ++i)
        {
            fr[i] = re[i];
            fi[i] = im[i];
        }
    }

    /** Non-null only when `deviceFieldUsmReady()`; used by graph Helmholtz Jacobi SYCL path. */
    Real *mutableDeviceFieldRealUsmForJacobi()
    {
        return deviceFieldUsmReady() ? dev_field_real_ : nullptr;
    }
    Real *mutableDeviceFieldImagUsmForJacobi()
    {
        return deviceFieldUsmReady() ? dev_field_imag_ : nullptr;
    }

    const Real *deviceFieldRealUsmRead() const
    {
        return deviceFieldUsmReady() ? dev_field_real_ : nullptr;
    }
    const Real *deviceFieldImagUsmRead() const
    {
        return deviceFieldUsmReady() ? dev_field_imag_ : nullptr;
    }

    void resize(size_t row_count)
    {
        if (row_count_ == row_count && field_real_buffer_ && field_imag_buffer_ && coefficient_buffer_ &&
            gradient_x_real_buffer_ && gradient_x_imag_buffer_ && gradient_y_real_buffer_ &&
            gradient_y_imag_buffer_ && gradient_z_real_buffer_ && gradient_z_imag_buffer_ &&
            harmonic_gradient_x_real_buffer_ && harmonic_gradient_x_imag_buffer_ &&
            harmonic_gradient_y_real_buffer_ && harmonic_gradient_y_imag_buffer_ &&
            harmonic_gradient_z_real_buffer_ && harmonic_gradient_z_imag_buffer_ &&
            laplace_real_buffer_ && laplace_imag_buffer_ && diagonal_scale_buffer_)
        {
            return;
        }

        releaseFieldUsm();

        row_count_ = row_count;
        field_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        field_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        coefficient_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_x_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_x_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_y_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_y_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_z_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        gradient_z_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_x_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_x_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_y_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_y_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_z_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        harmonic_gradient_z_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        laplace_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        laplace_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        diagonal_scale_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(row_count));
        field_source_ = nullptr;
        field_source_size_ = 0;
        coefficient_source_ = nullptr;
        coefficient_source_size_ = 0;
        if (matrixFreeSyclFieldValuesUseDeviceUsm())
        {
            dev_field_real_ = allocateDeviceOnly<Real>(row_count);
            dev_field_imag_ = allocateDeviceOnly<Real>(row_count);
        }
    }

    bool fieldMatches(const StdVec<Complex> &field) const
    {
        return field_source_ == field.data() && field_source_size_ == field.size();
    }

    void loadField(const StdVec<Complex> &field)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::write> field_real(*field_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> field_imag(*field_imag_buffer_);
        for (size_t i = 0; i != field.size(); ++i)
        {
            field_real[i] = field[i].real();
            field_imag[i] = field[i].imag();
        }
        field_source_ = field.data();
        field_source_size_ = field.size();
        ++field_upload_count_;
        if (deviceFieldUsmReady())
        {
            StdVec<Real> re(field.size()), im(field.size());
            for (size_t i = 0; i != field.size(); ++i)
            {
                re[i] = field_real[i];
                im[i] = field_imag[i];
            }
            copyToDevice(re.data(), dev_field_real_, field.size());
            copyToDevice(im.data(), dev_field_imag_, field.size());
        }
    }

    void loadFieldIfNeeded(const StdVec<Complex> &field)
    {
        if (fieldMatches(field))
        {
            ++field_upload_skip_count_;
            return;
        }
        loadField(field);
    }

    void markFieldCurrentForSource(const StdVec<Complex> &field)
    {
        field_source_ = field.data();
        field_source_size_ = field.size();
    }

    void invalidateFieldCache()
    {
        field_source_ = nullptr;
        field_source_size_ = 0;
    }

    void loadCoefficientIfNeeded(const StdVec<Real> &coefficient)
    {
        if (coefficient_source_ == coefficient.data() && coefficient_source_size_ == coefficient.size())
        {
            return;
        }
        sycl::host_accessor<Real, 1, sycl::access::mode::write> coefficient_values(*coefficient_buffer_);
        for (size_t i = 0; i != coefficient.size(); ++i)
        {
            coefficient_values[i] = coefficient[i];
        }
        coefficient_source_ = coefficient.data();
        coefficient_source_size_ = coefficient.size();
        ++coefficient_upload_count_;
    }

    void loadLaplaceState(const ScalarComplexHelmholtzResiduals &residuals)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::write> laplace_real(*laplace_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> laplace_imag(*laplace_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> diagonal_scale(*diagonal_scale_buffer_);
        for (size_t i = 0; i != residuals.laplace_term_.size(); ++i)
        {
            laplace_real[i] = residuals.laplace_term_[i].real();
            laplace_imag[i] = residuals.laplace_term_[i].imag();
            diagonal_scale[i] = residuals.diagonal_scale_[i];
        }
        ++laplace_upload_count_;
    }

    void readGradient(StdVec<Vec3c> &gradient)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_x_real(*gradient_x_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_x_imag(*gradient_x_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_y_real(*gradient_y_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_y_imag(*gradient_y_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_z_real(*gradient_z_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> gradient_z_imag(*gradient_z_imag_buffer_);
        for (size_t i = 0; i != gradient.size(); ++i)
        {
            gradient[i][0] = Complex(gradient_x_real[i], gradient_x_imag[i]);
            gradient[i][1] = Complex(gradient_y_real[i], gradient_y_imag[i]);
            gradient[i][2] = Complex(gradient_z_real[i], gradient_z_imag[i]);
        }
        ++gradient_download_count_;
    }

    /** Read standard gradient from primary buffers and harmonic weighted gradient from secondary; one diagnostic tick. */
    void readFusedStandardAndHarmonicGradients(StdVec<Vec3c> &standard_gradient, StdVec<Vec3c> &harmonic_gradient)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::read> sxr(*gradient_x_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> sxi(*gradient_x_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> syr(*gradient_y_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> syi(*gradient_y_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> szr(*gradient_z_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> szi(*gradient_z_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hxr(*harmonic_gradient_x_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hxi(*harmonic_gradient_x_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hyr(*harmonic_gradient_y_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hyi(*harmonic_gradient_y_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hzr(*harmonic_gradient_z_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> hzi(*harmonic_gradient_z_imag_buffer_);
        for (size_t i = 0; i != standard_gradient.size(); ++i)
        {
            standard_gradient[i][0] = Complex(sxr[i], sxi[i]);
            standard_gradient[i][1] = Complex(syr[i], syi[i]);
            standard_gradient[i][2] = Complex(szr[i], szi[i]);
            harmonic_gradient[i][0] = Complex(hxr[i], hxi[i]);
            harmonic_gradient[i][1] = Complex(hyr[i], hyi[i]);
            harmonic_gradient[i][2] = Complex(hzr[i], hzi[i]);
        }
        ++gradient_download_count_;
    }

    void readLaplaceState(ScalarComplexHelmholtzResiduals &residuals)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::read> laplace_real(*laplace_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> laplace_imag(*laplace_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> diagonal_scale(*diagonal_scale_buffer_);
        for (size_t i = 0; i != residuals.laplace_term_.size(); ++i)
        {
            residuals.laplace_term_[i] = Complex(laplace_real[i], laplace_imag[i]);
            residuals.diagonal_scale_[i] = diagonal_scale[i];
        }
        ++laplace_download_count_;
    }
};

/** USM staging for `tryMatrixFreeDivergenceAndGradientOfDivergenceOfVectorPotentialSycl`: six real/imag A
 *  components plus divergence real/imag on device when `matrixFreeSyclFieldValuesUseDeviceUsm()` is true. */
struct MatrixFreeSyclVectorPotentialDivergenceStaging
{
    size_t rows_ = 0;
    Real *d_axr_ = nullptr;
    Real *d_axi_ = nullptr;
    Real *d_ayr_ = nullptr;
    Real *d_ayi_ = nullptr;
    Real *d_azr_ = nullptr;
    Real *d_azi_ = nullptr;
    Real *d_div_r_ = nullptr;
    Real *d_div_i_ = nullptr;

    void releaseAll()
    {
        auto free_one = [](Real *&p) {
            if (p)
            {
                freeDeviceData(p);
                p = nullptr;
            }
        };
        free_one(d_axr_);
        free_one(d_axi_);
        free_one(d_ayr_);
        free_one(d_ayi_);
        free_one(d_azr_);
        free_one(d_azi_);
        free_one(d_div_r_);
        free_one(d_div_i_);
        rows_ = 0;
    }

    void ensureForRowCount(size_t n)
    {
        if (!matrixFreeSyclFieldValuesUseDeviceUsm())
        {
            releaseAll();
            return;
        }
        if (n == 0)
        {
            releaseAll();
            return;
        }
        if (rows_ == n && d_axr_ != nullptr)
        {
            return;
        }
        releaseAll();
        rows_ = n;
        d_axr_ = allocateDeviceOnly<Real>(n);
        d_axi_ = allocateDeviceOnly<Real>(n);
        d_ayr_ = allocateDeviceOnly<Real>(n);
        d_ayi_ = allocateDeviceOnly<Real>(n);
        d_azr_ = allocateDeviceOnly<Real>(n);
        d_azi_ = allocateDeviceOnly<Real>(n);
        d_div_r_ = allocateDeviceOnly<Real>(n);
        d_div_i_ = allocateDeviceOnly<Real>(n);
    }

    bool ready() const
    {
        return d_axr_ != nullptr && d_div_r_ != nullptr && rows_ > 0;
    }

    void uploadAFromHost(const StdVec<Real> &axr, const StdVec<Real> &axi, const StdVec<Real> &ayr,
                         const StdVec<Real> &ayi, const StdVec<Real> &azr, const StdVec<Real> &azi)
    {
        if (!ready() || axr.size() != rows_)
        {
            return;
        }
        copyToDevice(axr.data(), d_axr_, rows_);
        copyToDevice(axi.data(), d_axi_, rows_);
        copyToDevice(ayr.data(), d_ayr_, rows_);
        copyToDevice(ayi.data(), d_ayi_, rows_);
        copyToDevice(azr.data(), d_azr_, rows_);
        copyToDevice(azi.data(), d_azi_, rows_);
    }

    void readDivToComplexHost(StdVec<Complex> &out) const
    {
        if (!ready() || out.size() != rows_)
        {
            return;
        }
        StdVec<Real> re(rows_), im(rows_);
        copyFromDevice(re.data(), d_div_r_, rows_);
        copyFromDevice(im.data(), d_div_i_, rows_);
        for (size_t i = 0; i != rows_; ++i)
        {
            out[i] = Complex(re[i], im[i]);
        }
    }

    Real *mutableDivReal() { return d_div_r_; }
    Real *mutableDivImag() { return d_div_i_; }
    const Real *axrRead() const { return d_axr_; }
    const Real *axiRead() const { return d_axi_; }
    const Real *ayrRead() const { return d_ayr_; }
    const Real *ayiRead() const { return d_ayi_; }
    const Real *azrRead() const { return d_azr_; }
    const Real *aziRead() const { return d_azi_; }
    const Real *divRealRead() const { return d_div_r_; }
    const Real *divImagRead() const { return d_div_i_; }
};

struct MatrixFreeAPhiSyclWorkspace
{
    MatrixFreeSyclFlatGraphBuffers graph_buffers_;
    MatrixFreeSyclScalarValueBuffers scalar_value_buffers_;
    MatrixFreeSyclVectorPotentialDivergenceStaging vector_potential_divergence_staging_;
};

inline MatrixFreeAPhiSyclWorkspace &matrixFreeAPhiSyclWorkspace()
{
    static MatrixFreeAPhiSyclWorkspace workspace;
    return workspace;
}

inline bool isFlatGraphFieldCompatible(const MatrixFreePairwiseGraph &graph, size_t field_size)
{
    return graph.flat_.isValidFor(graph.rows_.size()) && field_size == graph.rows_.size();
}

inline void prepareMatrixFreeAPhiSyclGraphWorkspaceImpl(const MatrixFreePairwiseGraph &graph)
{
    if (!graph.flat_.isValidFor(graph.rows_.size()))
    {
        return;
    }
    matrixFreeAPhiSyclWorkspace().graph_buffers_.loadIfNeeded(graph.flat_, graph.rows_.size());
}

inline void accumulateScalarLaplaceResidualsFromGraphSyclWithWorkspace(const MatrixFreePairwiseGraph &graph,
                                                                      const StdVec<Complex> &field,
                                                                      const StdVec<Real> &diffusion_coefficient,
                                                                      ScalarComplexHelmholtzResiduals &residuals,
                                                                      MatrixFreeAPhiSyclWorkspace &sycl_workspace,
                                                                      bool laplace_state_is_zero = false)
{
    const size_t number_of_rows = graph.rows_.size();
    if (!isFlatGraphFieldCompatible(graph, field.size()) || diffusion_coefficient.size() != number_of_rows ||
        residuals.laplace_term_.size() != number_of_rows || residuals.diagonal_scale_.size() != number_of_rows)
    {
        return;
    }

    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);
    value_buffers.loadFieldIfNeeded(field);
    value_buffers.loadCoefficientIfNeeded(diffusion_coefficient);
    if (!laplace_state_is_zero)
    {
        value_buffers.loadLaplaceState(residuals);
    }

    const Real *field_usm_r = value_buffers.deviceFieldRealUsmRead();
    const Real *field_usm_i = value_buffers.deviceFieldImagUsmRead();
    const bool field_reads_from_usm = field_usm_r != nullptr && field_usm_i != nullptr;

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *edge_w_dev = graph_buffers.deviceTopologyDiffusionWeights();
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto laplace_r = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto laplace_i = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                const bool initialize_from_zero = laplace_state_is_zero;
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real local_laplace_real = initialize_from_zero ? Real(0.0) : laplace_r[index_i];
                    Real local_laplace_imag = initialize_from_zero ? Real(0.0) : laplace_i[index_i];
                    Real local_diagonal = initialize_from_zero ? Real(0.0) : diagonal[index_i];
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real coefficient_ij =
                            2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                            (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                        const Real pair_weight = coefficient_ij * edge_w_dev[edge];
                        local_laplace_real += pair_weight * (value_i_real - f_fr[index_j]);
                        local_laplace_imag += pair_weight * (value_i_imag - f_fi[index_j]);
                        local_diagonal += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                    }
                    laplace_r[index_i] = local_laplace_real;
                    laplace_i[index_i] = local_laplace_imag;
                    diagonal[index_i] = local_diagonal;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto laplace_r = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto laplace_i = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                const bool initialize_from_zero = laplace_state_is_zero;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real local_laplace_real = initialize_from_zero ? Real(0.0) : laplace_r[index_i];
                    Real local_laplace_imag = initialize_from_zero ? Real(0.0) : laplace_i[index_i];
                    Real local_diagonal = initialize_from_zero ? Real(0.0) : diagonal[index_i];
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real coefficient_ij =
                            2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                            (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                        const Real pair_weight = coefficient_ij * edge_w_dev[edge];
                        local_laplace_real += pair_weight * (value_i_real - f_real[index_j]);
                        local_laplace_imag += pair_weight * (value_i_imag - f_imag[index_j]);
                        local_diagonal += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                    }
                    laplace_r[index_i] = local_laplace_real;
                    laplace_i[index_i] = local_laplace_imag;
                    diagonal[index_i] = local_diagonal;
                });
            });
        }
    }
    else
    {
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto edge_weights = graph_buffers.diffusion_weights_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto laplace_r = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto laplace_i = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                const bool initialize_from_zero = laplace_state_is_zero;
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real local_laplace_real = initialize_from_zero ? Real(0.0) : laplace_r[index_i];
                    Real local_laplace_imag = initialize_from_zero ? Real(0.0) : laplace_i[index_i];
                    Real local_diagonal = initialize_from_zero ? Real(0.0) : diagonal[index_i];
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real coefficient_ij =
                            2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                            (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                        const Real pair_weight = coefficient_ij * edge_weights[edge];
                        local_laplace_real += pair_weight * (value_i_real - f_fr[index_j]);
                        local_laplace_imag += pair_weight * (value_i_imag - f_fi[index_j]);
                        local_diagonal += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                    }
                    laplace_r[index_i] = local_laplace_real;
                    laplace_i[index_i] = local_laplace_imag;
                    diagonal[index_i] = local_diagonal;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto diffusion_values = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto edge_weights = graph_buffers.diffusion_weights_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto laplace_r = value_buffers.laplace_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto laplace_i = value_buffers.laplace_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                auto diagonal = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read_write>(cgh);
                const bool initialize_from_zero = laplace_state_is_zero;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real local_laplace_real = initialize_from_zero ? Real(0.0) : laplace_r[index_i];
                    Real local_laplace_imag = initialize_from_zero ? Real(0.0) : laplace_i[index_i];
                    Real local_diagonal = initialize_from_zero ? Real(0.0) : diagonal[index_i];
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real coefficient_ij =
                            2.0 * diffusion_values[index_i] * diffusion_values[index_j] /
                            (diffusion_values[index_i] + diffusion_values[index_j] + TinyReal);
                        const Real pair_weight = coefficient_ij * edge_weights[edge];
                        local_laplace_real += pair_weight * (value_i_real - f_real[index_j]);
                        local_laplace_imag += pair_weight * (value_i_imag - f_imag[index_j]);
                        local_diagonal += pair_weight >= 0.0 ? pair_weight : -pair_weight;
                    }
                    laplace_r[index_i] = local_laplace_real;
                    laplace_i[index_i] = local_laplace_imag;
                    diagonal[index_i] = local_diagonal;
                });
            });
        }
    }

    value_buffers.readLaplaceState(residuals);
}

inline void accumulateScalarLaplaceResidualsFromGraphSycl(const MatrixFreePairwiseGraph &graph,
                                                         const StdVec<Complex> &field,
                                                         const StdVec<Real> &diffusion_coefficient,
                                                         ScalarComplexHelmholtzResiduals &residuals)
{
    accumulateScalarLaplaceResidualsFromGraphSyclWithWorkspace(graph, field, diffusion_coefficient, residuals,
                                                              matrixFreeAPhiSyclWorkspace());
}

inline void accumulateScalarLaplaceResidualsFromClearedGraphSycl(const MatrixFreePairwiseGraph &graph,
                                                                const StdVec<Complex> &field,
                                                                const StdVec<Real> &diffusion_coefficient,
                                                                ScalarComplexHelmholtzResiduals &residuals)
{
    accumulateScalarLaplaceResidualsFromGraphSyclWithWorkspace(graph, field, diffusion_coefficient, residuals,
                                                              matrixFreeAPhiSyclWorkspace(), true);
}
} // namespace

inline void prepareMatrixFreeAPhiSyclGraphWorkspace(const MatrixFreePairwiseGraph &graph)
{
    prepareMatrixFreeAPhiSyclGraphWorkspaceImpl(graph);
}

/** Call once per timestep (or after topology resize) from the driver that owns `graph` and particle count:
 *  uploads flat graph if needed and allocates scalar workspace rows so later kernels can skip resize churn. */
inline void prepareMatrixFreeAPhiSyclDeviceResources(const MatrixFreePairwiseGraph &graph, size_t scalar_row_count)
{
    prepareMatrixFreeAPhiSyclGraphWorkspace(graph);
    matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.resize(scalar_row_count);
    matrixFreeAPhiSyclWorkspace().vector_potential_divergence_staging_.ensureForRowCount(scalar_row_count);
}

inline size_t matrixFreeAPhiSyclGraphWorkspaceUploadCount()
{
    return matrixFreeAPhiSyclWorkspace().graph_buffers_.upload_count_;
}

inline size_t matrixFreeAPhiSyclValueFieldUploadCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.field_upload_count_;
}

inline size_t matrixFreeAPhiSyclValueCoefficientUploadCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.coefficient_upload_count_;
}

inline size_t matrixFreeAPhiSyclValueLaplaceUploadCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.laplace_upload_count_;
}

inline size_t matrixFreeAPhiSyclValueGradientDownloadCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.gradient_download_count_;
}

inline size_t matrixFreeAPhiSyclValueLaplaceDownloadCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.laplace_download_count_;
}

inline size_t matrixFreeAPhiSyclValueFieldUploadSkipCount()
{
    return matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.field_upload_skip_count_;
}

inline void invalidateMatrixFreeAPhiSyclValueFieldCache()
{
    matrixFreeAPhiSyclWorkspace().scalar_value_buffers_.invalidateFieldCache();
}

inline bool matrixFreeSyclGraphTopologyUsesDeviceUsm()
{
    static const bool use = []() {
        const char *e = std::getenv("EM_APHI_MATRIX_FREE_SYCL_GRAPH_TOPOLOGY_USM");
        return e != nullptr && e[0] == '1' && e[1] == '\0';
    }();
    return use;
}

inline void setMatrixFreeLaplaceResidualUseSycl(bool enable)
{
    g_matrix_free_laplace_residual_use_sycl = enable;
}

inline bool matrixFreeLaplaceResidualUseSycl()
{
    return g_matrix_free_laplace_residual_use_sycl;
}
#endif

inline void accumulateScalarLaplaceResidualsFromGraph(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                      const StdVec<Real> &diffusion_coefficient,
                                                      ScalarComplexHelmholtzResiduals &residuals)
{
#if SPHINXSYS_USE_SYCL
    if (matrixFreeLaplaceResidualUseSycl())
    {
        accumulateScalarLaplaceResidualsFromGraphSycl(graph, field, diffusion_coefficient, residuals);
        return;
    }
#endif
    if (graph.flat_.isValidFor(graph.rows_.size()))
    {
        for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
        {
            const Complex value_i = field[index_i];
            for (size_t edge = graph.flat_.row_offsets_[index_i]; edge != graph.flat_.row_offsets_[index_i + 1]; ++edge)
            {
                const size_t index_j = graph.flat_.column_indices_[edge];
                const Real coefficient_ij = harmonicMean(diffusion_coefficient[index_i], diffusion_coefficient[index_j]);
                const Real pair_weight = coefficient_ij * graph.flat_.diffusion_weights_[edge];
                accumulateDirectionalScalarLaplaceContribution(index_i, value_i, field[index_j], pair_weight, residuals);
            }
        }
        return;
    }

    for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
    {
        const Complex value_i = field[index_i];
        for (const MatrixFreePairwiseNeighborEntry &neighbor : graph.rows_[index_i].neighbors_)
        {
            const size_t index_j = neighbor.index_j_;
            const Real coefficient_ij = harmonicMean(diffusion_coefficient[index_i], diffusion_coefficient[index_j]);
            const Real pair_weight = coefficient_ij * neighbor.diffusion_weight_;
            accumulateDirectionalScalarLaplaceContribution(index_i, value_i, field[index_j], pair_weight, residuals);
        }
    }
}

inline void accumulateScalarLaplaceResidualsFromClearedGraph(const MatrixFreePairwiseGraph &graph,
                                                            const StdVec<Complex> &field,
                                                            const StdVec<Real> &diffusion_coefficient,
                                                            ScalarComplexHelmholtzResiduals &residuals)
{
#if SPHINXSYS_USE_SYCL
    if (matrixFreeLaplaceResidualUseSycl())
    {
        accumulateScalarLaplaceResidualsFromClearedGraphSycl(graph, field, diffusion_coefficient, residuals);
        return;
    }
#endif
    accumulateScalarLaplaceResidualsFromGraph(graph, field, diffusion_coefficient, residuals);
}

inline StdVec<Complex> applyScalarNegativeLaplaceFromGraph(const MatrixFreePairwiseGraph &graph,
                                                           const StdVec<Complex> &field,
                                                           const StdVec<Real> &diffusion_coefficient)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(field.size());
    residuals.clear();
    accumulateScalarLaplaceResidualsFromClearedGraph(graph, field, diffusion_coefficient, residuals);
    return residuals.laplace_term_;
}

#if SPHINXSYS_USE_SYCL
namespace
{
bool g_matrix_free_gradient_use_sycl = false;
bool g_matrix_free_harmonic_gradient_use_sycl = false;

inline StdVec<Vec3c> applyMatrixFreeGradientSyclWithWorkspace(const MatrixFreePairwiseGraph &graph,
                                                             const StdVec<Complex> &field,
                                                             MatrixFreeAPhiSyclWorkspace &sycl_workspace)
{
    const size_t number_of_rows = graph.rows_.size();
    StdVec<Vec3c> gradient(field.size(), Vec3c::Zero());
    if (!graph.flat_.isValidFor(number_of_rows) || field.size() != number_of_rows)
    {
        return gradient;
    }

    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);
    value_buffers.loadFieldIfNeeded(field);

    const Real *field_usm_r = value_buffers.deviceFieldRealUsmRead();
    const Real *field_usm_i = value_buffers.deviceFieldImagUsmRead();
    const bool field_reads_from_usm = field_usm_r != nullptr && field_usm_i != nullptr;

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real delta_real = f_fr[index_j] - value_i_real;
                        const Real delta_imag = f_fi[index_j] - value_i_imag;
                        sum_x_real += delta_real * gwx_dev[edge];
                        sum_x_imag += delta_imag * gwx_dev[edge];
                        sum_y_real += delta_real * gwy_dev[edge];
                        sum_y_imag += delta_imag * gwy_dev[edge];
                        sum_z_real += delta_real * gwz_dev[edge];
                        sum_z_imag += delta_imag * gwz_dev[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real delta_real = f_real[index_j] - value_i_real;
                        const Real delta_imag = f_imag[index_j] - value_i_imag;
                        sum_x_real += delta_real * gwx_dev[edge];
                        sum_x_imag += delta_imag * gwx_dev[edge];
                        sum_y_real += delta_real * gwy_dev[edge];
                        sum_y_imag += delta_imag * gwy_dev[edge];
                        sum_z_real += delta_real * gwz_dev[edge];
                        sum_z_imag += delta_imag * gwz_dev[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
    }
    else
    {
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real delta_real = f_fr[index_j] - value_i_real;
                        const Real delta_imag = f_fi[index_j] - value_i_imag;
                        sum_x_real += delta_real * gwx[edge];
                        sum_x_imag += delta_imag * gwx[edge];
                        sum_y_real += delta_real * gwy[edge];
                        sum_y_imag += delta_imag * gwy[edge];
                        sum_z_real += delta_real * gwz[edge];
                        sum_z_imag += delta_imag * gwz[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real delta_real = f_real[index_j] - value_i_real;
                        const Real delta_imag = f_imag[index_j] - value_i_imag;
                        sum_x_real += delta_real * gwx[edge];
                        sum_x_imag += delta_imag * gwx[edge];
                        sum_y_real += delta_real * gwy[edge];
                        sum_y_imag += delta_imag * gwy[edge];
                        sum_z_real += delta_real * gwz[edge];
                        sum_z_imag += delta_imag * gwz[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
    }

    value_buffers.readGradient(gradient);
    return gradient;
}

inline StdVec<Vec3c> applyMatrixFreeGradientSycl(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field)
{
    return applyMatrixFreeGradientSyclWithWorkspace(graph, field, matrixFreeAPhiSyclWorkspace());
}

inline StdVec<Vec3c> applyMatrixFreeHarmonicWeightedGradientSyclWithWorkspace(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
    const StdVec<Real> &edge_weight_coefficient, MatrixFreeAPhiSyclWorkspace &sycl_workspace)
{
    const size_t number_of_rows = graph.rows_.size();
    StdVec<Vec3c> gradient(field.size(), Vec3c::Zero());
    if (!graph.flat_.isValidFor(number_of_rows) || field.size() != number_of_rows ||
        edge_weight_coefficient.size() != number_of_rows)
    {
        return gradient;
    }

    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);
    value_buffers.loadFieldIfNeeded(field);
    value_buffers.loadCoefficientIfNeeded(edge_weight_coefficient);

    const Real *field_usm_r = value_buffers.deviceFieldRealUsmRead();
    const Real *field_usm_i = value_buffers.deviceFieldImagUsmRead();
    const bool field_reads_from_usm = field_usm_r != nullptr && field_usm_i != nullptr;

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real delta_real = coefficient_ij * (f_fr[index_j] - value_i_real);
                        const Real delta_imag = coefficient_ij * (f_fi[index_j] - value_i_imag);
                        sum_x_real += delta_real * gwx_dev[edge];
                        sum_x_imag += delta_imag * gwx_dev[edge];
                        sum_y_real += delta_real * gwy_dev[edge];
                        sum_y_imag += delta_imag * gwy_dev[edge];
                        sum_z_real += delta_real * gwz_dev[edge];
                        sum_z_imag += delta_imag * gwz_dev[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real delta_real = coefficient_ij * (f_real[index_j] - value_i_real);
                        const Real delta_imag = coefficient_ij * (f_imag[index_j] - value_i_imag);
                        sum_x_real += delta_real * gwx_dev[edge];
                        sum_x_imag += delta_imag * gwx_dev[edge];
                        sum_y_real += delta_real * gwy_dev[edge];
                        sum_y_imag += delta_imag * gwy_dev[edge];
                        sum_z_real += delta_real * gwz_dev[edge];
                        sum_z_imag += delta_imag * gwz_dev[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
    }
    else
    {
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real delta_real = coefficient_ij * (f_fr[index_j] - value_i_real);
                        const Real delta_imag = coefficient_ij * (f_fi[index_j] - value_i_imag);
                        sum_x_real += delta_real * gwx[edge];
                        sum_x_imag += delta_imag * gwx[edge];
                        sum_y_real += delta_real * gwy[edge];
                        sum_y_imag += delta_imag * gwy[edge];
                        sum_z_real += delta_real * gwz[edge];
                        sum_z_imag += delta_imag * gwz[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real delta_real = coefficient_ij * (f_real[index_j] - value_i_real);
                        const Real delta_imag = coefficient_ij * (f_imag[index_j] - value_i_imag);
                        sum_x_real += delta_real * gwx[edge];
                        sum_x_imag += delta_imag * gwx[edge];
                        sum_y_real += delta_real * gwy[edge];
                        sum_y_imag += delta_imag * gwy[edge];
                        sum_z_real += delta_real * gwz[edge];
                        sum_z_imag += delta_imag * gwz[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }
    }

    value_buffers.readGradient(gradient);
    return gradient;
}

inline StdVec<Vec3c> applyMatrixFreeHarmonicWeightedGradientSycl(const MatrixFreePairwiseGraph &graph,
                                                                const StdVec<Complex> &field,
                                                                const StdVec<Real> &edge_weight_coefficient)
{
    return applyMatrixFreeHarmonicWeightedGradientSyclWithWorkspace(graph, field, edge_weight_coefficient,
                                                                   matrixFreeAPhiSyclWorkspace());
}
} // namespace

inline void setMatrixFreeGradientUseSycl(bool enable)
{
    g_matrix_free_gradient_use_sycl = enable;
}

inline bool matrixFreeGradientUseSycl()
{
    return g_matrix_free_gradient_use_sycl;
}

inline void setMatrixFreeHarmonicGradientUseSycl(bool enable)
{
    g_matrix_free_harmonic_gradient_use_sycl = enable;
}

inline bool matrixFreeHarmonicGradientUseSycl()
{
    return g_matrix_free_harmonic_gradient_use_sycl;
}

/** One field upload, one kernel, one host read for both standard and harmonic weighted gradients of the same scalar field (e.g. phi). */
inline bool tryApplyPhiStandardAndHarmonicWeightedGradientsFusedSycl(const MatrixFreePairwiseGraph &graph,
                                                                      const StdVec<Complex> &phi,
                                                                      const StdVec<Real> &edge_weight_coefficient,
                                                                      StdVec<Vec3c> &standard_gradient,
                                                                      StdVec<Vec3c> &harmonic_gradient)
{
    if (!matrixFreeGradientUseSycl() || !matrixFreeHarmonicGradientUseSycl())
    {
        return false;
    }
    const size_t number_of_rows = graph.rows_.size();
    if (!graph.flat_.isValidFor(number_of_rows) || phi.size() != number_of_rows ||
        edge_weight_coefficient.size() != number_of_rows || standard_gradient.size() != number_of_rows ||
        harmonic_gradient.size() != number_of_rows)
    {
        return false;
    }

    MatrixFreeAPhiSyclWorkspace &sycl_workspace = matrixFreeAPhiSyclWorkspace();
    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);
    value_buffers.loadFieldIfNeeded(phi);
    value_buffers.loadCoefficientIfNeeded(edge_weight_coefficient);

    const Real *field_usm_r = value_buffers.deviceFieldRealUsmRead();
    const Real *field_usm_i = value_buffers.deviceFieldImagUsmRead();
    const bool field_reads_from_usm = field_usm_r != nullptr && field_usm_i != nullptr;

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto sxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto sxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxr = value_buffers.harmonic_gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxi = value_buffers.harmonic_gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyr = value_buffers.harmonic_gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyi = value_buffers.harmonic_gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzr = value_buffers.harmonic_gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzi = value_buffers.harmonic_gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real std_x_real = 0.0;
                    Real std_x_imag = 0.0;
                    Real std_y_real = 0.0;
                    Real std_y_imag = 0.0;
                    Real std_z_real = 0.0;
                    Real std_z_imag = 0.0;
                    Real har_x_real = 0.0;
                    Real har_x_imag = 0.0;
                    Real har_y_real = 0.0;
                    Real har_y_imag = 0.0;
                    Real har_z_real = 0.0;
                    Real har_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real delta_real = f_fr[index_j] - value_i_real;
                        const Real delta_imag = f_fi[index_j] - value_i_imag;
                        std_x_real += delta_real * gwx_dev[edge];
                        std_x_imag += delta_imag * gwx_dev[edge];
                        std_y_real += delta_real * gwy_dev[edge];
                        std_y_imag += delta_imag * gwy_dev[edge];
                        std_z_real += delta_real * gwz_dev[edge];
                        std_z_imag += delta_imag * gwz_dev[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real weighted_real = coefficient_ij * delta_real;
                        const Real weighted_imag = coefficient_ij * delta_imag;
                        har_x_real += weighted_real * gwx_dev[edge];
                        har_x_imag += weighted_imag * gwx_dev[edge];
                        har_y_real += weighted_real * gwy_dev[edge];
                        har_y_imag += weighted_imag * gwy_dev[edge];
                        har_z_real += weighted_real * gwz_dev[edge];
                        har_z_imag += weighted_imag * gwz_dev[edge];
                    }
                    sxr[index_i] = std_x_real;
                    sxi[index_i] = std_x_imag;
                    syr[index_i] = std_y_real;
                    syi[index_i] = std_y_imag;
                    szr[index_i] = std_z_real;
                    szi[index_i] = std_z_imag;
                    hxr[index_i] = har_x_real;
                    hxi[index_i] = har_x_imag;
                    hyr[index_i] = har_y_real;
                    hyi[index_i] = har_y_imag;
                    hzr[index_i] = har_z_real;
                    hzi[index_i] = har_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto sxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto sxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxr = value_buffers.harmonic_gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxi = value_buffers.harmonic_gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyr = value_buffers.harmonic_gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyi = value_buffers.harmonic_gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzr = value_buffers.harmonic_gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzi = value_buffers.harmonic_gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real std_x_real = 0.0;
                    Real std_x_imag = 0.0;
                    Real std_y_real = 0.0;
                    Real std_y_imag = 0.0;
                    Real std_z_real = 0.0;
                    Real std_z_imag = 0.0;
                    Real har_x_real = 0.0;
                    Real har_x_imag = 0.0;
                    Real har_y_real = 0.0;
                    Real har_y_imag = 0.0;
                    Real har_z_real = 0.0;
                    Real har_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real delta_real = f_real[index_j] - value_i_real;
                        const Real delta_imag = f_imag[index_j] - value_i_imag;
                        std_x_real += delta_real * gwx_dev[edge];
                        std_x_imag += delta_imag * gwx_dev[edge];
                        std_y_real += delta_real * gwy_dev[edge];
                        std_y_imag += delta_imag * gwy_dev[edge];
                        std_z_real += delta_real * gwz_dev[edge];
                        std_z_imag += delta_imag * gwz_dev[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real weighted_real = coefficient_ij * delta_real;
                        const Real weighted_imag = coefficient_ij * delta_imag;
                        har_x_real += weighted_real * gwx_dev[edge];
                        har_x_imag += weighted_imag * gwx_dev[edge];
                        har_y_real += weighted_real * gwy_dev[edge];
                        har_y_imag += weighted_imag * gwy_dev[edge];
                        har_z_real += weighted_real * gwz_dev[edge];
                        har_z_imag += weighted_imag * gwz_dev[edge];
                    }
                    sxr[index_i] = std_x_real;
                    sxi[index_i] = std_x_imag;
                    syr[index_i] = std_y_real;
                    syi[index_i] = std_y_imag;
                    szr[index_i] = std_z_real;
                    szi[index_i] = std_z_imag;
                    hxr[index_i] = har_x_real;
                    hxi[index_i] = har_x_imag;
                    hyr[index_i] = har_y_real;
                    hyi[index_i] = har_y_imag;
                    hzr[index_i] = har_z_real;
                    hzi[index_i] = har_z_imag;
                });
            });
        }
    }
    else
    {
        if (field_reads_from_usm)
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto sxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto sxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxr = value_buffers.harmonic_gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxi = value_buffers.harmonic_gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyr = value_buffers.harmonic_gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyi = value_buffers.harmonic_gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzr = value_buffers.harmonic_gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzi = value_buffers.harmonic_gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const f_fr = field_usm_r;
                const Real *const f_fi = field_usm_i;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_fr[index_i];
                    const Real value_i_imag = f_fi[index_i];
                    Real std_x_real = 0.0;
                    Real std_x_imag = 0.0;
                    Real std_y_real = 0.0;
                    Real std_y_imag = 0.0;
                    Real std_z_real = 0.0;
                    Real std_z_imag = 0.0;
                    Real har_x_real = 0.0;
                    Real har_x_imag = 0.0;
                    Real har_y_real = 0.0;
                    Real har_y_imag = 0.0;
                    Real har_z_real = 0.0;
                    Real har_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real delta_real = f_fr[index_j] - value_i_real;
                        const Real delta_imag = f_fi[index_j] - value_i_imag;
                        std_x_real += delta_real * gwx[edge];
                        std_x_imag += delta_imag * gwx[edge];
                        std_y_real += delta_real * gwy[edge];
                        std_y_imag += delta_imag * gwy[edge];
                        std_z_real += delta_real * gwz[edge];
                        std_z_imag += delta_imag * gwz[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real weighted_real = coefficient_ij * delta_real;
                        const Real weighted_imag = coefficient_ij * delta_imag;
                        har_x_real += weighted_real * gwx[edge];
                        har_x_imag += weighted_imag * gwx[edge];
                        har_y_real += weighted_real * gwy[edge];
                        har_y_imag += weighted_imag * gwy[edge];
                        har_z_real += weighted_real * gwz[edge];
                        har_z_imag += weighted_imag * gwz[edge];
                    }
                    sxr[index_i] = std_x_real;
                    sxi[index_i] = std_x_imag;
                    syr[index_i] = std_y_real;
                    syi[index_i] = std_y_imag;
                    szr[index_i] = std_z_real;
                    szi[index_i] = std_z_imag;
                    hxr[index_i] = har_x_real;
                    hxi[index_i] = har_x_imag;
                    hyr[index_i] = har_y_real;
                    hyi[index_i] = har_y_imag;
                    hzr[index_i] = har_z_real;
                    hzi[index_i] = har_z_imag;
                });
            });
        }
        else
        {
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto row_offsets = graph_buffers.row_offsets_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto column_indices = graph_buffers.column_indices_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_real = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto f_imag = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwx = graph_buffers.gradient_weight_x_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwy = graph_buffers.gradient_weight_y_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto gwz = graph_buffers.gradient_weight_z_buffer_->get_access<sycl::access::mode::read>(cgh);
                auto sxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto sxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto syi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto szi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxr = value_buffers.harmonic_gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hxi = value_buffers.harmonic_gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyr = value_buffers.harmonic_gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hyi = value_buffers.harmonic_gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzr = value_buffers.harmonic_gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto hzi = value_buffers.harmonic_gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real value_i_real = f_real[index_i];
                    const Real value_i_imag = f_imag[index_i];
                    Real std_x_real = 0.0;
                    Real std_x_imag = 0.0;
                    Real std_y_real = 0.0;
                    Real std_y_imag = 0.0;
                    Real std_z_real = 0.0;
                    Real std_z_imag = 0.0;
                    Real har_x_real = 0.0;
                    Real har_x_imag = 0.0;
                    Real har_y_real = 0.0;
                    Real har_y_imag = 0.0;
                    Real har_z_real = 0.0;
                    Real har_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real delta_real = f_real[index_j] - value_i_real;
                        const Real delta_imag = f_imag[index_j] - value_i_imag;
                        std_x_real += delta_real * gwx[edge];
                        std_x_imag += delta_imag * gwx[edge];
                        std_y_real += delta_real * gwy[edge];
                        std_y_imag += delta_imag * gwy[edge];
                        std_z_real += delta_real * gwz[edge];
                        std_z_imag += delta_imag * gwz[edge];
                        const Real coefficient_ij =
                            2.0 * coeff[index_i] * coeff[index_j] / (coeff[index_i] + coeff[index_j] + TinyReal);
                        const Real weighted_real = coefficient_ij * delta_real;
                        const Real weighted_imag = coefficient_ij * delta_imag;
                        har_x_real += weighted_real * gwx[edge];
                        har_x_imag += weighted_imag * gwx[edge];
                        har_y_real += weighted_real * gwy[edge];
                        har_y_imag += weighted_imag * gwy[edge];
                        har_z_real += weighted_real * gwz[edge];
                        har_z_imag += weighted_imag * gwz[edge];
                    }
                    sxr[index_i] = std_x_real;
                    sxi[index_i] = std_x_imag;
                    syr[index_i] = std_y_real;
                    syi[index_i] = std_y_imag;
                    szr[index_i] = std_z_real;
                    szi[index_i] = std_z_imag;
                    hxr[index_i] = har_x_real;
                    hxi[index_i] = har_x_imag;
                    hyr[index_i] = har_y_real;
                    hyi[index_i] = har_y_imag;
                    hzr[index_i] = har_z_real;
                    hzi[index_i] = har_z_imag;
                });
            });
        }
    }

    value_buffers.readFusedStandardAndHarmonicGradients(standard_gradient, harmonic_gradient);
    return true;
}

/** Writes standard-gradient pairwise scalar divergence of the six A components already uploaded into
 *  `a_stage` into `a_stage` div USM slots. Caller must have `a_stage.ready()`, `graph_buffers.loadIfNeeded`, and
 *  `a_stage.uploadAFromHost(...)`. */
inline void matrixFreeSubmitStandardDivergenceThreeComponentsToStagingUsm(MatrixFreeSyclFlatGraphBuffers &graph_buffers,
                                                                        MatrixFreeSyclVectorPotentialDivergenceStaging &a_stage,
                                                                        size_t number_of_rows)
{
    const Real *const axr_usm = a_stage.axrRead();
    const Real *const axi_usm = a_stage.axiRead();
    const Real *const ayr_usm = a_stage.ayrRead();
    const Real *const ayi_usm = a_stage.ayiRead();
    const Real *const azr_usm = a_stage.azrRead();
    const Real *const azi_usm = a_stage.aziRead();
    Real *const div_write_r = a_stage.mutableDivReal();
    Real *const div_write_i = a_stage.mutableDivImag();

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            const Real *const p_axr = axr_usm;
            const Real *const p_axi = axi_usm;
            const Real *const p_ayr = ayr_usm;
            const Real *const p_ayi = ayi_usm;
            const Real *const p_azr = azr_usm;
            const Real *const p_azi = azi_usm;
            Real *const p_divr = div_write_r;
            Real *const p_divi = div_write_i;

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                Real div_real = 0.0;
                Real div_imag = 0.0;
                for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                {
                    const size_t index_j = col_dev[edge];
                    div_real += gwx_dev[edge] * (p_axr[index_j] - p_axr[index_i]) +
                                gwy_dev[edge] * (p_ayr[index_j] - p_ayr[index_i]) +
                                gwz_dev[edge] * (p_azr[index_j] - p_azr[index_i]);
                    div_imag += gwx_dev[edge] * (p_axi[index_j] - p_axi[index_i]) +
                                gwy_dev[edge] * (p_ayi[index_j] - p_ayi[index_i]) +
                                gwz_dev[edge] * (p_azi[index_j] - p_azi[index_i]);
                }
                p_divr[index_i] = div_real;
                p_divi[index_i] = div_imag;
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
            const Real *const p_axr = axr_usm;
            const Real *const p_axi = axi_usm;
            const Real *const p_ayr = ayr_usm;
            const Real *const p_ayi = ayi_usm;
            const Real *const p_azr = azr_usm;
            const Real *const p_azi = azi_usm;
            Real *const p_divr = div_write_r;
            Real *const p_divi = div_write_i;

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                Real div_real = 0.0;
                Real div_imag = 0.0;
                for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                {
                    const size_t index_j = column_indices[edge];
                    div_real += gwx[edge] * (p_axr[index_j] - p_axr[index_i]) +
                                gwy[edge] * (p_ayr[index_j] - p_ayr[index_i]) +
                                gwz[edge] * (p_azr[index_j] - p_azr[index_i]);
                    div_imag += gwx[edge] * (p_axi[index_j] - p_axi[index_i]) +
                                gwy[edge] * (p_ayi[index_j] - p_ayi[index_i]) +
                                gwz[edge] * (p_azi[index_j] - p_azi[index_i]);
                }
                p_divr[index_i] = div_real;
                p_divi[index_i] = div_imag;
            });
        });
    }
}

/** Scalar divergence of three complex component fields with **standard** pairwise gradient weights (same as
 *  `computeGradientDivergenceOfVectorField` SYCL buffer path). When field-USM staging is ready, avoids wrapping
 *  six component arrays in temporary `sycl::buffer`. */
inline bool tryMatrixFreeGradientScalarDivergenceThreeComponentsWithStagingSycl(const MatrixFreePairwiseGraph &graph,
                                                                               MatrixFreeAPhiSyclWorkspace &sycl_workspace,
                                                                               const StdVec<Real> &fxr, const StdVec<Real> &fxi,
                                                                               const StdVec<Real> &fyr, const StdVec<Real> &fyi,
                                                                               const StdVec<Real> &fzr, const StdVec<Real> &fzi,
                                                                               StdVec<Complex> &divergence_out)
{
    if (!matrixFreeGradientUseSycl())
    {
        return false;
    }
    const size_t number_of_rows = graph.rows_.size();
    if (!graph.flat_.isValidFor(number_of_rows) || fxr.size() != number_of_rows || divergence_out.size() != number_of_rows)
    {
        return false;
    }

    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclVectorPotentialDivergenceStaging &a_stage = sycl_workspace.vector_potential_divergence_staging_;
    a_stage.ensureForRowCount(number_of_rows);
    if (!a_stage.ready())
    {
        return false;
    }
    a_stage.uploadAFromHost(fxr, fxi, fyr, fyi, fzr, fzi);
    matrixFreeSubmitStandardDivergenceThreeComponentsToStagingUsm(graph_buffers, a_stage, number_of_rows);

    a_stage.readDivToComplexHost(divergence_out);
    return true;
}

/** Scalar divergence with **harmonic** edge weights (same as `computeHarmonicDivergenceOfVectorField` SYCL buffer path).
 *  Uses `scalar_value_buffers_` for `coefficient_buffer_` reads; six vector components use vector-potential staging USM. */
inline bool tryMatrixFreeHarmonicScalarDivergenceThreeComponentsWithStagingSycl(
    const MatrixFreePairwiseGraph &graph, MatrixFreeAPhiSyclWorkspace &sycl_workspace, const StdVec<Real> &fxr,
    const StdVec<Real> &fxi, const StdVec<Real> &fyr, const StdVec<Real> &fyi, const StdVec<Real> &fzr,
    const StdVec<Real> &fzi, const StdVec<Real> &edge_weight_coefficient, StdVec<Complex> &divergence_out)
{
    if (!matrixFreeHarmonicGradientUseSycl())
    {
        return false;
    }
    const size_t number_of_rows = graph.rows_.size();
    if (!graph.flat_.isValidFor(number_of_rows) || fxr.size() != number_of_rows ||
        edge_weight_coefficient.size() != number_of_rows || divergence_out.size() != number_of_rows)
    {
        return false;
    }

    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);
    value_buffers.loadCoefficientIfNeeded(edge_weight_coefficient);

    MatrixFreeSyclVectorPotentialDivergenceStaging &a_stage = sycl_workspace.vector_potential_divergence_staging_;
    a_stage.ensureForRowCount(number_of_rows);
    if (!a_stage.ready())
    {
        return false;
    }
    a_stage.uploadAFromHost(fxr, fxi, fyr, fyi, fzr, fzi);
    const Real *const axr_usm = a_stage.axrRead();
    const Real *const axi_usm = a_stage.axiRead();
    const Real *const ayr_usm = a_stage.ayrRead();
    const Real *const ayi_usm = a_stage.ayiRead();
    const Real *const azr_usm = a_stage.azrRead();
    const Real *const azi_usm = a_stage.aziRead();
    Real *const div_write_r = a_stage.mutableDivReal();
    Real *const div_write_i = a_stage.mutableDivImag();

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            auto coeff = value_buffers.coefficient_buffer_->get_access<sycl::access::mode::read>(cgh);
            const Real *const p_axr = axr_usm;
            const Real *const p_axi = axi_usm;
            const Real *const p_ayr = ayr_usm;
            const Real *const p_ayi = ayi_usm;
            const Real *const p_azr = azr_usm;
            const Real *const p_azi = azi_usm;
            Real *const p_divr = div_write_r;
            Real *const p_divi = div_write_i;

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
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
                    div_real += weighted_x * (p_axr[index_j] - p_axr[index_i]) +
                                weighted_y * (p_ayr[index_j] - p_ayr[index_i]) +
                                weighted_z * (p_azr[index_j] - p_azr[index_i]);
                    div_imag += weighted_x * (p_axi[index_j] - p_axi[index_i]) +
                                weighted_y * (p_ayi[index_j] - p_ayi[index_i]) +
                                weighted_z * (p_azi[index_j] - p_azi[index_i]);
                }
                p_divr[index_i] = div_real;
                p_divi[index_i] = div_imag;
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
            const Real *const p_axr = axr_usm;
            const Real *const p_axi = axi_usm;
            const Real *const p_ayr = ayr_usm;
            const Real *const p_ayi = ayi_usm;
            const Real *const p_azr = azr_usm;
            const Real *const p_azi = azi_usm;
            Real *const p_divr = div_write_r;
            Real *const p_divi = div_write_i;

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
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
                    div_real += weighted_x * (p_axr[index_j] - p_axr[index_i]) +
                                weighted_y * (p_ayr[index_j] - p_ayr[index_i]) +
                                weighted_z * (p_azr[index_j] - p_azr[index_i]);
                    div_imag += weighted_x * (p_axi[index_j] - p_axi[index_i]) +
                                weighted_y * (p_ayi[index_j] - p_ayi[index_i]) +
                                weighted_z * (p_azi[index_j] - p_azi[index_i]);
                }
                p_divr[index_i] = div_real;
                p_divi[index_i] = div_imag;
            });
        });
    }

    a_stage.readDivToComplexHost(divergence_out);
    return true;
}

/** div(A) on device with the same weights as computeGradientDivergenceOfVectorField, then grad(div A) like
 *  applyMatrixFreeGradient on that scalar field. Avoids a host round-trip and re-upload of div(A) when both are needed.
 *  If divergence_out is null, div is not copied to host (gauge penalty path that only needs grad(div A)). */
inline bool tryMatrixFreeDivergenceAndGradientOfDivergenceOfVectorPotentialSycl(
    const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
    const StdVec<Complex> &field_z, StdVec<Complex> *divergence_out, StdVec<Vec3c> &gradient_of_divergence)
{
    if (!matrixFreeGradientUseSycl())
    {
        return false;
    }
    const size_t number_of_rows = graph.rows_.size();
    if (!graph.flat_.isValidFor(number_of_rows) || field_x.size() != number_of_rows ||
        field_y.size() != number_of_rows || field_z.size() != number_of_rows ||
        gradient_of_divergence.size() != number_of_rows)
    {
        return false;
    }
    if (divergence_out != nullptr && divergence_out->size() != number_of_rows)
    {
        return false;
    }

    StdVec<Real> fxr(number_of_rows, 0.0);
    StdVec<Real> fxi(number_of_rows, 0.0);
    StdVec<Real> fyr(number_of_rows, 0.0);
    StdVec<Real> fyi(number_of_rows, 0.0);
    StdVec<Real> fzr(number_of_rows, 0.0);
    StdVec<Real> fzi(number_of_rows, 0.0);
    for (size_t i = 0; i != number_of_rows; ++i)
    {
        fxr[i] = field_x[i].real();
        fxi[i] = field_x[i].imag();
        fyr[i] = field_y[i].real();
        fyi[i] = field_y[i].imag();
        fzr[i] = field_z[i].real();
        fzi[i] = field_z[i].imag();
    }

    MatrixFreeAPhiSyclWorkspace &sycl_workspace = matrixFreeAPhiSyclWorkspace();
    MatrixFreeSyclFlatGraphBuffers &graph_buffers = sycl_workspace.graph_buffers_;
    graph_buffers.loadIfNeeded(graph.flat_, number_of_rows);
    MatrixFreeSyclScalarValueBuffers &value_buffers = sycl_workspace.scalar_value_buffers_;
    value_buffers.resize(number_of_rows);

    MatrixFreeSyclVectorPotentialDivergenceStaging &a_stage = sycl_workspace.vector_potential_divergence_staging_;
    a_stage.ensureForRowCount(number_of_rows);
    if (a_stage.ready())
    {
        a_stage.uploadAFromHost(fxr, fxi, fyr, fyi, fzr, fzi);
        matrixFreeSubmitStandardDivergenceThreeComponentsToStagingUsm(graph_buffers, a_stage, number_of_rows);

        if (divergence_out != nullptr)
        {
            a_stage.readDivToComplexHost(*divergence_out);
        }

        const Real *const dr_usm = a_stage.divRealRead();
        const Real *const di_usm = a_stage.divImagRead();

        if (graph_buffers.deviceTopologyUsmReady())
        {
            const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
            const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
            const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
            const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
            const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
            matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const dr = dr_usm;
                const Real *const di = di_usm;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real d_i_re = dr[index_i];
                    const Real d_i_im = di[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                    {
                        const size_t index_j = col_dev[edge];
                        const Real delta_real = dr[index_j] - d_i_re;
                        const Real delta_imag = di[index_j] - d_i_im;
                        sum_x_real += delta_real * gwx_dev[edge];
                        sum_x_imag += delta_imag * gwx_dev[edge];
                        sum_y_real += delta_real * gwy_dev[edge];
                        sum_y_imag += delta_imag * gwy_dev[edge];
                        sum_z_real += delta_real * gwz_dev[edge];
                        sum_z_imag += delta_imag * gwz_dev[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
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
                auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
                auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
                const Real *const dr = dr_usm;
                const Real *const di = di_usm;

                cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                    const size_t index_i = id[0];
                    const Real d_i_re = dr[index_i];
                    const Real d_i_im = di[index_i];
                    Real sum_x_real = 0.0;
                    Real sum_x_imag = 0.0;
                    Real sum_y_real = 0.0;
                    Real sum_y_imag = 0.0;
                    Real sum_z_real = 0.0;
                    Real sum_z_imag = 0.0;
                    for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                    {
                        const size_t index_j = column_indices[edge];
                        const Real delta_real = dr[index_j] - d_i_re;
                        const Real delta_imag = di[index_j] - d_i_im;
                        sum_x_real += delta_real * gwx[edge];
                        sum_x_imag += delta_imag * gwx[edge];
                        sum_y_real += delta_real * gwy[edge];
                        sum_y_imag += delta_imag * gwy[edge];
                        sum_z_real += delta_real * gwz[edge];
                        sum_z_imag += delta_imag * gwz[edge];
                    }
                    gxr[index_i] = sum_x_real;
                    gxi[index_i] = sum_x_imag;
                    gyr[index_i] = sum_y_real;
                    gyi[index_i] = sum_y_imag;
                    gzr[index_i] = sum_z_real;
                    gzi[index_i] = sum_z_imag;
                });
            });
        }

        value_buffers.readGradient(gradient_of_divergence);
        return true;
    }

    sycl::buffer<Real, 1> buf_fxr{fxr.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> buf_fxi{fxi.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> buf_fyr{fyr.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> buf_fyi{fyi.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> buf_fzr{fzr.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> buf_fzi{fzi.data(), sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> div_real_buffer{sycl::range<1>(number_of_rows)};
    sycl::buffer<Real, 1> div_imag_buffer{sycl::range<1>(number_of_rows)};

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            auto axr = buf_fxr.get_access<sycl::access::mode::read>(cgh);
            auto axi = buf_fxi.get_access<sycl::access::mode::read>(cgh);
            auto ayr = buf_fyr.get_access<sycl::access::mode::read>(cgh);
            auto ayi = buf_fyi.get_access<sycl::access::mode::read>(cgh);
            auto azr = buf_fzr.get_access<sycl::access::mode::read>(cgh);
            auto azi = buf_fzi.get_access<sycl::access::mode::read>(cgh);
            auto divr = div_real_buffer.get_access<sycl::access::mode::write>(cgh);
            auto divi = div_imag_buffer.get_access<sycl::access::mode::write>(cgh);

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                Real div_real = 0.0;
                Real div_imag = 0.0;
                for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                {
                    const size_t index_j = col_dev[edge];
                    div_real += gwx_dev[edge] * (axr[index_j] - axr[index_i]) +
                                gwy_dev[edge] * (ayr[index_j] - ayr[index_i]) +
                                gwz_dev[edge] * (azr[index_j] - azr[index_i]);
                    div_imag += gwx_dev[edge] * (axi[index_j] - axi[index_i]) +
                                gwy_dev[edge] * (ayi[index_j] - ayi[index_i]) +
                                gwz_dev[edge] * (azi[index_j] - azi[index_i]);
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
            auto axr = buf_fxr.get_access<sycl::access::mode::read>(cgh);
            auto axi = buf_fxi.get_access<sycl::access::mode::read>(cgh);
            auto ayr = buf_fyr.get_access<sycl::access::mode::read>(cgh);
            auto ayi = buf_fyi.get_access<sycl::access::mode::read>(cgh);
            auto azr = buf_fzr.get_access<sycl::access::mode::read>(cgh);
            auto azi = buf_fzi.get_access<sycl::access::mode::read>(cgh);
            auto divr = div_real_buffer.get_access<sycl::access::mode::write>(cgh);
            auto divi = div_imag_buffer.get_access<sycl::access::mode::write>(cgh);

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                Real div_real = 0.0;
                Real div_imag = 0.0;
                for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                {
                    const size_t index_j = column_indices[edge];
                    div_real += gwx[edge] * (axr[index_j] - axr[index_i]) + gwy[edge] * (ayr[index_j] - ayr[index_i]) +
                                gwz[edge] * (azr[index_j] - azr[index_i]);
                    div_imag += gwx[edge] * (axi[index_j] - axi[index_i]) + gwy[edge] * (ayi[index_j] - ayi[index_i]) +
                                gwz[edge] * (azi[index_j] - azi[index_i]);
                }
                divr[index_i] = div_real;
                divi[index_i] = div_imag;
            });
        });
    }

    if (divergence_out != nullptr)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::read> divr(div_real_buffer);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> divi(div_imag_buffer);
        for (size_t i = 0; i != number_of_rows; ++i)
        {
            (*divergence_out)[i] = Complex(divr[i], divi[i]);
        }
    }

    if (graph_buffers.deviceTopologyUsmReady())
    {
        const size_t *row_dev = graph_buffers.deviceTopologyRowOffsets();
        const size_t *col_dev = graph_buffers.deviceTopologyColumnIndices();
        const Real *gwx_dev = graph_buffers.deviceTopologyGradientWeightX();
        const Real *gwy_dev = graph_buffers.deviceTopologyGradientWeightY();
        const Real *gwz_dev = graph_buffers.deviceTopologyGradientWeightZ();
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            auto dr = div_real_buffer.get_access<sycl::access::mode::read>(cgh);
            auto di = div_imag_buffer.get_access<sycl::access::mode::read>(cgh);
            auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                const Real d_i_re = dr[index_i];
                const Real d_i_im = di[index_i];
                Real sum_x_real = 0.0;
                Real sum_x_imag = 0.0;
                Real sum_y_real = 0.0;
                Real sum_y_imag = 0.0;
                Real sum_z_real = 0.0;
                Real sum_z_imag = 0.0;
                for (size_t edge = row_dev[index_i]; edge != row_dev[index_i + 1]; ++edge)
                {
                    const size_t index_j = col_dev[edge];
                    const Real delta_real = dr[index_j] - d_i_re;
                    const Real delta_imag = di[index_j] - d_i_im;
                    sum_x_real += delta_real * gwx_dev[edge];
                    sum_x_imag += delta_imag * gwx_dev[edge];
                    sum_y_real += delta_real * gwy_dev[edge];
                    sum_y_imag += delta_imag * gwy_dev[edge];
                    sum_z_real += delta_real * gwz_dev[edge];
                    sum_z_imag += delta_imag * gwz_dev[edge];
                }
                gxr[index_i] = sum_x_real;
                gxi[index_i] = sum_x_imag;
                gyr[index_i] = sum_y_real;
                gyi[index_i] = sum_y_imag;
                gzr[index_i] = sum_z_real;
                gzi[index_i] = sum_z_imag;
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
            auto dr = div_real_buffer.get_access<sycl::access::mode::read>(cgh);
            auto di = div_imag_buffer.get_access<sycl::access::mode::read>(cgh);
            auto gxr = value_buffers.gradient_x_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gxi = value_buffers.gradient_x_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gyr = value_buffers.gradient_y_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gyi = value_buffers.gradient_y_imag_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gzr = value_buffers.gradient_z_real_buffer_->get_access<sycl::access::mode::write>(cgh);
            auto gzi = value_buffers.gradient_z_imag_buffer_->get_access<sycl::access::mode::write>(cgh);

            cgh.parallel_for(sycl::range<1>(number_of_rows), [=](sycl::id<1> id) {
                const size_t index_i = id[0];
                const Real d_i_re = dr[index_i];
                const Real d_i_im = di[index_i];
                Real sum_x_real = 0.0;
                Real sum_x_imag = 0.0;
                Real sum_y_real = 0.0;
                Real sum_y_imag = 0.0;
                Real sum_z_real = 0.0;
                Real sum_z_imag = 0.0;
                for (size_t edge = row_offsets[index_i]; edge != row_offsets[index_i + 1]; ++edge)
                {
                    const size_t index_j = column_indices[edge];
                    const Real delta_real = dr[index_j] - d_i_re;
                    const Real delta_imag = di[index_j] - d_i_im;
                    sum_x_real += delta_real * gwx[edge];
                    sum_x_imag += delta_imag * gwx[edge];
                    sum_y_real += delta_real * gwy[edge];
                    sum_y_imag += delta_imag * gwy[edge];
                    sum_z_real += delta_real * gwz[edge];
                    sum_z_imag += delta_imag * gwz[edge];
                }
                gxr[index_i] = sum_x_real;
                gxi[index_i] = sum_x_imag;
                gyr[index_i] = sum_y_real;
                gyi[index_i] = sum_y_imag;
                gzr[index_i] = sum_z_real;
                gzi[index_i] = sum_z_imag;
            });
        });
    }

    value_buffers.readGradient(gradient_of_divergence);
    return true;
}
#endif

inline StdVec<Vec3c> applyMatrixFreeGradient(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field)
{
    StdVec<Vec3c> gradient(field.size(), Vec3c::Zero());
#if SPHINXSYS_USE_SYCL
    if (matrixFreeGradientUseSycl())
    {
        return applyMatrixFreeGradientSycl(graph, field);
    }
#endif
    if (graph.flat_.isValidFor(graph.rows_.size()))
    {
        for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
        {
            const Complex value_i = field[index_i];
            for (size_t edge = graph.flat_.row_offsets_[index_i]; edge != graph.flat_.row_offsets_[index_i + 1]; ++edge)
            {
                const Complex delta = field[graph.flat_.column_indices_[edge]] - value_i;
                for (int axis = 0; axis != Dimensions; ++axis)
                {
                    gradient[index_i][axis] += delta * graph.flat_.gradient_weights_[edge][axis];
                }
            }
        }
        return gradient;
    }

    for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
    {
        const Complex value_i = field[index_i];
        for (const MatrixFreePairwiseNeighborEntry &neighbor : graph.rows_[index_i].neighbors_)
        {
            const Complex delta = field[neighbor.index_j_] - value_i;
            for (int axis = 0; axis != Dimensions; ++axis)
            {
                gradient[index_i][axis] += delta * neighbor.gradient_weight_[axis];
            }
        }
    }
    return gradient;
}

inline StdVec<Vec3c> applyMatrixFreeHarmonicWeightedGradient(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                             const StdVec<Real> &edge_weight_coefficient)
{
    StdVec<Vec3c> gradient(field.size(), Vec3c::Zero());
#if SPHINXSYS_USE_SYCL
    if (matrixFreeHarmonicGradientUseSycl())
    {
        return applyMatrixFreeHarmonicWeightedGradientSycl(graph, field, edge_weight_coefficient);
    }
#endif
    if (graph.flat_.isValidFor(graph.rows_.size()))
    {
        for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
        {
            const Complex value_i = field[index_i];
            for (size_t edge = graph.flat_.row_offsets_[index_i]; edge != graph.flat_.row_offsets_[index_i + 1]; ++edge)
            {
                const size_t index_j = graph.flat_.column_indices_[edge];
                const Real coefficient_ij = harmonicMean(edge_weight_coefficient[index_i], edge_weight_coefficient[index_j]);
                const Complex delta = coefficient_ij * (field[index_j] - value_i);
                for (int axis = 0; axis != Dimensions; ++axis)
                {
                    gradient[index_i][axis] += delta * graph.flat_.gradient_weights_[edge][axis];
                }
            }
        }
        return gradient;
    }

    for (size_t index_i = 0; index_i != graph.rows_.size(); ++index_i)
    {
        const Complex value_i = field[index_i];
        for (const MatrixFreePairwiseNeighborEntry &neighbor : graph.rows_[index_i].neighbors_)
        {
            const size_t index_j = neighbor.index_j_;
            const Real coefficient_ij = harmonicMean(edge_weight_coefficient[index_i], edge_weight_coefficient[index_j]);
            const Complex delta = coefficient_ij * (field[index_j] - value_i);
            for (int axis = 0; axis != Dimensions; ++axis)
            {
                gradient[index_i][axis] += delta * neighbor.gradient_weight_[axis];
            }
        }
    }
    return gradient;
}

inline void accumulateScalarDivergenceOfGradientResidualsFromGraph(const MatrixFreePairwiseGraph &graph,
                                                                   const StdVec<Complex> &field,
                                                                   ScalarComplexHelmholtzResiduals &residuals)
{
    const size_t number_of_particles = field.size();
    const StdVec<Vec3c> gradient = applyMatrixFreeGradient(graph, field);

    StdVec<Complex> gradient_x(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> gradient_y(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> gradient_z(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        gradient_x[i] = gradient[i][0];
        gradient_y[i] = gradient[i][1];
        gradient_z[i] = gradient[i][2];
    }

    const StdVec<Vec3c> grad_gradient_x = applyMatrixFreeGradient(graph, gradient_x);
    const StdVec<Vec3c> grad_gradient_y = applyMatrixFreeGradient(graph, gradient_y);
    const StdVec<Vec3c> grad_gradient_z = applyMatrixFreeGradient(graph, gradient_z);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.laplace_term_[i] += grad_gradient_x[i][0] + grad_gradient_y[i][1] + grad_gradient_z[i][2];
    }

    ScalarComplexHelmholtzResiduals diagonal_proxy;
    diagonal_proxy.resize(number_of_particles);
    diagonal_proxy.clear();
    accumulateScalarLaplaceResidualsFromClearedGraph(graph, field, StdVec<Real>(number_of_particles, 1.0),
                                                    diagonal_proxy);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.diagonal_scale_[i] += diagonal_proxy.diagonal_scale_[i];
    }
}

inline StdVec<Complex> applyScalarDivergenceOfGradientFromGraph(const MatrixFreePairwiseGraph &graph,
                                                                const StdVec<Complex> &field)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(field.size());
    residuals.clear();
    accumulateScalarDivergenceOfGradientResidualsFromGraph(graph, field, residuals);
    return residuals.laplace_term_;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_HPP
