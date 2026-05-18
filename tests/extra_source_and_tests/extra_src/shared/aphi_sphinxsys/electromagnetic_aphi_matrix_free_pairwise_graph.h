#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct MatrixFreePairwiseNeighborEntry
{
    size_t index_j_;
    Vecd gradient_weight_;
    Real diffusion_weight_;
    bool is_contact_;

    MatrixFreePairwiseNeighborEntry(size_t index_j, const Vecd &gradient_weight, Real diffusion_weight, bool is_contact);
};

struct MatrixFreePairwiseRow
{
    StdVec<MatrixFreePairwiseNeighborEntry> neighbors_;
};

struct MatrixFreeFlatPairwiseGraph
{
    StdVec<size_t> row_offsets_;
    StdVec<size_t> column_indices_;
    StdVec<Vecd> gradient_weights_;
    StdVec<Real> diffusion_weights_;
    StdVec<uint8_t> contact_flags_;

    void clear();
    bool isValidFor(size_t number_of_rows) const;
};

struct MatrixFreePairwiseGraph
{
    StdVec<MatrixFreePairwiseRow> rows_;
    MatrixFreeFlatPairwiseGraph flat_;
};

MatrixFreePairwiseGraph buildMatrixFreePairwiseGraph(const MatrixFreeAPhiDiscreteView &discrete_view,
                                                     const MatrixFreeAPhiParameters &parameters);

void refreshMatrixFreeFlatPairwiseGraph(MatrixFreePairwiseGraph &graph);

#if SPHINXSYS_USE_SYCL
void prepareMatrixFreeAPhiSyclGraphWorkspace(const MatrixFreePairwiseGraph &graph);
size_t matrixFreeAPhiSyclGraphWorkspaceUploadCount();
size_t matrixFreeAPhiSyclValueFieldUploadCount();
size_t matrixFreeAPhiSyclValueCoefficientUploadCount();
size_t matrixFreeAPhiSyclValueLaplaceUploadCount();
size_t matrixFreeAPhiSyclValueGradientDownloadCount();
size_t matrixFreeAPhiSyclValueLaplaceDownloadCount();
size_t matrixFreeAPhiSyclValueFieldUploadSkipCount();
void invalidateMatrixFreeAPhiSyclValueFieldCache();
/** Upload flat graph CSR-like buffers and pre-size scalar SYCL workspace rows (host data still staged per kernel). */
void prepareMatrixFreeAPhiSyclDeviceResources(const MatrixFreePairwiseGraph &graph, size_t scalar_row_count);
/** True when `EM_APHI_MATRIX_FREE_SYCL_GRAPH_TOPOLOGY_USM=1`: flat graph CSR + per-edge gradient/diffusion weights mirrored to device USM for all graph-based SYCL paths here (Laplace, gradients, fused phi, div(A)/grad(div A), vector divergences in residuals, Helmholtz residual kernel); field components remain in SYCL buffers. */
bool matrixFreeSyclGraphTopologyUsesDeviceUsm();
/** When `EM_APHI_MATRIX_FREE_SYCL_FIELD_VALUES_USM=1`, device USM mirrors for scalar field real/imag: Jacobi (and related) may write USM; Laplace / gradient / fused phi / graph Helmholtz residual kernels may read field from USM when allocated. SYCL field buffers are refreshed when syncing to host or when code paths require buffer-backed reads. */
bool matrixFreeSyclFieldValuesUseDeviceUsm();
#else
inline void prepareMatrixFreeAPhiSyclGraphWorkspace(const MatrixFreePairwiseGraph &) {}
inline size_t matrixFreeAPhiSyclGraphWorkspaceUploadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueFieldUploadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueCoefficientUploadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueLaplaceUploadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueGradientDownloadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueLaplaceDownloadCount()
{
    return 0;
}
inline size_t matrixFreeAPhiSyclValueFieldUploadSkipCount()
{
    return 0;
}
inline void invalidateMatrixFreeAPhiSyclValueFieldCache() {}
inline void prepareMatrixFreeAPhiSyclDeviceResources(const MatrixFreePairwiseGraph &, size_t) {}
inline bool matrixFreeSyclGraphTopologyUsesDeviceUsm()
{
    return false;
}
inline bool matrixFreeSyclFieldValuesUseDeviceUsm()
{
    return false;
}
#endif

void accumulateScalarLaplaceResidualsFromGraph(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                               const StdVec<Real> &diffusion_coefficient,
                                               ScalarComplexHelmholtzResiduals &residuals);
void accumulateScalarLaplaceResidualsFromClearedGraph(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                      const StdVec<Real> &diffusion_coefficient,
                                                      ScalarComplexHelmholtzResiduals &residuals);

#if SPHINXSYS_USE_SYCL
void setMatrixFreeLaplaceResidualUseSycl(bool enable);
bool matrixFreeLaplaceResidualUseSycl();
#else
inline void setMatrixFreeLaplaceResidualUseSycl(bool) {}
inline bool matrixFreeLaplaceResidualUseSycl()
{
    return false;
}
#endif

StdVec<Complex> applyScalarNegativeLaplaceFromGraph(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                    const StdVec<Real> &diffusion_coefficient);

StdVec<Vec3c> applyMatrixFreeGradient(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field);

#if SPHINXSYS_USE_SYCL
void setMatrixFreeGradientUseSycl(bool enable);
bool matrixFreeGradientUseSycl();
#else
inline void setMatrixFreeGradientUseSycl(bool) {}
inline bool matrixFreeGradientUseSycl()
{
    return false;
}
#endif

StdVec<Vec3c> applyMatrixFreeHarmonicWeightedGradient(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &field,
                                                      const StdVec<Real> &edge_weight_coefficient);

#if SPHINXSYS_USE_SYCL
void setMatrixFreeHarmonicGradientUseSycl(bool enable);
bool matrixFreeHarmonicGradientUseSycl();
#else
inline void setMatrixFreeHarmonicGradientUseSycl(bool) {}
inline bool matrixFreeHarmonicGradientUseSycl()
{
    return false;
}
#endif

void accumulateScalarDivergenceOfGradientResidualsFromGraph(const MatrixFreePairwiseGraph &graph,
                                                            const StdVec<Complex> &field,
                                                            ScalarComplexHelmholtzResiduals &residuals);

StdVec<Complex> applyScalarDivergenceOfGradientFromGraph(const MatrixFreePairwiseGraph &graph,
                                                         const StdVec<Complex> &field);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.hpp"

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_PAIRWISE_GRAPH_H
