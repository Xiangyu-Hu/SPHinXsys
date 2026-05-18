#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct MatrixFreeAPhiSphNativeContext;

enum class MatrixFreeAPhiDiscreteOperatorSemanticsKind
{
    LegacyGraphCompatibility,
    SphNativeParticleKernel
};

struct MatrixFreeAPhiResidualOperatorSemantics
{
    MatrixFreeAPhiDiscreteOperatorSemanticsKind particle_kernel;
    MatrixFreeAPhiDiscreteOperatorSemanticsKind legacy_solver_gradient;

    bool usesParticleKernel() const
    {
        return particle_kernel == MatrixFreeAPhiDiscreteOperatorSemanticsKind::SphNativeParticleKernel;
    }
};

struct MatrixFreeAPhiSources
{
    StdVec<Complex> source_ax;
    StdVec<Complex> source_ay;
    StdVec<Complex> source_az;
    StdVec<Complex> source_phi;
};

struct MatrixFreeAPhiResiduals
{
    StdVec<Complex> residual_ax;
    StdVec<Complex> residual_ay;
    StdVec<Complex> residual_az;
    StdVec<Complex> residual_phi;
    StdVec<Complex> divergence_a;
    StdVec<Complex> divergence_j;

    Real residual_ax_l2;
    Real residual_ay_l2;
    Real residual_az_l2;
    Real residual_phi_l2;
    Real divergence_a_l2;
    Real divergence_j_l2;

    size_t interface_edge_count;
    size_t high_contrast_edge_count;
    size_t high_contrast_particle_count;
    Real max_interface_sigma_ratio;
    Real high_contrast_sigma_grad_phi_l2;
    Real high_contrast_div_sigma_a_l2;
    Real high_contrast_residual_a_l2;
    Real high_contrast_residual_phi_l2;

    MatrixFreeAPhiResiduals();
    void resize(size_t size);
    void clear();
};

MatrixFreeAPhiResiduals evaluateMatrixFreeAPhiResiduals(const MatrixFreePairwiseGraph &graph,
                                                        const MatrixFreeAPhiFields &fields,
                                                        const StdVec<Real> &electrical_conductivity,
                                                        const StdVec<Real> &magnetic_reluctivity,
                                                        const MatrixFreeAPhiParameters &parameters,
                                                        const MatrixFreeAPhiSources &sources,
                                                        bool enable_gauge_penalty,
                                                        Real gauge_penalty_coefficient,
                                                        MatrixFreeAPhiSphNativeContext *sph_native_context = nullptr);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_residuals.hpp"

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_APHI_RESIDUALS_H
