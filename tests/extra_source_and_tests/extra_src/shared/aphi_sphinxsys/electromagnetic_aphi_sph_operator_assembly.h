#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_H
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_diagnostics.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_divergence.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_gradient.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_laplace.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_weighted_gradient.h"
#include <complex>
#include <memory>

namespace SPH
{
template <typename...>
class Inner;

namespace electromagnetics
{
namespace sph
{

/**
 * Lazily-initialized SPH-native operator assembly for matrix-free A-phi.
 * This keeps operator ownership and lightweight diagnostics close to the CK
 * execution layer so higher-level contexts remain thin orchestration shells.
 */
class SPHAPhiOperatorAssembly
{
  public:
    explicit SPHAPhiOperatorAssembly(RealBody &body, Inner<> &ck_inner_relation);

    SPHComplexScalarGradientOperator &gradientOperator();
    SPHComplexScalarNegativeLaplaceOperator &laplaceOperator();
    const matrix_free::MatrixFreeAPhiSphNativeDiagnostics &diagnostics() const;

    StdVec<Vec3c> computeStandardGradient(const StdVec<Complex> &field);
    void accumulateScalarLaplaceHelmholtzResiduals(const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient,
                                                   matrix_free::ScalarComplexHelmholtzResiduals &residuals);
    void accumulateScalarDivergenceOfGradientHelmholtzResiduals(const StdVec<Complex> &field,
                                                                matrix_free::ScalarComplexHelmholtzResiduals &residuals);
    StdVec<Vec3c> computeSigmaGradPhi(const StdVec<Complex> &phi, const StdVec<Real> &sigma);
    /** Harmonic-weighted gradient with caller-supplied edge weights; used as a reusable coupled-assembly primitive. */
    StdVec<Vec3c> computeHarmonicWeightedGradient(const StdVec<Complex> &field, const StdVec<Real> &edge_weight_coefficient);
    StdVec<Complex> computeGradientDivergenceOfVectorField(const StdVec<Complex> &field_x,
                                                           const StdVec<Complex> &field_y,
                                                           const StdVec<Complex> &field_z);
    StdVec<Complex> computeHarmonicDivergenceOfVectorField(const StdVec<Complex> &field_x,
                                                           const StdVec<Complex> &field_y,
                                                           const StdVec<Complex> &field_z,
                                                           const StdVec<Real> &sigma);

  private:
    RealBody &body_;
    Inner<> &ck_inner_relation_;
    matrix_free::MatrixFreeAPhiSphNativeDiagnostics diagnostics_;
    std::unique_ptr<SPHComplexScalarGradientOperator> gradient_;
    std::unique_ptr<SPHComplexScalarNegativeLaplaceOperator> laplace_;
    std::unique_ptr<SPHComplexScalarGradientOperator> standard_gradient_;
    std::unique_ptr<SPHComplexHarmonicWeightedGradientOperator> harmonic_gradient_;
    std::unique_ptr<SPHComplexStandardVectorDivergenceOperator> standard_divergence_;
};

using SPHAPhiCoupledOperators = SPHAPhiOperatorAssembly;

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_assembly.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_H
