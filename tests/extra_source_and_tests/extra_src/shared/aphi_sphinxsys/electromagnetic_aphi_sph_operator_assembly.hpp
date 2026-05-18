#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_HPP
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_assembly.h"

#include <chrono>

namespace SPH
{
namespace electromagnetics
{
namespace sph
{

inline SPHAPhiOperatorAssembly::SPHAPhiOperatorAssembly(RealBody &body, Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexScalarGradientOperator &SPHAPhiOperatorAssembly::gradientOperator()
{
    if (!gradient_)
    {
        gradient_ = std::make_unique<SPHComplexScalarGradientOperator>(body_, ck_inner_relation_);
    }
    return *gradient_;
}

inline SPHComplexScalarNegativeLaplaceOperator &SPHAPhiOperatorAssembly::laplaceOperator()
{
    if (!laplace_)
    {
        laplace_ = std::make_unique<SPHComplexScalarNegativeLaplaceOperator>(body_, ck_inner_relation_,
                                                                             SPHScalarNegativeLaplaceParameters{});
    }
    return *laplace_;
}

inline const matrix_free::MatrixFreeAPhiSphNativeDiagnostics &SPHAPhiOperatorAssembly::diagnostics() const
{
    return diagnostics_;
}

inline StdVec<Vec3c> SPHAPhiOperatorAssembly::computeStandardGradient(const StdVec<Complex> &field)
{
    const auto start = std::chrono::steady_clock::now();
    StdVec<Vec3c> gradient = gradientOperator().computeFromField(field);
    diagnostics_.standard_gradient_calls_ += 1;
    diagnostics_.standard_gradient_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return gradient;
}

inline void SPHAPhiOperatorAssembly::accumulateScalarLaplaceHelmholtzResiduals(
    const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient,
    matrix_free::ScalarComplexHelmholtzResiduals &residuals)
{
    const auto start = std::chrono::steady_clock::now();
    laplaceOperator().accumulateHelmholtzLaplaceResiduals(field, diffusion_coefficient, residuals);
    diagnostics_.laplace_helmholtz_calls_ += 1;
    diagnostics_.laplace_helmholtz_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

inline void SPHAPhiOperatorAssembly::accumulateScalarDivergenceOfGradientHelmholtzResiduals(
    const StdVec<Complex> &field, matrix_free::ScalarComplexHelmholtzResiduals &residuals)
{
    const auto start = std::chrono::steady_clock::now();
    const size_t number_of_particles = field.size();
    const StdVec<Vec3c> gradient = computeStandardGradient(field);

    StdVec<Complex> gradient_x(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> gradient_y(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> gradient_z(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        gradient_x[i] = gradient[i][0];
        gradient_y[i] = gradient[i][1];
        gradient_z[i] = gradient[i][2];
    }

    const StdVec<Vec3c> grad_gradient_x = computeStandardGradient(gradient_x);
    const StdVec<Vec3c> grad_gradient_y = computeStandardGradient(gradient_y);
    const StdVec<Vec3c> grad_gradient_z = computeStandardGradient(gradient_z);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.laplace_term_[i] += grad_gradient_x[i][0] + grad_gradient_y[i][1] + grad_gradient_z[i][2];
    }

    matrix_free::ScalarComplexHelmholtzResiduals diagonal_proxy;
    diagonal_proxy.resize(number_of_particles);
    diagonal_proxy.clear();
    StdVec<Real> unit_diffusion(number_of_particles, Real(1));
    accumulateScalarLaplaceHelmholtzResiduals(field, unit_diffusion, diagonal_proxy);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.diagonal_scale_[i] += diagonal_proxy.diagonal_scale_[i];
    }
    diagnostics_.divergence_of_gradient_calls_ += 1;
    diagnostics_.divergence_of_gradient_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

inline StdVec<Vec3c> SPHAPhiOperatorAssembly::computeHarmonicWeightedGradient(const StdVec<Complex> &field,
                                                                              const StdVec<Real> &edge_weight_coefficient)
{
    const auto start = std::chrono::steady_clock::now();
    if (!harmonic_gradient_)
    {
        harmonic_gradient_ = std::make_unique<SPHComplexHarmonicWeightedGradientOperator>(body_, ck_inner_relation_);
    }
    StdVec<Vec3c> gradient = harmonic_gradient_->computeFromField(field, edge_weight_coefficient);
    diagnostics_.harmonic_weighted_gradient_calls_ += 1;
    diagnostics_.harmonic_weighted_gradient_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return gradient;
}

inline StdVec<Vec3c> SPHAPhiOperatorAssembly::computeSigmaGradPhi(const StdVec<Complex> &phi, const StdVec<Real> &sigma)
{
    const auto start = std::chrono::steady_clock::now();
    StdVec<Vec3c> sigma_grad_phi = computeHarmonicWeightedGradient(phi, sigma);
    diagnostics_.sigma_grad_phi_calls_ += 1;
    diagnostics_.sigma_grad_phi_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return sigma_grad_phi;
}

inline StdVec<Complex> SPHAPhiOperatorAssembly::computeGradientDivergenceOfVectorField(const StdVec<Complex> &field_x,
                                                                                       const StdVec<Complex> &field_y,
                                                                                       const StdVec<Complex> &field_z)
{
    const auto start = std::chrono::steady_clock::now();
    if (!standard_divergence_)
    {
        standard_divergence_ = std::make_unique<SPHComplexStandardVectorDivergenceOperator>(body_, ck_inner_relation_);
    }
    StdVec<Complex> divergence = standard_divergence_->computeFromComponents(field_x, field_y, field_z);
    diagnostics_.standard_divergence_calls_ += 1;
    diagnostics_.standard_divergence_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return divergence;
}

inline StdVec<Complex> SPHAPhiOperatorAssembly::computeHarmonicDivergenceOfVectorField(const StdVec<Complex> &field_x,
                                                                                       const StdVec<Complex> &field_y,
                                                                                       const StdVec<Complex> &field_z,
                                                                                       const StdVec<Real> &sigma)
{
    const auto start = std::chrono::steady_clock::now();
    const StdVec<Vec3c> grad_x = computeHarmonicWeightedGradient(field_x, sigma);
    const StdVec<Vec3c> grad_y = computeHarmonicWeightedGradient(field_y, sigma);
    const StdVec<Vec3c> grad_z = computeHarmonicWeightedGradient(field_z, sigma);

    StdVec<Complex> divergence(field_x.size(), Complex(0.0, 0.0));
    for (size_t i = 0; i != field_x.size(); ++i)
    {
        divergence[i] = grad_x[i][0] + grad_y[i][1] + grad_z[i][2];
    }
    diagnostics_.harmonic_divergence_calls_ += 1;
    diagnostics_.harmonic_divergence_seconds_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return divergence;
}

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_ASSEMBLY_HPP
