#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_HPP
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_gradient.h"

#include "general_gradient.h"
#include "kernel_correction_ck.h"

#include <stdexcept>

namespace SPH
{
namespace electromagnetics
{
namespace sph
{
namespace
{
constexpr const char *kScratchRealVariableName = "APhiSphGradScratchReal";
constexpr const char *kScratchImagVariableName = "APhiSphGradScratchImag";
constexpr const char *kScratchRealGradientName = "APhiSphGradScratchRealGradient";
constexpr const char *kScratchImagGradientName = "APhiSphGradScratchImagGradient";
} // namespace

struct SPHComplexScalarGradientOperatorDynamics
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    explicit SPHComplexScalarGradientOperatorDynamics(RealBody &body, Inner<> &inner_ck_relation)
        : body_(body), inner_ck_(inner_ck_relation),
          update_cell_linked_list_(body_),
          update_inner_relation_(inner_ck_),
          linear_correction_matrix_(DynamicsArgs(inner_ck_, 0.0)),
          gradient_real_(DynamicsArgs(inner_ck_, std::string(kScratchRealVariableName))),
          gradient_imag_(DynamicsArgs(inner_ck_, std::string(kScratchImagVariableName)))
    {
    }

    void execOperators()
    {
        update_cell_linked_list_.exec();
        update_inner_relation_.exec();
        linear_correction_matrix_.exec();
        gradient_real_.exec();
        gradient_imag_.exec();
    }

    RealBody &body_;
    Inner<> &inner_ck_;
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_;
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_;
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix_;
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> gradient_real_;
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> gradient_imag_;
};

inline SPHComplexScalarGradientOperator::SPHComplexScalarGradientOperator(RealBody &body, Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexScalarGradientOperator::~SPHComplexScalarGradientOperator() = default;

inline void SPHComplexScalarGradientOperator::ensureInitialized()
{
    if (initialized_)
    {
        return;
    }

    BaseParticles &particles = body_.getBaseParticles();
    particles.registerStateVariable<Real>(kScratchRealVariableName, Real(0));
    particles.registerStateVariable<Real>(kScratchImagVariableName, Real(0));
    dynamics_ = std::make_unique<SPHComplexScalarGradientOperatorDynamics>(body_, ck_inner_relation_);
    initialized_ = true;
}

inline void SPHComplexScalarGradientOperator::copyFieldToParticleScalars(const StdVec<Complex> &field)
{
    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (field.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexScalarGradientOperator field size mismatch");
    }

    Real *scratch_real = particles.getVariableDataByName<Real>(kScratchRealVariableName);
    Real *scratch_imag = particles.getVariableDataByName<Real>(kScratchImagVariableName);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        scratch_real[i] = field[i].real();
        scratch_imag[i] = field[i].imag();
    }
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kScratchRealVariableName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kScratchImagVariableName)->finalizeLoadIn(execution::par_device);
#endif
    body_.setNewlyUpdated();
}

inline StdVec<Vec3c> SPHComplexScalarGradientOperator::computeFromField(const StdVec<Complex> &field)
{
    ensureInitialized();
    copyFieldToParticleScalars(field);
    dynamics_->execOperators();

    BaseParticles &particles = body_.getBaseParticles();
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Vecd>(kScratchRealGradientName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Vecd>(kScratchImagGradientName)->prepareForOutput(execution::par_device);
#endif
    const size_t number_of_particles = particles.TotalRealParticles();
    const Vecd *gradient_real = particles.getVariableDataByName<Vecd>(kScratchRealGradientName);
    const Vecd *gradient_imag = particles.getVariableDataByName<Vecd>(kScratchImagGradientName);

    StdVec<Vec3c> gradient(number_of_particles, Vec3c::Zero());
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            gradient[i][axis] = Complex(gradient_real[i][axis], gradient_imag[i][axis]);
        }
    }
    return gradient;
}

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_HPP
