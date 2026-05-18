#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_HPP
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_laplace.h"

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_coefficient_policy.h"
#include "interaction_ck.hpp"

#include <cmath>

#include <stdexcept>

namespace SPH
{
namespace electromagnetics
{
namespace sph
{
namespace
{
constexpr const char *kLaplaceScratchRealName = "APhiSphLaplaceScratchReal";
constexpr const char *kLaplaceScratchImagName = "APhiSphLaplaceScratchImag";
constexpr const char *kLaplaceOutRealName = "APhiSphLaplaceOutReal";
constexpr const char *kLaplaceOutImagName = "APhiSphLaplaceOutImag";
constexpr const char *kLaplaceDiffusionName = "APhiSphLaplaceDiffusion";
constexpr const char *kLaplaceDiagonalScaleName = "APhiSphLaplaceDiagonalScale";

} // namespace

template <typename... RelationTypes>
class APhiScalarNegativeLaplaceCK;

template <template <typename...> class RelationType, typename... Parameters>
class APhiScalarNegativeLaplaceCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit APhiScalarNegativeLaplaceCK(BaseRelationType &base_relation, Real pair_weight_regularization)
        : Interaction<RelationType<Parameters...>>(base_relation),
          pair_weight_regularization_(pair_weight_regularization),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_scratch_real_(this->particles_->template getVariableByName<Real>(kLaplaceScratchRealName)),
          dv_scratch_imag_(this->particles_->template getVariableByName<Real>(kLaplaceScratchImagName)),
          dv_laplace_real_(this->particles_->template getVariableByName<Real>(kLaplaceOutRealName)),
          dv_laplace_imag_(this->particles_->template getVariableByName<Real>(kLaplaceOutImagName)),
          dv_diffusion_(this->particles_->template getVariableByName<Real>(kLaplaceDiffusionName)),
          dv_diagonal_scale_(this->particles_->template getVariableByName<Real>(kLaplaceDiagonalScaleName))
    {
    }

    virtual ~APhiScalarNegativeLaplaceCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, APhiScalarNegativeLaplaceCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args)
            : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              scratch_real_(encloser.dv_scratch_real_->DelegatedData(ex_policy)),
              scratch_imag_(encloser.dv_scratch_imag_->DelegatedData(ex_policy)),
              laplace_real_(encloser.dv_laplace_real_->DelegatedData(ex_policy)),
              laplace_imag_(encloser.dv_laplace_imag_->DelegatedData(ex_policy)),
              diffusion_(encloser.dv_diffusion_->DelegatedData(ex_policy)),
              diagonal_scale_(encloser.dv_diagonal_scale_->DelegatedData(ex_policy)),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_) {}

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *scratch_real_;
        Real *scratch_imag_;
        Real *laplace_real_;
        Real *laplace_imag_;
        Real *diffusion_;
        Real *diagonal_scale_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_scratch_real_;
    DiscreteVariable<Real> *dv_scratch_imag_;
    DiscreteVariable<Real> *dv_laplace_real_;
    DiscreteVariable<Real> *dv_laplace_imag_;
    DiscreteVariable<Real> *dv_diffusion_;
    DiscreteVariable<Real> *dv_diagonal_scale_;
};

template <typename... Parameters>
class APhiScalarNegativeLaplaceCK<Inner<Parameters...>> : public APhiScalarNegativeLaplaceCK<Base, Inner<Parameters...>>
{
  public:
    explicit APhiScalarNegativeLaplaceCK(Inner<Parameters...> &inner_relation, Real pair_weight_regularization)
        : APhiScalarNegativeLaplaceCK<Base, Inner<Parameters...>>(inner_relation, pair_weight_regularization) {}
    virtual ~APhiScalarNegativeLaplaceCK() {}
};

template <template <typename...> class RelationType, typename... Parameters>
void APhiScalarNegativeLaplaceCK<Base, RelationType<Parameters...>>::InteractKernel::interact(size_t index_i, Real)
{
    const Real value_real_i = scratch_real_[index_i];
    const Real value_imag_i = scratch_imag_[index_i];
    Real laplace_real_i = 0.0;
    Real laplace_imag_i = 0.0;
    Real diagonal_scale_i = 0.0;

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j];
        const Real coefficient_ij = harmonicMeanCoefficient(diffusion_[index_i], diffusion_[index_j]);
        const Real diffusion_weight = pairwiseDiffusionWeight(dW_ijV_j, r_ij, pair_weight_regularization_, reference_smoothing_length_);
        const Real pair_weight = coefficient_ij * diffusion_weight;
        laplace_real_i += pair_weight * (value_real_i - scratch_real_[index_j]);
        laplace_imag_i += pair_weight * (value_imag_i - scratch_imag_[index_j]);
        diagonal_scale_i += std::abs(pair_weight);
    }

    laplace_real_[index_i] = laplace_real_i;
    laplace_imag_[index_i] = laplace_imag_i;
    diagonal_scale_[index_i] = diagonal_scale_i;
}

struct SPHComplexScalarNegativeLaplaceOperatorDynamics
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    explicit SPHComplexScalarNegativeLaplaceOperatorDynamics(RealBody &body, Inner<> &inner_ck_relation,
                                                             Real pair_weight_regularization)
        : body_(body), inner_ck_(inner_ck_relation), pair_weight_regularization_(pair_weight_regularization),
          update_cell_linked_list_(body_), update_inner_relation_(inner_ck_),
          negative_laplace_(inner_ck_, pair_weight_regularization_)
    {
    }

    void initializeRelation()
    {
        update_cell_linked_list_.exec();
        update_inner_relation_.exec();
    }

    void execOperators()
    {
        negative_laplace_.exec();
    }

    RealBody &body_;
    Inner<> &inner_ck_;
    Real pair_weight_regularization_;
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_;
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_;
    InteractionDynamicsCK<MainExecutionPolicy, APhiScalarNegativeLaplaceCK<Inner<>>> negative_laplace_;
};

inline SPHComplexScalarNegativeLaplaceOperator::SPHComplexScalarNegativeLaplaceOperator(
    RealBody &body, Inner<> &ck_inner_relation, SPHScalarNegativeLaplaceParameters parameters)
    : body_(body), parameters_(parameters), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexScalarNegativeLaplaceOperator::~SPHComplexScalarNegativeLaplaceOperator() = default;

inline void SPHComplexScalarNegativeLaplaceOperator::ensureInitialized()
{
    if (initialized_)
    {
        return;
    }

    BaseParticles &particles = body_.getBaseParticles();
    particles.registerStateVariable<Real>(kLaplaceScratchRealName, Real(0));
    particles.registerStateVariable<Real>(kLaplaceScratchImagName, Real(0));
    particles.registerStateVariable<Real>(kLaplaceOutRealName, Real(0));
    particles.registerStateVariable<Real>(kLaplaceOutImagName, Real(0));
    particles.registerStateVariable<Real>(kLaplaceDiagonalScaleName, Real(0));
    particles.registerStateVariable<Real>(kLaplaceDiffusionName, Real(1));
    dynamics_ = std::make_unique<SPHComplexScalarNegativeLaplaceOperatorDynamics>(
        body_, ck_inner_relation_, parameters_.pair_weight_regularization);
    initialized_ = true;
}

inline void SPHComplexScalarNegativeLaplaceOperator::ensureRelationUpdated()
{
    if (!relation_initialized_)
    {
        dynamics_->initializeRelation();
        relation_initialized_ = true;
    }
}

inline void SPHComplexScalarNegativeLaplaceOperator::copyFieldAndDiffusionToParticles(
    const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient)
{
    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (field.size() != number_of_particles || diffusion_coefficient.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexScalarNegativeLaplaceOperator size mismatch");
    }

    Real *scratch_real = particles.getVariableDataByName<Real>(kLaplaceScratchRealName);
    Real *scratch_imag = particles.getVariableDataByName<Real>(kLaplaceScratchImagName);
    Real *diffusion = particles.getVariableDataByName<Real>(kLaplaceDiffusionName);
    const bool diffusion_changed =
        last_diffusion_data_ != diffusion_coefficient.data() || last_diffusion_size_ != diffusion_coefficient.size();
    last_diffusion_changed_ = diffusion_changed;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        scratch_real[i] = field[i].real();
        scratch_imag[i] = field[i].imag();
        if (diffusion_changed)
        {
            diffusion[i] = diffusion_coefficient[i];
        }
    }
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kLaplaceScratchRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kLaplaceScratchImagName)->finalizeLoadIn(execution::par_device);
    if (diffusion_changed)
    {
        particles.getVariableByName<Real>(kLaplaceDiffusionName)->finalizeLoadIn(execution::par_device);
    }
#endif
    if (diffusion_changed)
    {
        last_diffusion_data_ = diffusion_coefficient.data();
        last_diffusion_size_ = diffusion_coefficient.size();
        diagonal_cache_valid_ = false;
    }
}

inline StdVec<Complex> SPHComplexScalarNegativeLaplaceOperator::apply(const StdVec<Complex> &field,
                                                                      const StdVec<Real> &diffusion_coefficient)
{
    ensureInitialized();
    ensureRelationUpdated();
    copyFieldAndDiffusionToParticles(field, diffusion_coefficient);
    dynamics_->execOperators();

    BaseParticles &particles = body_.getBaseParticles();
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kLaplaceOutRealName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Real>(kLaplaceOutImagName)->prepareForOutput(execution::par_device);
#endif
    const size_t number_of_particles = particles.TotalRealParticles();
    const Real *laplace_real = particles.getVariableDataByName<Real>(kLaplaceOutRealName);
    const Real *laplace_imag = particles.getVariableDataByName<Real>(kLaplaceOutImagName);

    StdVec<Complex> laplace(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        laplace[i] = Complex(laplace_real[i], laplace_imag[i]);
    }
    return laplace;
}

inline void SPHComplexScalarNegativeLaplaceOperator::accumulateHelmholtzLaplaceResiduals(
    const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient, matrix_free::ScalarComplexHelmholtzResiduals &residuals)
{
    ensureInitialized();
    ensureRelationUpdated();
    copyFieldAndDiffusionToParticles(field, diffusion_coefficient);

    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (residuals.laplace_term_.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexScalarNegativeLaplaceOperator residual size mismatch");
    }

    Real *laplace_real = particles.getVariableDataByName<Real>(kLaplaceOutRealName);
    Real *laplace_imag = particles.getVariableDataByName<Real>(kLaplaceOutImagName);
    Real *diagonal_scale = particles.getVariableDataByName<Real>(kLaplaceDiagonalScaleName);

    dynamics_->execOperators();

#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kLaplaceOutRealName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Real>(kLaplaceOutImagName)->prepareForOutput(execution::par_device);
    if (!diagonal_cache_valid_ || last_diffusion_changed_)
    {
        particles.getVariableByName<Real>(kLaplaceDiagonalScaleName)->prepareForOutput(execution::par_device);
    }
#endif
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        residuals.laplace_term_[i] = Complex(laplace_real[i], laplace_imag[i]);
        if (!diagonal_cache_valid_ || last_diffusion_changed_)
        {
            residuals.diagonal_scale_[i] = diagonal_scale[i];
        }
        else
        {
            residuals.diagonal_scale_[i] = diagonal_scale_cache_[i];
        }
    }
    if (!diagonal_cache_valid_ || last_diffusion_changed_)
    {
        diagonal_scale_cache_.assign(diagonal_scale, diagonal_scale + number_of_particles);
        diagonal_cache_valid_ = true;
    }
}

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_HPP
