#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_HPP
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_weighted_gradient.h"

#include "aphi_sphinxsys/electromagnetic_aphi_sph_coefficient_policy.h"
#include "interaction_ck.hpp"

#include <stdexcept>

namespace SPH
{
namespace electromagnetics
{
namespace sph
{
namespace
{
constexpr const char *kHarmonicScratchRealName = "APhiSphHarmonicScratchReal";
constexpr const char *kHarmonicScratchImagName = "APhiSphHarmonicScratchImag";
constexpr const char *kHarmonicEdgeWeightName = "APhiSphHarmonicEdgeWeight";
constexpr const char *kHarmonicGradRealName = "APhiSphHarmonicGradReal";
constexpr const char *kHarmonicGradImagName = "APhiSphHarmonicGradImag";
} // namespace

template <typename... RelationTypes>
class APhiHarmonicWeightedGradientCK;

template <template <typename...> class RelationType, typename... Parameters>
class APhiHarmonicWeightedGradientCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit APhiHarmonicWeightedGradientCK(BaseRelationType &base_relation)
        : Interaction<RelationType<Parameters...>>(base_relation),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_scratch_real_(this->particles_->template getVariableByName<Real>(kHarmonicScratchRealName)),
          dv_scratch_imag_(this->particles_->template getVariableByName<Real>(kHarmonicScratchImagName)),
          dv_edge_weight_(this->particles_->template getVariableByName<Real>(kHarmonicEdgeWeightName)),
          dv_gradient_real_(this->particles_->template getVariableByName<Vecd>(kHarmonicGradRealName)),
          dv_gradient_imag_(this->particles_->template getVariableByName<Vecd>(kHarmonicGradImagName))
    {
    }

    virtual ~APhiHarmonicWeightedGradientCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, APhiHarmonicWeightedGradientCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args)
            : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              scratch_real_(encloser.dv_scratch_real_->DelegatedData(ex_policy)),
              scratch_imag_(encloser.dv_scratch_imag_->DelegatedData(ex_policy)),
              edge_weight_(encloser.dv_edge_weight_->DelegatedData(ex_policy)),
              gradient_real_(encloser.dv_gradient_real_->DelegatedData(ex_policy)),
              gradient_imag_(encloser.dv_gradient_imag_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *scratch_real_;
        Real *scratch_imag_;
        Real *edge_weight_;
        Vecd *gradient_real_;
        Vecd *gradient_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_scratch_real_;
    DiscreteVariable<Real> *dv_scratch_imag_;
    DiscreteVariable<Real> *dv_edge_weight_;
    DiscreteVariable<Vecd> *dv_gradient_real_;
    DiscreteVariable<Vecd> *dv_gradient_imag_;
};

template <typename... Parameters>
class APhiHarmonicWeightedGradientCK<Inner<Parameters...>> : public APhiHarmonicWeightedGradientCK<Base, Inner<Parameters...>>
{
  public:
    explicit APhiHarmonicWeightedGradientCK(Inner<Parameters...> &inner_relation)
        : APhiHarmonicWeightedGradientCK<Base, Inner<Parameters...>>(inner_relation) {}
    virtual ~APhiHarmonicWeightedGradientCK() {}
};

template <template <typename...> class RelationType, typename... Parameters>
void APhiHarmonicWeightedGradientCK<Base, RelationType<Parameters...>>::InteractKernel::interact(size_t index_i, Real)
{
    const Real value_real_i = scratch_real_[index_i];
    const Real value_imag_i = scratch_imag_[index_i];
    Vecd gradient_real_i = Vecd::Zero();
    Vecd gradient_imag_i = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real coefficient_ij = harmonicMeanCoefficient(edge_weight_[index_i], edge_weight_[index_j]);
        const Real delta_real = coefficient_ij * (scratch_real_[index_j] - value_real_i);
        const Real delta_imag = coefficient_ij * (scratch_imag_[index_j] - value_imag_i);
        const Vecd gradient_weight = this->dW_ij(index_i, index_j) * this->Vol_[index_j] * this->e_ij(index_i, index_j);
        gradient_real_i += gradient_weight * delta_real;
        gradient_imag_i += gradient_weight * delta_imag;
    }

    gradient_real_[index_i] = gradient_real_i;
    gradient_imag_[index_i] = gradient_imag_i;
}

struct SPHComplexHarmonicWeightedGradientOperatorDynamics
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    explicit SPHComplexHarmonicWeightedGradientOperatorDynamics(RealBody &body, Inner<> &inner_ck_relation)
        : body_(body), inner_ck_(inner_ck_relation), update_cell_linked_list_(body_), update_inner_relation_(inner_ck_),
          harmonic_weighted_gradient_(inner_ck_)
    {
    }

    void execOperators()
    {
        update_cell_linked_list_.exec();
        update_inner_relation_.exec();
        harmonic_weighted_gradient_.exec();
    }

    RealBody &body_;
    Inner<> &inner_ck_;
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_;
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_;
    InteractionDynamicsCK<MainExecutionPolicy, APhiHarmonicWeightedGradientCK<Inner<>>> harmonic_weighted_gradient_;
};

inline SPHComplexHarmonicWeightedGradientOperator::SPHComplexHarmonicWeightedGradientOperator(RealBody &body,
                                                                                              Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexHarmonicWeightedGradientOperator::~SPHComplexHarmonicWeightedGradientOperator() = default;

inline void SPHComplexHarmonicWeightedGradientOperator::ensureInitialized()
{
    if (initialized_)
    {
        return;
    }

    BaseParticles &particles = body_.getBaseParticles();
    particles.registerStateVariable<Real>(kHarmonicScratchRealName, Real(0));
    particles.registerStateVariable<Real>(kHarmonicScratchImagName, Real(0));
    particles.registerStateVariable<Real>(kHarmonicEdgeWeightName, Real(1));
    particles.registerStateVariable<Vecd>(kHarmonicGradRealName, Vecd::Zero().eval());
    particles.registerStateVariable<Vecd>(kHarmonicGradImagName, Vecd::Zero().eval());
    dynamics_ = std::make_unique<SPHComplexHarmonicWeightedGradientOperatorDynamics>(body_, ck_inner_relation_);
    initialized_ = true;
}

inline void SPHComplexHarmonicWeightedGradientOperator::copyFieldAndEdgeWeightsToParticles(
    const StdVec<Complex> &field, const StdVec<Real> &edge_weight_coefficient)
{
    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (field.size() != number_of_particles || edge_weight_coefficient.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexHarmonicWeightedGradientOperator size mismatch");
    }

    Real *scratch_real = particles.getVariableDataByName<Real>(kHarmonicScratchRealName);
    Real *scratch_imag = particles.getVariableDataByName<Real>(kHarmonicScratchImagName);
    Real *edge_weight = particles.getVariableDataByName<Real>(kHarmonicEdgeWeightName);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        scratch_real[i] = field[i].real();
        scratch_imag[i] = field[i].imag();
        edge_weight[i] = edge_weight_coefficient[i];
    }
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kHarmonicScratchRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kHarmonicScratchImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kHarmonicEdgeWeightName)->finalizeLoadIn(execution::par_device);
#endif
    body_.setNewlyUpdated();
}

inline StdVec<Vec3c> SPHComplexHarmonicWeightedGradientOperator::computeFromField(const StdVec<Complex> &field,
                                                                                  const StdVec<Real> &edge_weight_coefficient)
{
    ensureInitialized();
    copyFieldAndEdgeWeightsToParticles(field, edge_weight_coefficient);
    dynamics_->execOperators();

    BaseParticles &particles = body_.getBaseParticles();
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Vecd>(kHarmonicGradRealName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Vecd>(kHarmonicGradImagName)->prepareForOutput(execution::par_device);
#endif
    const size_t number_of_particles = particles.TotalRealParticles();
    const Vecd *gradient_real = particles.getVariableDataByName<Vecd>(kHarmonicGradRealName);
    const Vecd *gradient_imag = particles.getVariableDataByName<Vecd>(kHarmonicGradImagName);

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

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_HPP
