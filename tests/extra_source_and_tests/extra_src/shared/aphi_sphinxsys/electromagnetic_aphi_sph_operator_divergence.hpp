#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_HPP
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_divergence.h"

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
constexpr const char *kDivScratchXRealName = "APhiSphDivScratchXReal";
constexpr const char *kDivScratchXImagName = "APhiSphDivScratchXImag";
constexpr const char *kDivScratchYRealName = "APhiSphDivScratchYReal";
constexpr const char *kDivScratchYImagName = "APhiSphDivScratchYImag";
constexpr const char *kDivScratchZRealName = "APhiSphDivScratchZReal";
constexpr const char *kDivScratchZImagName = "APhiSphDivScratchZImag";
constexpr const char *kDivHarmonicEdgeWeightName = "APhiSphDivHarmonicEdgeWeight";
constexpr const char *kDivOutRealName = "APhiSphDivOutReal";
constexpr const char *kDivOutImagName = "APhiSphDivOutImag";
} // namespace

template <typename... RelationTypes>
class APhiStandardVectorDivergenceCK;

template <template <typename...> class RelationType, typename... Parameters>
class APhiStandardVectorDivergenceCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit APhiStandardVectorDivergenceCK(BaseRelationType &base_relation)
        : Interaction<RelationType<Parameters...>>(base_relation),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_scratch_x_real_(this->particles_->template getVariableByName<Real>(kDivScratchXRealName)),
          dv_scratch_x_imag_(this->particles_->template getVariableByName<Real>(kDivScratchXImagName)),
          dv_scratch_y_real_(this->particles_->template getVariableByName<Real>(kDivScratchYRealName)),
          dv_scratch_y_imag_(this->particles_->template getVariableByName<Real>(kDivScratchYImagName)),
          dv_scratch_z_real_(this->particles_->template getVariableByName<Real>(kDivScratchZRealName)),
          dv_scratch_z_imag_(this->particles_->template getVariableByName<Real>(kDivScratchZImagName)),
          dv_divergence_real_(this->particles_->template getVariableByName<Real>(kDivOutRealName)),
          dv_divergence_imag_(this->particles_->template getVariableByName<Real>(kDivOutImagName))
    {
    }

    virtual ~APhiStandardVectorDivergenceCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, APhiStandardVectorDivergenceCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args)
            : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              scratch_x_real_(encloser.dv_scratch_x_real_->DelegatedData(ex_policy)),
              scratch_x_imag_(encloser.dv_scratch_x_imag_->DelegatedData(ex_policy)),
              scratch_y_real_(encloser.dv_scratch_y_real_->DelegatedData(ex_policy)),
              scratch_y_imag_(encloser.dv_scratch_y_imag_->DelegatedData(ex_policy)),
              scratch_z_real_(encloser.dv_scratch_z_real_->DelegatedData(ex_policy)),
              scratch_z_imag_(encloser.dv_scratch_z_imag_->DelegatedData(ex_policy)),
              divergence_real_(encloser.dv_divergence_real_->DelegatedData(ex_policy)),
              divergence_imag_(encloser.dv_divergence_imag_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *scratch_x_real_;
        Real *scratch_x_imag_;
        Real *scratch_y_real_;
        Real *scratch_y_imag_;
        Real *scratch_z_real_;
        Real *scratch_z_imag_;
        Real *divergence_real_;
        Real *divergence_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_scratch_x_real_;
    DiscreteVariable<Real> *dv_scratch_x_imag_;
    DiscreteVariable<Real> *dv_scratch_y_real_;
    DiscreteVariable<Real> *dv_scratch_y_imag_;
    DiscreteVariable<Real> *dv_scratch_z_real_;
    DiscreteVariable<Real> *dv_scratch_z_imag_;
    DiscreteVariable<Real> *dv_divergence_real_;
    DiscreteVariable<Real> *dv_divergence_imag_;
};

template <typename... Parameters>
class APhiStandardVectorDivergenceCK<Inner<Parameters...>> : public APhiStandardVectorDivergenceCK<Base, Inner<Parameters...>>
{
  public:
    explicit APhiStandardVectorDivergenceCK(Inner<Parameters...> &inner_relation)
        : APhiStandardVectorDivergenceCK<Base, Inner<Parameters...>>(inner_relation) {}
    virtual ~APhiStandardVectorDivergenceCK() {}
};

template <template <typename...> class RelationType, typename... Parameters>
void APhiStandardVectorDivergenceCK<Base, RelationType<Parameters...>>::InteractKernel::interact(size_t index_i, Real)
{
    const Real x_real_i = scratch_x_real_[index_i];
    const Real x_imag_i = scratch_x_imag_[index_i];
    const Real y_real_i = scratch_y_real_[index_i];
    const Real y_imag_i = scratch_y_imag_[index_i];
    const Real z_real_i = scratch_z_real_[index_i];
    const Real z_imag_i = scratch_z_imag_[index_i];
    Real div_real = 0.0;
    Real div_imag = 0.0;

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd gradient_weight = this->dW_ij(index_i, index_j) * this->Vol_[index_j] * this->e_ij(index_i, index_j);
        div_real += gradient_weight[0] * (scratch_x_real_[index_j] - x_real_i) +
                    gradient_weight[1] * (scratch_y_real_[index_j] - y_real_i) +
                    gradient_weight[2] * (scratch_z_real_[index_j] - z_real_i);
        div_imag += gradient_weight[0] * (scratch_x_imag_[index_j] - x_imag_i) +
                    gradient_weight[1] * (scratch_y_imag_[index_j] - y_imag_i) +
                    gradient_weight[2] * (scratch_z_imag_[index_j] - z_imag_i);
    }

    divergence_real_[index_i] = div_real;
    divergence_imag_[index_i] = div_imag;
}

template <typename... RelationTypes>
class APhiHarmonicVectorDivergenceCK;

template <template <typename...> class RelationType, typename... Parameters>
class APhiHarmonicVectorDivergenceCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit APhiHarmonicVectorDivergenceCK(BaseRelationType &base_relation)
        : Interaction<RelationType<Parameters...>>(base_relation),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_scratch_x_real_(this->particles_->template getVariableByName<Real>(kDivScratchXRealName)),
          dv_scratch_x_imag_(this->particles_->template getVariableByName<Real>(kDivScratchXImagName)),
          dv_scratch_y_real_(this->particles_->template getVariableByName<Real>(kDivScratchYRealName)),
          dv_scratch_y_imag_(this->particles_->template getVariableByName<Real>(kDivScratchYImagName)),
          dv_scratch_z_real_(this->particles_->template getVariableByName<Real>(kDivScratchZRealName)),
          dv_scratch_z_imag_(this->particles_->template getVariableByName<Real>(kDivScratchZImagName)),
          dv_edge_weight_(this->particles_->template getVariableByName<Real>(kDivHarmonicEdgeWeightName)),
          dv_divergence_real_(this->particles_->template getVariableByName<Real>(kDivOutRealName)),
          dv_divergence_imag_(this->particles_->template getVariableByName<Real>(kDivOutImagName))
    {
    }

    virtual ~APhiHarmonicVectorDivergenceCK() {}

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, APhiHarmonicVectorDivergenceCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args)
            : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              scratch_x_real_(encloser.dv_scratch_x_real_->DelegatedData(ex_policy)),
              scratch_x_imag_(encloser.dv_scratch_x_imag_->DelegatedData(ex_policy)),
              scratch_y_real_(encloser.dv_scratch_y_real_->DelegatedData(ex_policy)),
              scratch_y_imag_(encloser.dv_scratch_y_imag_->DelegatedData(ex_policy)),
              scratch_z_real_(encloser.dv_scratch_z_real_->DelegatedData(ex_policy)),
              scratch_z_imag_(encloser.dv_scratch_z_imag_->DelegatedData(ex_policy)),
              edge_weight_(encloser.dv_edge_weight_->DelegatedData(ex_policy)),
              divergence_real_(encloser.dv_divergence_real_->DelegatedData(ex_policy)),
              divergence_imag_(encloser.dv_divergence_imag_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *scratch_x_real_;
        Real *scratch_x_imag_;
        Real *scratch_y_real_;
        Real *scratch_y_imag_;
        Real *scratch_z_real_;
        Real *scratch_z_imag_;
        Real *edge_weight_;
        Real *divergence_real_;
        Real *divergence_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_scratch_x_real_;
    DiscreteVariable<Real> *dv_scratch_x_imag_;
    DiscreteVariable<Real> *dv_scratch_y_real_;
    DiscreteVariable<Real> *dv_scratch_y_imag_;
    DiscreteVariable<Real> *dv_scratch_z_real_;
    DiscreteVariable<Real> *dv_scratch_z_imag_;
    DiscreteVariable<Real> *dv_edge_weight_;
    DiscreteVariable<Real> *dv_divergence_real_;
    DiscreteVariable<Real> *dv_divergence_imag_;
};

template <typename... Parameters>
class APhiHarmonicVectorDivergenceCK<Inner<Parameters...>> : public APhiHarmonicVectorDivergenceCK<Base, Inner<Parameters...>>
{
  public:
    explicit APhiHarmonicVectorDivergenceCK(Inner<Parameters...> &inner_relation)
        : APhiHarmonicVectorDivergenceCK<Base, Inner<Parameters...>>(inner_relation) {}
    virtual ~APhiHarmonicVectorDivergenceCK() {}
};

template <template <typename...> class RelationType, typename... Parameters>
void APhiHarmonicVectorDivergenceCK<Base, RelationType<Parameters...>>::InteractKernel::interact(size_t index_i, Real)
{
    const Real x_real_i = scratch_x_real_[index_i];
    const Real x_imag_i = scratch_x_imag_[index_i];
    const Real y_real_i = scratch_y_real_[index_i];
    const Real y_imag_i = scratch_y_imag_[index_i];
    const Real z_real_i = scratch_z_real_[index_i];
    const Real z_imag_i = scratch_z_imag_[index_i];
    Real div_real = 0.0;
    Real div_imag = 0.0;

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real coefficient_ij = harmonicMeanCoefficient(edge_weight_[index_i], edge_weight_[index_j]);
        const Vecd gradient_weight = this->dW_ij(index_i, index_j) * this->Vol_[index_j] * this->e_ij(index_i, index_j);
        div_real += coefficient_ij * (gradient_weight[0] * (scratch_x_real_[index_j] - x_real_i) +
                                      gradient_weight[1] * (scratch_y_real_[index_j] - y_real_i) +
                                      gradient_weight[2] * (scratch_z_real_[index_j] - z_real_i));
        div_imag += coefficient_ij * (gradient_weight[0] * (scratch_x_imag_[index_j] - x_imag_i) +
                                      gradient_weight[1] * (scratch_y_imag_[index_j] - y_imag_i) +
                                      gradient_weight[2] * (scratch_z_imag_[index_j] - z_imag_i));
    }

    divergence_real_[index_i] = div_real;
    divergence_imag_[index_i] = div_imag;
}

struct SPHComplexStandardVectorDivergenceOperatorDynamics
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    explicit SPHComplexStandardVectorDivergenceOperatorDynamics(RealBody &body, Inner<> &inner_ck_relation)
        : body_(body), inner_ck_(inner_ck_relation), update_cell_linked_list_(body_), update_inner_relation_(inner_ck_),
          vector_divergence_(inner_ck_)
    {
    }

    void execOperators()
    {
        update_cell_linked_list_.exec();
        update_inner_relation_.exec();
        vector_divergence_.exec();
    }

    RealBody &body_;
    Inner<> &inner_ck_;
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_;
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_;
    InteractionDynamicsCK<MainExecutionPolicy, APhiStandardVectorDivergenceCK<Inner<>>> vector_divergence_;
};

struct SPHComplexHarmonicVectorDivergenceOperatorDynamics
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    explicit SPHComplexHarmonicVectorDivergenceOperatorDynamics(RealBody &body, Inner<> &inner_ck_relation)
        : body_(body), inner_ck_(inner_ck_relation), update_cell_linked_list_(body_), update_inner_relation_(inner_ck_),
          vector_divergence_(inner_ck_)
    {
    }

    void execOperators()
    {
        update_cell_linked_list_.exec();
        update_inner_relation_.exec();
        vector_divergence_.exec();
    }

    RealBody &body_;
    Inner<> &inner_ck_;
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_;
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_;
    InteractionDynamicsCK<MainExecutionPolicy, APhiHarmonicVectorDivergenceCK<Inner<>>> vector_divergence_;
};

inline SPHComplexStandardVectorDivergenceOperator::SPHComplexStandardVectorDivergenceOperator(RealBody &body,
                                                                                              Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexStandardVectorDivergenceOperator::~SPHComplexStandardVectorDivergenceOperator() = default;

inline void SPHComplexStandardVectorDivergenceOperator::ensureInitialized()
{
    if (initialized_)
    {
        return;
    }

    BaseParticles &particles = body_.getBaseParticles();
    particles.registerStateVariable<Real>(kDivScratchXRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchXImagName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchYRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchYImagName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchZRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchZImagName, Real(0));
    particles.registerStateVariable<Real>(kDivOutRealName, Real(0));
    particles.registerStateVariable<Real>(kDivOutImagName, Real(0));
    dynamics_ = std::make_unique<SPHComplexStandardVectorDivergenceOperatorDynamics>(body_, ck_inner_relation_);
    initialized_ = true;
}

inline void SPHComplexStandardVectorDivergenceOperator::copyComponentsToParticles(const StdVec<Complex> &field_x,
                                                                                const StdVec<Complex> &field_y,
                                                                                const StdVec<Complex> &field_z)
{
    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (field_x.size() != number_of_particles || field_y.size() != number_of_particles || field_z.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexStandardVectorDivergenceOperator size mismatch");
    }

    Real *scratch_x_real = particles.getVariableDataByName<Real>(kDivScratchXRealName);
    Real *scratch_x_imag = particles.getVariableDataByName<Real>(kDivScratchXImagName);
    Real *scratch_y_real = particles.getVariableDataByName<Real>(kDivScratchYRealName);
    Real *scratch_y_imag = particles.getVariableDataByName<Real>(kDivScratchYImagName);
    Real *scratch_z_real = particles.getVariableDataByName<Real>(kDivScratchZRealName);
    Real *scratch_z_imag = particles.getVariableDataByName<Real>(kDivScratchZImagName);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        scratch_x_real[i] = field_x[i].real();
        scratch_x_imag[i] = field_x[i].imag();
        scratch_y_real[i] = field_y[i].real();
        scratch_y_imag[i] = field_y[i].imag();
        scratch_z_real[i] = field_z[i].real();
        scratch_z_imag[i] = field_z[i].imag();
    }
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kDivScratchXRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchXImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchYRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchYImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchZRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchZImagName)->finalizeLoadIn(execution::par_device);
#endif
    body_.setNewlyUpdated();
}

inline StdVec<Complex> SPHComplexStandardVectorDivergenceOperator::computeFromComponents(const StdVec<Complex> &field_x,
                                                                                         const StdVec<Complex> &field_y,
                                                                                         const StdVec<Complex> &field_z)
{
    ensureInitialized();
    copyComponentsToParticles(field_x, field_y, field_z);
    dynamics_->execOperators();

    BaseParticles &particles = body_.getBaseParticles();
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kDivOutRealName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Real>(kDivOutImagName)->prepareForOutput(execution::par_device);
#endif
    const size_t number_of_particles = particles.TotalRealParticles();
    const Real *divergence_real = particles.getVariableDataByName<Real>(kDivOutRealName);
    const Real *divergence_imag = particles.getVariableDataByName<Real>(kDivOutImagName);

    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = Complex(divergence_real[i], divergence_imag[i]);
    }
    return divergence;
}

inline SPHComplexHarmonicVectorDivergenceOperator::SPHComplexHarmonicVectorDivergenceOperator(RealBody &body,
                                                                                              Inner<> &ck_inner_relation)
    : body_(body), ck_inner_relation_(ck_inner_relation)
{
}

inline SPHComplexHarmonicVectorDivergenceOperator::~SPHComplexHarmonicVectorDivergenceOperator() = default;

inline void SPHComplexHarmonicVectorDivergenceOperator::ensureInitialized()
{
    if (initialized_)
    {
        return;
    }

    BaseParticles &particles = body_.getBaseParticles();
    particles.registerStateVariable<Real>(kDivScratchXRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchXImagName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchYRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchYImagName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchZRealName, Real(0));
    particles.registerStateVariable<Real>(kDivScratchZImagName, Real(0));
    particles.registerStateVariable<Real>(kDivHarmonicEdgeWeightName, Real(1));
    particles.registerStateVariable<Real>(kDivOutRealName, Real(0));
    particles.registerStateVariable<Real>(kDivOutImagName, Real(0));
    dynamics_ = std::make_unique<SPHComplexHarmonicVectorDivergenceOperatorDynamics>(body_, ck_inner_relation_);
    initialized_ = true;
}

inline void SPHComplexHarmonicVectorDivergenceOperator::copyComponentsAndEdgeWeightsToParticles(
    const StdVec<Complex> &field_x, const StdVec<Complex> &field_y, const StdVec<Complex> &field_z,
    const StdVec<Real> &edge_weight_coefficient)
{
    BaseParticles &particles = body_.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    if (field_x.size() != number_of_particles || field_y.size() != number_of_particles || field_z.size() != number_of_particles ||
        edge_weight_coefficient.size() != number_of_particles)
    {
        throw std::runtime_error("SPHComplexHarmonicVectorDivergenceOperator size mismatch");
    }

    Real *scratch_x_real = particles.getVariableDataByName<Real>(kDivScratchXRealName);
    Real *scratch_x_imag = particles.getVariableDataByName<Real>(kDivScratchXImagName);
    Real *scratch_y_real = particles.getVariableDataByName<Real>(kDivScratchYRealName);
    Real *scratch_y_imag = particles.getVariableDataByName<Real>(kDivScratchYImagName);
    Real *scratch_z_real = particles.getVariableDataByName<Real>(kDivScratchZRealName);
    Real *scratch_z_imag = particles.getVariableDataByName<Real>(kDivScratchZImagName);
    Real *edge_weight = particles.getVariableDataByName<Real>(kDivHarmonicEdgeWeightName);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        scratch_x_real[i] = field_x[i].real();
        scratch_x_imag[i] = field_x[i].imag();
        scratch_y_real[i] = field_y[i].real();
        scratch_y_imag[i] = field_y[i].imag();
        scratch_z_real[i] = field_z[i].real();
        scratch_z_imag[i] = field_z[i].imag();
        edge_weight[i] = edge_weight_coefficient[i];
    }
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kDivScratchXRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchXImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchYRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchYImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchZRealName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivScratchZImagName)->finalizeLoadIn(execution::par_device);
    particles.getVariableByName<Real>(kDivHarmonicEdgeWeightName)->finalizeLoadIn(execution::par_device);
#endif
    body_.setNewlyUpdated();
}

inline StdVec<Complex> SPHComplexHarmonicVectorDivergenceOperator::computeFromComponents(
    const StdVec<Complex> &field_x, const StdVec<Complex> &field_y, const StdVec<Complex> &field_z,
    const StdVec<Real> &edge_weight_coefficient)
{
    ensureInitialized();
    copyComponentsAndEdgeWeightsToParticles(field_x, field_y, field_z, edge_weight_coefficient);
    dynamics_->execOperators();

    BaseParticles &particles = body_.getBaseParticles();
#if SPHINXSYS_USE_SYCL
    particles.getVariableByName<Real>(kDivOutRealName)->prepareForOutput(execution::par_device);
    particles.getVariableByName<Real>(kDivOutImagName)->prepareForOutput(execution::par_device);
#endif
    const size_t number_of_particles = particles.TotalRealParticles();
    const Real *divergence_real = particles.getVariableDataByName<Real>(kDivOutRealName);
    const Real *divergence_imag = particles.getVariableDataByName<Real>(kDivOutImagName);

    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = Complex(divergence_real[i], divergence_imag[i]);
    }
    return divergence;
}

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_HPP
