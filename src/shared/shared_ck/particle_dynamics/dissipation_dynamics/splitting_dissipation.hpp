#ifndef SPLITTING_DISSIPATION_HPP
#define SPLITTING_DISSIPATION_HPP

#include "splitting_dissipation.h"

namespace SPH
{
//=================================================================================================//
template <typename DissipationType, template <typename...> class RelationType, typename... Parameters>
Dissipation<Base, DissipationType, RelationType<Parameters...>>::
    Dissipation(RelationType<Parameters...> &relation, const std::string &variable_name)
    : Interaction<RelationType<Parameters...>>(relation),
      dissipation_model_(DynamicCast<DissipationType>(this, this->particles_->getBaseMaterial())),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)) {}
//=================================================================================================//
template <typename DissipationType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
Dissipation<Base, DissipationType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : Interaction<RelationType<Parameters...>>::InteractKernel(
          ex_policy, encloser, std::forward<Args>(args)...),
      dis_coeff_(ex_policy, encloser.dissipation_model_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      zero_error_(ReduceReference<ReduceSum<DataType>>::value) {}
//=================================================================================================//
template <typename DampingType, typename... Parameters>
ConservativeDamping<Inner<Splitting, DampingType, Parameters...>>::
    ConservativeDamping(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDampingType(inner_relation, variable_name) {}
//=================================================================================================//
template <typename DampingType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ConservativeDamping<Inner<Splitting, DampingType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDampingType::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <typename DampingType, typename... Parameters>
void ConservativeDamping<Inner<Splitting, DampingType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    // compute the error and parameters
    Real parameter_a = 0.0;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    int neighbor_j = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) *
                                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                                  this->vec_r_ij(index_i, index_j).norm();

        parameter_a += parameter_b[neighbor_j];
        ++neighbor_j;
    }
    parameter_a = (parameter_a - 1.0) / this->Vol_[index_i];
    DataType error = this->zero_error_;
    Real parameter_c = 0.0;
    neighbor_j = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        DataType difference = this->variable_[index_i] - this->variable_[index_j];
        error += difference * parameter_b[neighbor_j];
        parameter_b[neighbor_j] += parameter_a * this->Vol_[index_j];
        parameter_c += parameter_b[neighbor_j] * parameter_b[neighbor_j] * getSquaredNorm(difference);
        ++neighbor_j;
    }

    // update the variable
    DataType parameter_k = error / (parameter_c + TinyReal);
    neighbor_j = 0;
    DataType ttl_out_flux = this->zero_error_;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        DataType increment = parameter_k * parameter_b[neighbor_j] *
                             getSquaredNorm(this->variable_[index_i] - this->variable_[index_j]);
        this->variable_[index_j] += increment;
        ttl_out_flux += increment * this->Vol_[index_j];
        ++neighbor_j;
    }
    this->variable_[index_i] -= ttl_out_flux / this->Vol_[index_i];
}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    ProjectionDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDissipationType::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void ProjectionDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    // compute the error and parameters
    DataType error = this->zero_error_;
    Real parameter_a = 0.0;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    Real parameter_c = 0.0;
    Real ttl_Vol = this->Vol_[index_i];
    DataType ttl_variable = this->variable_[index_i] * this->Vol_[index_i];
    int neighbor_j = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) *
                                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                                  this->vec_r_ij(index_i, index_j).norm();
        error -= (this->variable_[index_i] - this->variable_[index_j]) * parameter_b[neighbor_j];

        parameter_a += parameter_b[neighbor_j];
        parameter_c += parameter_b[neighbor_j] * parameter_b[neighbor_j];
        ttl_Vol += this->Vol_[index_j];
        ttl_variable += this->variable_[index_j] * this->Vol_[index_j];
        ++neighbor_j;
    }
    parameter_a = (parameter_a - 1.0);

    // update the variable
    DataType parameter_k = error / (parameter_c + parameter_a * parameter_a + TinyReal);
    this->variable_[index_i] += parameter_k * parameter_a;
    neighbor_j = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        this->variable_[index_j] -= parameter_k * parameter_b[neighbor_j];
        ++neighbor_j;
    }
}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::
    PairwiseDissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDissipationType::InteractKernel(ex_policy, encloser),
      capacity_(encloser.dissipation_model_) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    int neighbor_j = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n) // forward sweep
    {
        UnsignedInt index_j = this->neighbor_index_[n];

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        parameter_b[neighbor_j] = this->dis_coeff_(index_i, index_j) * this->dW_ij(index_i, index_j) *
                                  this->Vol_[index_i] * this->Vol_[index_j] * dt /
                                  this->vec_r_ij(index_i, index_j).norm();

        DataType increment = parameter_b[neighbor_j] * pair_difference /
                             (capacity_(index_i) * capacity_(index_j) -
                              parameter_b[neighbor_j] * (capacity_(index_i) + capacity_(index_j)));
        this->variable_[index_i] += increment * capacity_(index_j);
        this->variable_[index_j] -= increment * capacity_(index_i);
        ++neighbor_j;
    }

    --neighbor_j;
    for (UnsignedInt n = this->LastNeighbor(index_i); n != this->FirstNeighbor(index_i); --n) // backward sweep
    {
        UnsignedInt index_j = this->neighbor_index_[n - 1];

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        DataType increment = parameter_b[neighbor_j] * pair_difference /
                             (capacity_(index_i) * capacity_(index_j) -
                              parameter_b[neighbor_j] * (capacity_(index_i) + capacity_(index_j)));
        this->variable_[index_i] += increment * capacity_(index_j);
        this->variable_[index_j] -= increment * capacity_(index_i);
        --neighbor_j;
    }
}
//=================================================================================================//
} // namespace SPH
#endif // SPLITTING_DISSIPATION_HPP
