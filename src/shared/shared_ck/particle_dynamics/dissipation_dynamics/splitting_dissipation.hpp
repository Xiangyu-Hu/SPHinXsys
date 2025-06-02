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
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        UnsignedInt neighbor_j = n - this->FirstNeighbor(index_i);

        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) *
                                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt /
                                  this->vec_r_ij(index_i, index_j).norm();
        error -= (this->variable_[index_i] - this->variable_[index_j]) * parameter_b[neighbor_j];

        parameter_a += parameter_b[neighbor_j];
        parameter_c += parameter_b[neighbor_j] * parameter_b[neighbor_j];
    }
    parameter_a = (parameter_a - 1.0);

    // update the variable
    DataType parameter_k = error / (parameter_c + parameter_a * parameter_a + TinyReal);
    this->variable_[index_i] += parameter_k * parameter_a;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        UnsignedInt neighbor_j = n - this->FirstNeighbor(index_i);

        this->variable_[index_j] -= parameter_k * parameter_b[neighbor_j];
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
      inverse_capacity_(encloser.dissipation_model_) {}
//=================================================================================================//
template <typename DissipationType, typename... Parameters>
void PairwiseDissipation<Inner<Splitting, DissipationType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Real split_dt = dt * 0.5;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    // forward sweep
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        UnsignedInt neighbor_j = n - this->FirstNeighbor(index_i);

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        parameter_b[neighbor_j] = 2.0 * this->dis_coeff_(index_i, index_j) * this->dW_ij(index_i, index_j) *
                                  this->Vol_[index_i] * this->Vol_[index_j] * split_dt /
                                  this->vec_r_ij(index_i, index_j).norm();

        DataType increment = 0.5 * parameter_b[neighbor_j] * pair_difference /
                             (this->Vol_[index_i] * this->Vol_[index_j] -
                              parameter_b[neighbor_j] * (this->Vol_[index_i] + this->Vol_[index_j]));
        this->variable_[index_i] += increment * inverse_capacity_(index_i) * this->Vol_[index_j];
        this->variable_[index_j] -= increment * inverse_capacity_(index_j) * this->Vol_[index_i];
    }
    // backward sweep
    for (UnsignedInt n = this->LastNeighbor(index_i); n != this->FirstNeighbor(index_i); --n)
    {
        UnsignedInt index_j = this->neighbor_index_[n - 1];
        UnsignedInt neighbor_j = n - 1 - this->FirstNeighbor(index_i);

        DataType pair_difference = (this->variable_[index_i] - this->variable_[index_j]);
        DataType increment = 0.5 * parameter_b[neighbor_j] * pair_difference /
                             (this->Vol_[index_i] * this->Vol_[index_j] -
                              parameter_b[neighbor_j] * (this->Vol_[index_i] + this->Vol_[index_j]));
        this->variable_[index_i] += increment * inverse_capacity_(index_i) * this->Vol_[index_j];
        this->variable_[index_j] -= increment * inverse_capacity_(index_j) * this->Vol_[index_i];
    }
}
//=================================================================================================//
} // namespace SPH
#endif // SPLITTING_DISSIPATION_HPP
