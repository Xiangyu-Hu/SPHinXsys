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
      variable_(encloser.dv_variable_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DissipationType, typename SourceType, typename... Parameters>
Dissipation<Inner<Splitting, DissipationType, SourceType, Parameters...>>::
    Dissipation(Inner<Parameters...> &inner_relation, const std::string &variable_name)
    : BaseDissipationType(inner_relation, variable_name), source_model_(this->particles_) {}
//=================================================================================================//
template <typename DissipationType, typename SourceType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Dissipation<Inner<Splitting, DissipationType, SourceType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : Dissipation<Base, DissipationType, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      source_(ex_policy, encloser.source_model_) {}
//=================================================================================================//
template <typename DissipationType, typename SourceType, typename... Parameters>
void Dissipation<Inner<Splitting, DissipationType, SourceType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    // compute the error and parameters
    ErrorAndParameters<DataType> error_and_parameters;
    std::array<Real, MaximumNeighborhoodSize> parameter_b;
    int neighbor_k = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        DataType variable_diff = this->variable_[index_i] - this->variable_[index_j];

        // linear projection
        parameter_b[neighbor_k] = 2.0 * this->dis_coeff_(index_i, index_j) *
                                  this->dW_ij(index_i, index_j) * this->Vol_[index_j] * dt / r_ij;

        error_and_parameters.error_ -= variable_diff * parameter_b[neighbor_k];
        error_and_parameters.a_ += parameter_b[neighbor_k];
        error_and_parameters.c_ += parameter_b[neighbor_k] * parameter_b[neighbor_k];
        ++neighbor_k;
    }
    error_and_parameters.a_ -= 1.0;
    error_and_parameters.error_ -= source_(index_i) * dt;

    // update the variable
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    DataType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;

    neighbor_k = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        this->variable_[index_j] -= parameter_k * parameter_b[neighbor_k];
        ++neighbor_k;
    }
}
//=================================================================================================//
} // namespace SPH
#endif // SPLITTING_DISSIPATION_HPP
