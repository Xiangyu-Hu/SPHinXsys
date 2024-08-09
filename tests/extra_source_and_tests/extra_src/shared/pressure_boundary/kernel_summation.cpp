#include "kernel_summation.hpp"

namespace SPH
{
//=================================================================================================//
NablaWV<Inner<>>::
    NablaWV(BaseInnerRelation &inner_relation)
    : NablaWV<DataDelegateInner>(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
void NablaWV<Inner<>>::interaction(size_t index_i, Real dt)
{
    kernel_sum_[index_i] = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        kernel_sum_[index_i] += inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    }
}
//=================================================================================================//
void NablaWV<Contact<>>::interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            kernel_sum_[index_i] += contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
        }
    }
}
//=================================================================================================//
} // namespace SPH