#include "kernel_summation.hpp"

namespace SPH
{
//=================================================================================================//
NablaWV<Inner<>>::
    NablaWV(BaseInnerRelation &inner_relation)
    : NablaWV<GeneralDataDelegateInner>(inner_relation) {}
//=================================================================================================//
void NablaWV<Inner<>>::interaction(size_t index_i, Real dt)
{
    kernel_sum_[index_i] = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        kernel_sum_[index_i] += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
    }
}
//=================================================================================================//
void NablaWV<Contact<>>::interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            kernel_sum_[index_i] += contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
    }
}
//=================================================================================================//
} // namespace SPH