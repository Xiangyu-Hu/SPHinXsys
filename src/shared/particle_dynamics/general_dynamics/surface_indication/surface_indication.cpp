#include "surface_indication.hpp"

namespace SPH
{
//=================================================================================================//
FreeSurfaceIndication<Inner<>>::
    FreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : FreeSurfaceIndication<DataDelegateInner>(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
void FreeSurfaceIndication<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        pos_div -= inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.r_ij_[n];
    }
    pos_div_[index_i] = pos_div;
}
//=================================================================================================//
void FreeSurfaceIndication<Inner<>>::update(size_t index_i, Real dt)
{
    indicator_[index_i] = 1;
    if (pos_div_[index_i] > threshold_by_dimensions_ && !isVeryNearFreeSurface(index_i))
        indicator_[index_i] = 0;
}
//=================================================================================================//
bool FreeSurfaceIndication<Inner<>>::isVeryNearFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        /** Two layer particles.*/
        if (pos_div_[inner_neighborhood.j_[n]] < threshold_by_dimensions_ &&
            inner_neighborhood.r_ij_[n] < smoothing_length_)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
FreeSurfaceIndication<Inner<SpatialTemporal>>::
    FreeSurfaceIndication(BaseInnerRelation &inner_relation)
    : FreeSurfaceIndication<Inner<>>(inner_relation),
      previous_surface_indicator_(particles_->registerStateVariable<int>("PreviousSurfaceIndicator", 1))
{
    particles_->addEvolvingVariable<int>("PreviousSurfaceIndicator");
}
//=================================================================================================//
void FreeSurfaceIndication<Inner<SpatialTemporal>>::interaction(size_t index_i, Real dt)
{
    FreeSurfaceIndication<Inner<>>::interaction(index_i, dt);

    if (pos_div_[index_i] < threshold_by_dimensions_ &&
        previous_surface_indicator_[index_i] != 1 &&
        !isNearPreviousFreeSurface(index_i))
        pos_div_[index_i] = 2.0 * threshold_by_dimensions_;
}
//=================================================================================================//
bool FreeSurfaceIndication<Inner<SpatialTemporal>>::isNearPreviousFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (previous_surface_indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
void FreeSurfaceIndication<Inner<SpatialTemporal>>::update(size_t index_i, Real dt)
{
    FreeSurfaceIndication<Inner<>>::update(index_i, dt);

    previous_surface_indicator_[index_i] = indicator_[index_i];
}
//=================================================================================================//
void FreeSurfaceIndication<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            pos_div -= contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.r_ij_[n];
        }
    }
    pos_div_[index_i] += pos_div;
}
//=================================================================================================//
FreeSurfaceIndication<Contact<NonWetting>>::FreeSurfaceIndication(BaseContactRelation &contact_relation)
    : FreeSurfaceIndication<DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_phi_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("Phi"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void FreeSurfaceIndication<Contact<NonWetting>>::interaction(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *wetting_k = contact_phi_[k];
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            pos_div -= wetting_k[index_j] * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.r_ij_[n];
        }
    }
    pos_div_[index_i] += pos_div;
};
//=================================================================================================//
} // namespace SPH
