#ifndef SURFACE_INDICATION_CK_HPP
#define SURFACE_INDICATION_CK_HPP

#include "base_particles.hpp"
#include "surface_indication_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class BaseRelationType>
FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::
    FreeSurfaceIndicationCK(BaseRelationType &base_relation)
    : Interaction<RelationType<Parameters...>>(base_relation),
      dv_indicator_(this->particles_->template registerStateVariable<int>("Indicator")),
      dv_pos_div_(this->particles_->template registerStateVariable<Real>("PositionDivergence")),
      dv_threshold_by_dimensions_(0.75 * Dimensions),
      dv_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   FreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      pos_div_(encloser.dv_pos_div_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      threshold_by_dimensions_(encloser.dv_threshold_by_dimensions_),
      smoothing_length_(encloser.dv_smoothing_length_) {}
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>::
    FreeSurfaceIndicationCK(Inner<Parameters...> &inner_relation)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>(inner_relation),
      dv_previous_surface_indicator_(
          this->particles_->template registerStateVariable<int>("PreviousSurfaceIndicator"))
{
    this->particles_->template addEvolvingVariable<int>("PreviousSurfaceIndicator");
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] = pos_div;
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy,
                 FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
      outer_(&encloser) {}
//=================================================================================================//
template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, Parameters...>>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    // Detect if near surface based on neighbors
    bool is_near_surface = false;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        if ((this->pos_div_[index_j] < this->threshold_by_dimensions_) &&
            (r_ij < this->smoothing_length_))
        {
            is_near_surface = true;
            break;
        }
    }

    // Check if near a previously marked surface
    bool is_near_previous_surface = false;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        if (this->previous_surface_indicator_[index_j] == 1)
        {
            is_near_previous_surface = true;
            break;
        }
    }

    // If the particle’s pos_div is under threshold but isolated (no surface neighbors),
    // push it above threshold to avoid false positives.
    if ((this->pos_div_[index_i] < this->threshold_by_dimensions_) &&
        previous_surface_indicator_[index_i] != 1 && !is_near_previous_surface)
    {
        this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
    }

    // Determine new indicator: surface (1) or not (0)
    int new_indicator = 1;
    if ((this->pos_div_[index_i] > this->threshold_by_dimensions_) && !is_near_surface)
    {
        new_indicator = 0;
    }

    // Update both current and previous indicators
    this->indicator_[index_i] = new_indicator;
    this->previous_surface_indicator_[index_i] = new_indicator;
}
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Contact<Parameters...>>::
    FreeSurfaceIndicationCK(Contact<Parameters...> &contact_relation)
    : FreeSurfaceIndicationCK<Base, Contact<Parameters...>>(contact_relation) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   FreeSurfaceIndicationCK<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : FreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void FreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] += pos_div;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // SURFACE_INDICATION_CK_HPP
