#ifndef SURFACE_INDICATION_CK_HPP
#define SURFACE_INDICATION_CK_HPP

#include "base_particles.hpp"
#include "surface_indication_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//=================================================================================================//
// FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class BaseRelationType>
FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::FreeSurfaceIndicationCK(BaseRelationType &base_relation)
    : Interaction<RelationType<Parameters...>>(base_relation),
      dv_indicator_(this->particles_->template registerStateVariableOnly<int>("Indicator")),
      dv_pos_div_(this->particles_->template registerStateVariableOnly<Real>("PositionDivergence")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_previous_surface_indicator_(
          this->particles_->template registerStateVariableOnly<int>("PreviousSurfaceIndicator")),
      dv_threshold_by_dimensions_(0.75 * Dimensions),
      dv_smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength())
{
}

template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Base, RelationType<Parameters...>> &encloser,
    Args &&...args)
    : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      pos_div_(encloser.dv_pos_div_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      threshold_by_dimensions_(encloser.dv_threshold_by_dimensions_),
      smoothing_length_(encloser.dv_smoothing_length_)
{
}

template <template <typename...> class RelationType, typename... Parameters>
void FreeSurfaceIndicationCK<Base, RelationType<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] = pos_div;
}

//=================================================================================================//
// FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>::FreeSurfaceIndicationCK(
    Relation<Inner<Parameters...>> &inner_relation)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>(inner_relation),
      dv_previous_surface_indicator_(
          this->particles_->template getVariableByName<int>("PreviousSurfaceIndicator"))
{
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>::InteractKernel::interact(
    size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] = pos_div;
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, Internal, Parameters...>>::UpdateKernel::update(
    size_t index_i, Real dt)
{
    // int new_indicator = 1;
    if (this->pos_div_[index_i] > this->threshold_by_dimensions_)
    {
        this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
    }
    // this->indicator_[index_i] = new_indicator;
    // this->previous_surface_indicator_[index_i] = new_indicator;
}

//=================================================================================================//
// FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>>
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>>::FreeSurfaceIndicationCK(
    Relation<Inner<Parameters...>> &inner_relation)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>(inner_relation)
{
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser)
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<Base, Internal, Parameters...>>::InteractKernel::interact(
    size_t index_i, Real dt)
{
    // Execute the base interaction without any update.
    FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel::interact(index_i, dt);
}

//=================================================================================================//
// FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>::FreeSurfaceIndicationCK(
    Relation<Inner<Parameters...>> &inner_relation)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>(inner_relation),
      dv_previous_surface_indicator_(
          this->particles_->template getVariableByName<int>("PreviousSurfaceIndicator")),
      dv_is_near_surface_indicator_(
          this->particles_->template registerStateVariableOnly<int>("IsNearSurfaceIndicator"))
{
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
      is_near_surface_indicator_(
          encloser.dv_is_near_surface_indicator_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>::InteractKernel::interact(
    size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] = pos_div;
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>> &encloser)
    : FreeSurfaceIndicationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
      is_near_surface_indicator_(
          encloser.dv_is_near_surface_indicator_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Inner<WithUpdate, SpatialTemporal, Parameters...>>::UpdateKernel::update(
    size_t index_i, Real dt)
{
    bool is_near_surface = false;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        /** Two layer particles.*/
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        if ((this->pos_div_[index_j] < this->threshold_by_dimensions_) &&
            (r_ij < this->smoothing_length_))
        {
            is_near_surface = true;
            break;
        }
    }
    this->is_near_surface_indicator_[index_i] = is_near_surface;

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

    if (this->pos_div_[index_i] < this->threshold_by_dimensions_ &&
        this->previous_surface_indicator_[index_i] != 1 &&
        !is_near_previous_surface)
    {
        this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
    }

    int new_indicator = 1;
    if ((this->pos_div_[index_i] > this->threshold_by_dimensions_) &&
        !is_near_surface)
    {
        new_indicator = 0;
    }
    this->indicator_[index_i] = new_indicator;
    this->previous_surface_indicator_[index_i] = new_indicator;
}

//=================================================================================================//
// FreeSurfaceIndicationCK<Contact<Parameters...>>
//=================================================================================================//
template <typename... Parameters>
FreeSurfaceIndicationCK<Contact<Parameters...>>::FreeSurfaceIndicationCK(
    Relation<Contact<Parameters...>> &contact_relation)
    : FreeSurfaceIndicationCK<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k < this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy>
FreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy,
    FreeSurfaceIndicationCK<Contact<Parameters...>> &encloser,
    size_t contact_index)
    : FreeSurfaceIndicationCK<Base, Contact<Parameters...>>::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
void FreeSurfaceIndicationCK<Contact<Parameters...>>::InteractKernel::interact(
    size_t index_i, Real dt)
{
    Real pos_div = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        pos_div -= this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j] * r_ij;
    }
    this->pos_div_[index_i] += pos_div;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
SurfaceIndicationByAlignedBoxCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy))
{
}

void SurfaceIndicationByAlignedBoxCK::UpdateKernel::update(size_t index_i, Real dt)
{
    int indicator = 0;
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        indicator = 1;
    }
    indicator_[index_i] = indicator;
}
} // namespace fluid_dynamics
} // namespace SPH

#endif // SURFACE_INDICATION_CK_HPP
