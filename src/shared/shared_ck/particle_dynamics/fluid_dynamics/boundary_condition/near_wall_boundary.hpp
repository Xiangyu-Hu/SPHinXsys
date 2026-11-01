#ifndef NEAR_WALL_BOUNDARY_HPP
#define NEAR_WALL_BOUNDARY_HPP

#include "near_wall_boundary.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <typename... Parameters>
DistanceFromWallCK<Contact<WithInitialization, Parameters...>>::
    DistanceFromWallCK(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation),
      spacing_ref_(this->getSPHAdaptation().ReferenceSpacing()),
      distance_default_(100.0 * spacing_ref_),
      h_ref_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_n_(this->particles_->template registerStateVariable<Vecd>("NormalDirection")),
      dv_phi_(this->particles_->template registerStateVariable<Real>("SignedDistance", distance_default_)),
      dv_wall_pos_(this->contact_particles_->template getVariableByName<Vecd>("Position")),
      dv_wall_n_(this->contact_particles_->template getVariableByName<Vecd>("NormalDirection")),
      dv_wall_phi_(this->contact_particles_->template getVariableByName<Real>("SignedDistance")) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DistanceFromWallCK<Contact<WithInitialization, Parameters...>>::InitializeKernel::
    InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : distance_default_(encloser.distance_default_), phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DistanceFromWallCK<Contact<WithInitialization, Parameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    phi_[index_i] = -distance_default_; // bulk region far from wall
    n_[index_i] = Vecd::Zero();
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DistanceFromWallCK<Contact<WithInitialization, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      spacing_ref_square_(encloser.spacing_ref_ * encloser.spacing_ref_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      wall_pos_(encloser.dv_wall_pos_->DelegatedData(ex_policy)),
      wall_n_(encloser.dv_wall_n_->DelegatedData(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      wall_phi_(encloser.dv_wall_phi_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DistanceFromWallCK<Contact<WithInitialization, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{

    Vecd normal_sum = Vecd::Zero();
    Real phi_sum = 0.0;
    Real total_weight = TinyReal;
    bool is_near_wall = false;
    UnsignedInt count_neighbor = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
        Real r_ij_square = vec_r_ij.squaredNorm();

        if (r_ij_square < spacing_ref_square_)
        {
            is_near_wall = true;
        }
        Real W_ij = this->W_ij(index_i, index_j);
        const Vecd &wall_normal = wall_n_[index_j];
        Real projection_ij = vec_r_ij.dot(wall_normal);
        Real weight_ij = W_ij * W_ij * projection_ij * projection_ij / r_ij_square;
        normal_sum += -wall_normal * weight_ij;
        phi_sum += -(wall_phi_[index_j] + projection_ij) * weight_ij;
        total_weight += weight_ij;
        count_neighbor++;
    }

    if (is_near_wall && count_neighbor > Dimensions)
    {
        phi_[index_i] = phi_sum / total_weight;
        n_[index_i] = normal_sum.normalized();
    }
}
//=================================================================================================//
template <typename... Parameters>
BoundingFromWallCK<Contact<OneLevel, Parameters...>>::BoundingFromWallCK(
    Contact<Parameters...> &wall_contact_relation)
    : DistanceFromWallCK<Contact<WithInitialization, Parameters...>>(wall_contact_relation),
      distance_min_(0.25 * this->spacing_ref_) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
BoundingFromWallCK<Contact<OneLevel, Parameters...>>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : spacing_ref_(encloser.spacing_ref_),
      distance_min_(encloser.distance_min_),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void BoundingFromWallCK<Contact<OneLevel, Parameters...>>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (phi_[index_i] > -distance_min_)
    {
        pos_[index_i] += -(phi_[index_i] + 0.5 * spacing_ref_) * n_[index_i]; // flip near wall distance
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // NEAR_WALL_BOUNDARY_HPP
