#ifndef GEOMETRIC_DYNAMICS_HPP
#define GEOMETRIC_DYNAMICS_HPP

#include "execution_policy.h"
#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NormalFromBodyShapeCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, NormalFromBodyShapeCK &encloser)
    : HostKernel(ex_policy, encloser),
      initial_shape_(encloser.initial_shape_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedData(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      phi0_(encloser.dv_phi0_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
SurfaceIndicationFromBodyShape::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, SurfaceIndicationFromBodyShape &encloser)
    : HostKernel(ex_policy, encloser),
      initial_shape_(encloser.initial_shape_),
      spacing_ref_(encloser.spacing_ref_),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
NormalDirectionFromParticlesCK<Base, RelationType<Parameters...>>::
    NormalDirectionFromParticlesCK(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      sph_body_(&identifier.getSPHBody()),
      initial_shape_(&sph_body_->getInitialShape()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_n_(this->particles_->template registerStateVariableOnly<Vecd>("NormalDirection")),
      dv_n0_(this->particles_->template registerStateVariableOnly<Vecd>("InitialNormalDirection", dv_n_)),
      dv_phi_(this->particles_->template registerStateVariableOnly<Real>("SignedDistance")),
      dv_phi0_(this->particles_->template registerStateVariableOnly<Real>("InitialSignedDistance", dv_phi_)),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
NormalDirectionFromParticlesCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   NormalDirectionFromParticlesCK<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      initial_shape_(encloser.initial_shape_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedData(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      phi0_(encloser.dv_phi0_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
NormalDirectionFromParticlesCK<Inner<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   NormalDirectionFromParticlesCK<Inner<Parameters...>> &encloser)
    : NormalDirectionFromParticlesCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <typename... Parameters>
void NormalDirectionFromParticlesCK<Inner<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd normal_direction = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        normal_direction -= this->dW_ij(index_i, index_j)  * this->Vol_[index_j] * this->e_ij(index_i, index_j) ;
    }
    normal_direction = normal_direction / (normal_direction.norm() + TinyReal);
    this->n_[index_i] = normal_direction;
    this->n0_[index_i] = normal_direction;
    Real signed_distance = this->initial_shape_->findSignedDistance(this->pos_[index_i]);
    this->phi_[index_i] = signed_distance;
    this->phi0_[index_i] = signed_distance;
}
//=================================================================================================//
} // namespace SPH
#endif // GEOMETRIC_DYNAMICS_HPP
