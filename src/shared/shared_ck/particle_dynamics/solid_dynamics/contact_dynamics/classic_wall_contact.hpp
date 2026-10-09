// SPDX-License-Identifier: Apache-2.0
#ifndef CLASSIC_WALL_CONTACT_HPP
#define CLASSIC_WALL_CONTACT_HPP

#include "classic_wall_contact.h"

#include "base_body_part.h"
#include "base_particles.hpp"
#include "kernel_wendland_c2.h"

#include <stdexcept>
#include <typeinfo>

namespace SPH
{
namespace solid_dynamics
{
template <typename... Parameters>
template <class DynamicsIdentifier>
ClassicWallContactForceCK<Contact<WithUpdate, Parameters...>>::
    ClassicWallContactForceCK(DynamicsIdentifier &identifier, BodySurfaceLayer &surface)
    : BaseInteractionType(identifier),
      ForcePriorCK(this->particles_, "RepulsionForce" + identifier.Name()),
      stiffness_(DynamicCast<SolidContact>(this, this->sph_body_->getMatterMaterial()).ContactStiffness())
{
    if (&surface.getSPHBody() != this->sph_body_)
        throw std::invalid_argument("Classic wall contact surface must belong to the source body");
    for (SPHAdaptation *adaptation : {this->sph_adaptation_, this->contact_adaptation_})
        if (typeid(*adaptation) != typeid(SPHAdaptation) ||
            typeid(*adaptation->getKernel()) != typeid(KernelWendlandC2))
            throw std::invalid_argument("Classic wall contact requires uniform isotropic Wendland C2 adaptations");
    auto *wall_count = this->particles_->template registerSingleVariable<int>("ClassicWallContactCount", 0);
    if (wall_count->getValue() != 0)
        throw std::invalid_argument("Classic wall contact supports only one wall per source body");
    *wall_count->Data() = 1;
    dv_contact_factor_ = this->particles_->template registerStateVariable<Real>("ClassicWallContactFactor");
    dv_surface_mask_ = this->particles_->template registerStateVariable<int>("ClassicWallContactSurface");
    int *surface_mask = dv_surface_mask_->Data();
    for (size_t i = 0; i < this->particles_->TotalRealParticles(); ++i)
        surface_mask[i] = 0;
    for (size_t i : surface.LoopRange())
        surface_mask[i] = 1;
    this->particles_->template addEvolvingVariable<int>(dv_surface_mask_);

    const Real h = Real(0.5) * (this->sph_adaptation_->ReferenceSmoothingLength() +
                               this->contact_adaptation_->ReferenceSmoothingLength());
    const Real spacing = Real(0.5) * (this->sph_adaptation_->ReferenceSpacing() +
                                     this->contact_adaptation_->ReferenceSpacing());
    const KernelWendlandC2 kernel(h);
    const Vecd zero = Vecd::Zero();
    inv_h_ = Real(1) / h;
    factor_W_ = kernel.W0(zero);
    factor_dW_ = inv_h_ * factor_W_;
    offset_W_ = kernel.W(spacing, zero);
    cutoff_ = kernel.CutOffRadius();
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ClassicWallContactForceCK<Contact<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteractionType::InteractKernel(ex_policy, encloser),
      inv_h_(encloser.inv_h_), factor_W_(encloser.factor_W_), factor_dW_(encloser.factor_dW_),
      offset_W_(encloser.offset_W_), cutoff_(encloser.cutoff_), stiffness_(encloser.stiffness_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_->DelegatedData(ex_policy)),
      contact_factor_(encloser.dv_contact_factor_->DelegatedData(ex_policy)),
      contact_force_(encloser.getCurrentForce()->DelegatedData(ex_policy)),
      surface_mask_(encloser.dv_surface_mask_->DelegatedData(ex_policy)) {}

template <typename... Parameters>
void ClassicWallContactForceCK<Contact<WithUpdate, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    if (surface_mask_[index_i] == 0)
    {
        contact_factor_[index_i] = 0;
        contact_force_[index_i] = Vecd::Zero();
        return;
    }
    Real sigma = 0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real distance = this->vec_r_ij(index_i, index_j).norm();
        if (distance < cutoff_)
        {
            const Real q = distance * inv_h_;
            const Real t = Real(1) - Real(0.5) * q;
            const Real t2 = t * t;
            const Real W = factor_W_ * (t2 * t2 * (Real(1) + Real(2) * q));
            sigma += SMAX(W - offset_W_, Real(0)) * contact_Vol_[index_j];
        }
    }
    contact_factor_[index_i] = sigma;
    const Real pressure = sigma * stiffness_;
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd displacement = this->vec_r_ij(index_i, index_j);
        const Real distance = displacement.norm();
        if (distance < cutoff_)
        {
            const Real q = distance * inv_h_;
            const Real t = q - Real(2);
            const Real dW = factor_dW_ * (Real(0.625) * t * t * t * q);
            const Vecd direction = displacement / (distance + TinyReal);
            force -= Real(2) * pressure * direction * dW * contact_Vol_[index_j];
        }
    }
    contact_force_[index_i] = force * Vol_[index_i];
}
} // namespace solid_dynamics
} // namespace SPH
#endif // CLASSIC_WALL_CONTACT_HPP
