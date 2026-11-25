/**
 * @file spdb_contact_penalty_ck.h
 * @brief CK solid–solid contact along e_ij (ContactForce formula), for SYCL/GPU.
 * @details Library RepulsionForceCK uses average normals and fails on flush faces.
 *          This matches host ContactForce: sigma_star along e_ij / nablaW, no normals.
 *          Extra-test only; does not modify src/.
 */
#ifndef SPDB_CONTACT_PENALTY_CK_H
#define SPDB_CONTACT_PENALTY_CK_H

#include "force_prior_ck.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace solid_dynamics
{

template <typename...>
class ContactPenaltyForceCK;

template <typename... Parameters>
class ContactPenaltyForceCK<Base, Contact<Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
  public:
    template <class DynamicsIdentifier>
    explicit ContactPenaltyForceCK(DynamicsIdentifier &identifier)
        : Interaction<Contact<Parameters...>>(identifier),
          ForcePriorCK(this->particles_, "RepulsionForce" + identifier.Name()),
          solid_contact_(DynamicCast<SolidContact>(this, this->sph_body_->getMatterMaterial())),
          stiffness_(solid_contact_.ContactStiffness()),
          dv_repulsion_factor_(this->particles_->template getVariableByName<Real>(
              "RepulsionFactor" + identifier.Name())),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_repulsion_force_(ForcePriorCK::getCurrentForce()) {}
    virtual ~ContactPenaltyForceCK() {};

  protected:
    SolidContact &solid_contact_;
    Real stiffness_;
    DiscreteVariable<Real> *dv_repulsion_factor_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_repulsion_force_;
};

template <typename... Parameters>
class ContactPenaltyForceCK<Contact<WithUpdate, Parameters...>>
    : public ContactPenaltyForceCK<Base, Contact<Parameters...>>
{
    using BaseInteractionType = ContactPenaltyForceCK<Base, Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier, typename... Args>
    ContactPenaltyForceCK(DynamicsIdentifier &identifier, Args &&...args)
        : BaseInteractionType(identifier, std::forward<Args>(args)...)
    {
        Solid &solid = DynamicCast<Solid>(this, this->contact_body_->getMatterMaterial());
        contact_stiffness_ = solid.ContactStiffness();
        dv_contact_repulsion_factor_ =
            this->contact_particles_->template getVariableByName<Real>(
                "RepulsionFactor" + identifier.ReverseName());
    }
    virtual ~ContactPenaltyForceCK() {};

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseInteractionType::InteractKernel(ex_policy, encloser),
              stiffness_ave_(2.0 * encloser.stiffness_ * encloser.contact_stiffness_ /
                             (encloser.stiffness_ + encloser.contact_stiffness_)),
              repulsion_factor_(encloser.dv_repulsion_factor_->DelegatedData(ex_policy)),
              contact_repulsion_factor_(encloser.dv_contact_repulsion_factor_->DelegatedData(ex_policy)),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              contact_Vol_(encloser.dv_contact_Vol_->DelegatedData(ex_policy)),
              repulsion_force_(encloser.dv_repulsion_force_->DelegatedData(ex_policy)) {}

        void interact(size_t index_i, Real dt = 0.0)
        {
            Vecd force = Vecd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                UnsignedInt index_j = this->neighbor_index_[n];
                // Host ContactForce: 0.5*(sigma_i+sigma_j) * stiffness_ave * e_ij * dW * Vol_j
                // with nablaW = dW * e_ij, so e_ij*dW == nablaW.
                const Real sigma_star = 0.5 * (repulsion_factor_[index_i] + contact_repulsion_factor_[index_j]);
                force -= 2.0 * stiffness_ave_ * sigma_star * this->nablaW_ij(index_i, index_j) *
                         contact_Vol_[index_j];
            }
            repulsion_force_[index_i] = force * Vol_[index_i];
        }

      protected:
        Real stiffness_ave_;
        Real *repulsion_factor_, *contact_repulsion_factor_;
        Real *Vol_, *contact_Vol_;
        Vecd *repulsion_force_;
    };

  protected:
    Real contact_stiffness_;
    DiscreteVariable<Real> *dv_contact_repulsion_factor_;
};

} // namespace solid_dynamics
} // namespace SPH

#endif // SPDB_CONTACT_PENALTY_CK_H
