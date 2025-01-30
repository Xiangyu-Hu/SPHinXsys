#include "kernel_summation_ck.h"
namespace SPH
{

namespace fluid_dynamics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
NablaWV<Base, RelationType<Parameters...>>::
    NablaWV(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_kernel_sum_(this->particles_->template registerStateVariableOnly<Vecd>("KernelSummation"))
{
}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
NablaWV<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   NablaWV<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      kernel_sum_(encloser.dv_kernel_sum_->DelegatedData(ex_policy))
{
}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
NablaWV<Inner<RegularizationType, Parameters...>>::
    NablaWV(Relation<Inner<Parameters...>> &inner_relation)
    : NablaWV<Base, Inner<Parameters...>>(inner_relation),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
template <class ExecutionPolicy>
NablaWV<Inner<RegularizationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   NablaWV<Inner<RegularizationType, Parameters...>> &encloser)
    : NablaWV<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
void NablaWV<Inner<RegularizationType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd kernel_sum = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        kernel_sum += dW_ijV_j * e_ij;
    }
    this->kernel_sum_[index_i] = kernel_sum;
}
//=================================================================================================//
template <typename... Parameters>
NablaWV<Contact<Parameters...>>::
    NablaWV(Relation<Contact<Parameters...>> &contact_relation)
    : NablaWV<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        dv_contact_wall_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
NablaWV<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   NablaWV<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : NablaWV<Base, Contact<Parameters...>>::
          InteractKernel(ex_policy, encloser, contact_index),
      contact_wall_Vol_(encloser.dv_contact_wall_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void NablaWV<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd kernel_sum = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_wall_Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        kernel_sum += dW_ijV_j * e_ij;
    }
    // std::cout << "kernel_sum: " << kernel_sum << std::endl;
    this->kernel_sum_[index_i] += kernel_sum;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH