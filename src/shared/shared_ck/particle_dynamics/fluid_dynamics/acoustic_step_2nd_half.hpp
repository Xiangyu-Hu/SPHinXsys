#ifndef ACOUSTIC_STEP_2ND_HALF_HPP
#define ACOUSTIC_STEP_2ND_HALF_HPP

#include "acoustic_step_2nd_half.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(DynamicsIdentifier &identifier)
    : AcousticStep<Interaction<Inner<Parameters...>>>(identifier),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_(ex_policy, encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    force_[index_i] = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        Vecd vel_ave = riemann_.AverageV(index_i, index_j, vel_[index_i], vel_[index_j]);
        compression_rate_[index_i] +=
            2.0 * (vel_[index_i] - vel_ave).dot(correction_(index_i) * e_ij) * dW_ijV_j;
        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        force_[index_i] +=
            riemann_.DissipativePJump(index_i, index_j, u_jump) * dW_ijV_j * Vol_[index_i] * e_ij;
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : eos_(ex_policy, encloser.fluid_), rho_(encloser.dv_rho_->DelegatedDataView(ex_policy)),
      compression_(encloser.dv_compression_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    compression_[index_i] += 0.5 * dt * compression_rate_[index_i];
    rho_[index_i] = compression_[index_i] * eos_.getReferenceDensity(index_i);
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), Interaction<Wall>(identifier),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_(ex_policy, encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataView(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedDataView(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        Vecd vel_diff = 2.0 * (vel_[index_i] - wall_vel_ave_[index_j]);
        compression_rate_[index_i] += vel_diff.dot(correction_(index_i) * e_ij) * dW_ijV_j;
        Vecd face_to_fluid_n = SGN(e_ij.dot(wall_n_[index_j])) * wall_n_[index_j];
        Real u_jump = vel_diff.dot(face_to_fluid_n);
        force_[index_i] += riemann_.DissipativePJump(index_i, index_j, u_jump) * dW_ijV_j *
                           Vol_[index_i] * face_to_fluid_n;
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), kernel_correction_(this->particles_)
{
    SourceFluidType &source_fluid =
        DynamicCast<SourceFluidType>(this, this->sph_body_->getMatterMaterial());
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        TargetFluidType &target_fluid =
            DynamicCast<TargetFluidType>(this, this->contact_bodies_[k]->getMatterMaterial());
        riemann_solvers_.push_back(RiemannSolverType(source_fluid, target_fluid));
        dv_contact_vel_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_(ex_policy, encloser.riemann_solvers_[contact_index]),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataView(ex_policy)),
      contact_vel_(encloser.dv_contact_vel_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        Vecd vel_ave = riemann_.AverageV(index_i, index_j, vel_[index_i], contact_vel_[index_j]);
        compression_rate_[index_i] +=
            2.0 * (vel_[index_i] - vel_ave).dot(correction_(index_i) * e_ij) * dW_ijV_j;
        Real u_jump = (vel_[index_i] - contact_vel_[index_j]).dot(e_ij);
        force_[index_i] +=
            riemann_.DissipativePJump(index_i, index_j, u_jump) * dW_ijV_j * Vol_[index_i] * e_ij;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_2ND_HALF_HPP
