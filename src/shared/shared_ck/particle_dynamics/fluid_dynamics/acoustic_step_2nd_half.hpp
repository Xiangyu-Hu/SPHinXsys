#ifndef ACOUSTIC_STEP_2ND_HALF_HPP
#define ACOUSTIC_STEP_2ND_HALF_HPP

#include "acoustic_step_2nd_half.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Relation<Inner<Parameters...>> &inner_relation)
    : AcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
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
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(corrected_e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * corrected_e_ij;
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation * Vol_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Relation<Contact<Parameters...>> &wall_contact_relation)
    : AcousticStep<Interaction<Contact<Wall, Parameters...>>>(wall_contact_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedDataField(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedDataField(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Vecd vel_in_wall = 2.0 * wall_vel_ave_[index_j] - vel_[index_i];
        density_change_rate += (vel_[index_i] - vel_in_wall).dot(corrected_e_ij) * dW_ijV_j;
        Real u_jump = 2.0 * (vel_[index_i] - wall_vel_ave_[index_j]).dot(wall_n_[index_j]);
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * wall_n_[index_j];
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];
}
//=================================================================================================//


//add
//step2 - innner
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStep2ndHalf(Relation<Inner<Parameters...>> &inner_relation)
    : PlasticAcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      velocity_gradient_(encloser.dv_velocity_gradient_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();

    Matd velocity_gradient = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(corrected_e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * corrected_e_ij;

        //add
        velocity_gradient -= (vel_[index_i] - vel_[index_j]) * dW_ijV_j * corrected_e_ij.transpose();
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation * Vol_[index_i];
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedDataField(ex_policy)), 
      strain_tensor_3D_(encloser.dv_strain_tensor_3D_->DelegatedDataField(ex_policy)),
      stress_rate_3D_(encloser.dv_stress_rate_3D_->DelegatedDataField(ex_policy)),
      strain_rate_3D_(encloser.dv_strain_rate_3D_->DelegatedDataField(ex_policy)),
      plastic_func_(encloser.plastic_continuum_)
      {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;


    //add
    Mat3d velocity_gradient = upgradeToMat3d_sycl(velocity_gradient_[index_i]);
    Mat3d stress_tensor_rate_3D_ = plastic_func_.ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]);
    stress_rate_3D_[index_i] += stress_tensor_rate_3D_; //ERROR

    //--------------------------------------Wrong Here--------------------------------------------------------//
    //I can't find the reason
    //If stress_rate_3D_[index_i] = Mat3d::Zero(); It is correct
    //if Real stress_tensor_rate_3D_ = 
    //  plastic_func_.modified_ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]); 
    //  rho_[index_i] = stress_tensor_rate_3D_; also correct
    //--------------------------------------Wrong Here--------------------------------------------------------//

    //--------------------------------------Error Output--------------------------------------------------------//
    /*
    UR CUDA ERROR:  
	Value:           700  
	Name:            CUDA_ERROR_ILLEGAL_ADDRESS  
	Description:     an illegal memory access was encountered  
	Function:        wait  
	Source Location: /tmp/tmp.7vgJ2wJCWQ/intel-llvm-mirror/sycl/plugins/unified_runtime/ur/adapters/cuda/event.cpp:140  
        terminate called after throwing an instance of 'sycl::_V1::runtime_error'  
        what():  Native API failed. Native API returns: -999 (Unknown PI error) -999 (Unknown PI error)  
        Aborted (core dumped) 
    */
    //--------------------------------------Error Output--------------------------------------------------------//
    
    // stress_tensor_3D_[index_i] += stress_rate_3D_[index_i] * dt;
    // /*return mapping*/
    // stress_tensor_3D_[index_i] = plastic_func_.ReturnMapping(stress_tensor_3D_[index_i]);
    // strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    // strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt;
}
//step2 - wall
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStep2ndHalf(Relation<Contact<Parameters...>> &wall_contact_relation)
    : PlasticAcousticStep<Interaction<Contact<Wall, Parameters...>>>(wall_contact_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedDataField(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedDataField(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedDataField(ex_policy)),
      velocity_gradient_(encloser.dv_velocity_gradient_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    Matd velocity_gradient = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Vecd vel_in_wall = 2.0 * wall_vel_ave_[index_j] - vel_[index_i];
        density_change_rate += (vel_[index_i] - vel_in_wall).dot(corrected_e_ij) * dW_ijV_j;
        Real u_jump = 2.0 * (vel_[index_i] - wall_vel_ave_[index_j]).dot(wall_n_[index_j]);
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * wall_n_[index_j];

        //add
        velocity_gradient -= (vel_[index_i] - vel_in_wall) * dW_ijV_j * corrected_e_ij.transpose();
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];

    //add
    velocity_gradient_[index_i] += velocity_gradient;
}
//=================================================================================================//

} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_2ND_HALF_HPP
