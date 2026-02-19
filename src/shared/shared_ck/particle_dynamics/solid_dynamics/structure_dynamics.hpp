#ifndef STRUCTURE_DYNAMICS_HPP
#define STRUCTURE_DYNAMICS_HPP

#include "structure_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
AcousticTimeStepCK::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : h_min_(encloser.h_min_), c0_(encloser.c0_),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Real AcousticTimeStepCK::ReduceKernel::reduce(size_t index_i, Real dt)
{
    Real force_norm = (force_[index_i] + force_prior_[index_i]).norm();
    Real acc_scale = math::sqrt(4.0 * h_min_ * force_norm / mass_[index_i]);
    return SMAX(c0_ + vel_[index_i].norm(), acc_scale);
}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    StructureIntegration1stHalf(Inner<Parameters...> &inner_relation, Real numerical_damping_factor)
    : BaseInteraction(inner_relation), StructureIntegrationVariables(this->particles_),
      material_(DynamicCast<MaterialType>(this, this->particles_->getBaseMaterial())),
      adaptation_(DynamicCast<Adaptation>(this, this->sph_body_->getSPHAdaptation())),
      kernel_correction_(this->particles_), h_ref_(adaptation_.ReferenceSmoothingLength()),
      numerical_damping_factor_(numerical_damping_factor) {}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : constitute_(ex_policy, encloser.material_), correction_(ex_policy, encloser.kernel_correction_),
      rho0_(encloser.material_.ReferenceDensity()), G_(encloser.material_.ShearModulus()),
      h_ref_(encloser.h_ref_), numerical_damping_factor_(encloser.numerical_damping_factor_),
      h_ratio_(ex_policy, encloser.adaptation_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)),
      inverse_F_(encloser.dv_inverse_F_->DelegatedData(ex_policy)),
      dF_dt_(encloser.dv_dF_dt_->DelegatedData(ex_policy)),
      scaling_matrix_(encloser.dv_scaling_matrix_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
void StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Real one_over_J = 1.0 / J;
    rho_[index_i] = rho0_ * one_over_J;

    Matd normalized_be = constitute_.ElasticLeftCauchy(F_[index_i], index_i, dt);
    inverse_F_[index_i] = F_[index_i].inverse();
    Matd inverse_F_T = inverse_F_[index_i].transpose();
    scaling_matrix_[index_i] = normalized_be * inverse_F_T;
    Real isotropic_stress = G_ * normalized_be.trace() * OneOverDimensions;
    // Note that as we use small numerical damping here, the time step size (CFL number) may need to be decreased.
    stress_on_particle_[index_i] =
        inverse_F_T * (constitute_.VolumetricKirchhoff(J) - isotropic_stress) * correction_(index_i) +
        numerical_damping_factor_ *
            constitute_.NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], h_ref_ / h_ratio_(index_i), index_i) *
            inverse_F_T;
}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      G_(encloser.material_.ShearModulus()),
      Vol0_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      scaling_matrix_(encloser.dv_scaling_matrix_->DelegatedData(ex_policy)),
      inverse_F_(encloser.dv_inverse_F_->DelegatedData(ex_policy)),
      stress_on_particle_(encloser.dv_stress_on_particle_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
void StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd sum = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd nablaW_ijV_j = this->nablaW_ij(index_i, index_j) * Vol0_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        Vecd pair_distance = pos_[index_i] - pos_[index_j];
        Matd pair_scaling = scaling_matrix_[index_i] + scaling_matrix_[index_j];
        Matd pair_inverse_F = 0.5 * (inverse_F_[index_i] + inverse_F_[index_j]);
        Vecd e_ij_difference = pair_inverse_F * pair_distance / r_ij - e_ij;
        Real e_ij_difference_norm = e_ij_difference.norm();

        Real limiter = SMIN(10.0 * SMAX(e_ij_difference_norm - 0.05, 0.0), 1.0);
        Vecd corrected_e_ij = e_ij + limiter * e_ij_difference;
        Matd correction_stress = G_ * pair_scaling * corrected_e_ij * corrected_e_ij.transpose();
        sum += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) + correction_stress) * nablaW_ijV_j;
    }

    force_[index_i] = sum * Vol0_[index_i];
};
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
void StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
template <typename... Parameters>
StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::StructureIntegration2ndHalf(
    Inner<Parameters...> &inner_relation)
    : BaseInteraction(inner_relation), StructureIntegrationVariables(this->particles_) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::InitializeKernel::
    InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::InitializeKernel::
    initialize(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol0_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)),
      dF_dt_(encloser.dv_dF_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Matd sum = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd nablaW_ijV_j = this->nablaW_ij(index_i, index_j) * Vol0_[index_j];
        sum -= (vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
    }
    dF_dt_[index_i] = sum * B_[index_i];
};
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : F_(encloser.dv_F_->DelegatedData(ex_policy)),
      dF_dt_(encloser.dv_dF_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateElasticNormalDirectionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : n_(encloser.dv_n_->DelegatedData(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedData(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      phi0_(encloser.dv_phi0_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void UpdateElasticNormalDirectionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    n_[index_i] = polarRotation(F_[index_i]) * n0_[index_i];
    // Nanson's relation is used to update the distance to surface
    Vecd current_normal = F_[index_i].inverse().transpose() * n0_[index_i];
    phi_[index_i] = phi0_[index_i] / (current_normal.norm() + SqrtEps); // todo: check this
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
UpdateAnisotropicMeasure::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : scaling_(encloser.dv_scaling_->DelegatedData(ex_policy)),
      scaling0_(encloser.dv_scaling0_->DelegatedData(ex_policy)),
      orientation_(encloser.dv_orientation_->DelegatedData(ex_policy)),
      orientation0_(encloser.dv_orientation0_->DelegatedData(ex_policy)),
      F_(encloser.dv_F_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void UpdateAnisotropicMeasure::UpdateKernel::update(size_t index_i, Real dt)
{
    Matd rotation = Matd::Identity();
    Matd stretch = Matd::Identity();
    polarDecomposition(F_[index_i], rotation, stretch);
    orientation_[index_i] = rotation * orientation0_[index_i];
    scaling_[index_i] = stretch * scaling0_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // STRUCTURE_DYNAMICS_HPP
