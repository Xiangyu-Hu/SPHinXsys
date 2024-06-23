#pragma once

#include "continuum_integration.h"

#include "base_particles.hpp"
namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class FluidDynamicsType>
BaseIntegration1stHalf<FluidDynamicsType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
    : FluidDynamicsType(inner_relation),
      acc_shear_(*this->particles_->template getVariableDataByName<Vecd>("AccelerationByShear")) {}
//=================================================================================================//
template <class FluidDynamicsType>
void BaseIntegration1stHalf<FluidDynamicsType>::update(size_t index_i, Real dt)
{
    this->vel_[index_i] += ((this->force_prior_[index_i] + this->force_[index_i]) / this->mass_[index_i] + this->acc_shear_[index_i]) * dt;
}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BasePlasticIntegration<DataDelegationType>::BasePlasticIntegration(BaseRelationType &base_relation)
    : fluid_dynamics::BaseIntegration<DataDelegationType>(base_relation),
      plastic_continuum_(DynamicCast<PlasticContinuum>(this, this->particles_->getBaseMaterial())),
      stress_tensor_3D_(*this->particles_->template registerSharedVariable<Mat3d>("StressTensor3D")),
      strain_tensor_3D_(*this->particles_->template registerSharedVariable<Mat3d>("StrainTensor3D")),
      stress_rate_3D_(*this->particles_->template registerSharedVariable<Mat3d>("StressRate3D")),
      strain_rate_3D_(*this->particles_->template registerSharedVariable<Mat3d>("StrainRate3D")),
      elastic_strain_tensor_3D_(*this->particles_->template registerSharedVariable<Mat3d>("ElasticStrainTensor3D")),
      elastic_strain_rate_3D_(*this->particles_->template registerSharedVariable<Mat3d>("ElasticStrainRate3D")),
      velocity_gradient_(*this->particles_->template registerSharedVariable<Matd>("VelocityGradient"))
{
    this->particles_->template addVariableToSort<Mat3d>("ElasticStrainTensor3D");
    this->particles_->template addVariableToSort<Mat3d>("ElasticStrainRate3D");
    this->particles_->template addVariableToSort<Mat3d>("StrainTensor3D");
    this->particles_->template addVariableToSort<Mat3d>("StressTensor3D");
    this->particles_->template addVariableToSort<Mat3d>("StrainRate3D");
    this->particles_->template addVariableToSort<Mat3d>("StressRate3D");
}
//=================================================================================================//
template <class RiemannSolverType>
PlasticIntegration1stHalf<Inner<>, RiemannSolverType>::
    PlasticIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BasePlasticIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(plastic_continuum_, plastic_continuum_) {}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration1stHalf<Inner<>, RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    p_[index_i] = -stress_tensor_3D_[index_i].trace() / 3;
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType>
Vecd PlasticIntegration1stHalf<Inner<>, RiemannSolverType>::computeNonConservativeForce(size_t index_i)
{
    Vecd force = force_prior_[index_i] * rho_[index_i];
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        force += mass_[index_i] * (p_[index_i] - p_[index_j]) * dW_ijV_j * e_ij;
    }
    return force / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration1stHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    Real rho_i = rho_[index_i];
    Matd stress_tensor_i = degradeToMatd(stress_tensor_3D_[index_i]);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        Matd stress_tensor_j = degradeToMatd(stress_tensor_3D_[index_j]);
        force += mass_[index_i] * rho_[index_j] * ((stress_tensor_i + stress_tensor_j) / (rho_i * rho_[index_j])) * nablaW_ijV_j;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
    }
    force_[index_i] += force;
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration1stHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
template <class RiemannSolverType>
PlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>::
    PlasticIntegration1stHalf(BaseContactRelation &wall_contact_relation)
    : BaseIntegrationWithWall(wall_contact_relation),
      riemann_solver_(plastic_continuum_, plastic_continuum_) {}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd force_prior_i = computeNonConservativeForce(index_i);
    Vecd force = force_prior_i;
    Real rho_dissipation(0);
    Matd stress_tensor_i = degradeToMatd(stress_tensor_3D_[index_i]);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &wall_acc_ave_k = *(wall_acc_ave_[k]);
        StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Real face_wall_external_acceleration = (force_prior_i / mass_[index_i] - wall_acc_ave_k[index_j]).dot(-e_ij);
            Real p_in_wall = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            force += 2 * mass_[index_i] * stress_tensor_i * dW_ijV_j * wall_neighborhood.e_ij_[n];
            rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    force_[index_i] += force / rho_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
Vecd PlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>::computeNonConservativeForce(size_t index_i)
{
    return this->force_prior_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>::PlasticIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : BasePlasticIntegration<DataDelegateInner>(inner_relation),
      riemann_solver_(plastic_continuum_, plastic_continuum_, 20.0 * (Real)Dimensions),
      acc_deviatoric_plastic_strain_(*particles_->registerSharedVariable<Real>("AccDeviatoricPlasticStrain")),
      vertical_stress_(*particles_->registerSharedVariable<Real>("VerticalStress")),
      Vol_(*particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableDataByName<Real>("Mass")),
      E_(plastic_continuum_.getYoungsModulus()), nu_(plastic_continuum_.getPoissonRatio())
{
    particles_->addVariableToSort<Real>("AccDeviatoricPlasticStrain");
    particles_->addVariableToSort<Real>("VerticalStress");
}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    Matd velocity_gradient = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += mass_[index_i] * riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
        velocity_gradient -= (vel_[index_i] - vel_[index_j]) * dW_ijV_j * e_ij.transpose();
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation / rho_[index_i];
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
template <class RiemannSolverType>
Real PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>::getDeviatoricPlasticStrain(Mat3d &strain_tensor)
{
    Mat3d deviatoric_strain_tensor = strain_tensor - (1.0 / (Real)Dimensions) * strain_tensor.trace() * Mat3d::Identity();
    Real sum = (deviatoric_strain_tensor.cwiseProduct(deviatoric_strain_tensor)).sum();
    return sqrt(sum * 2.0 / 3.0);
}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    Vol_[index_i] = mass_[index_i] / rho_[index_i];
    Mat3d velocity_gradient = upgradeToMat3d(velocity_gradient_[index_i]);
    Mat3d stress_tensor_rate_3D_ = plastic_continuum_.ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]);
    stress_rate_3D_[index_i] += stress_tensor_rate_3D_;
    stress_tensor_3D_[index_i] += stress_rate_3D_[index_i] * dt;
    /*return mapping*/
    Mat3d stress_tensor_ = plastic_continuum_.ReturnMapping(stress_tensor_3D_[index_i]);
    stress_tensor_3D_[index_i] = stress_tensor_;
    vertical_stress_[index_i] = stress_tensor_3D_[index_i](1, 1);
    strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt;
    /*calculate elastic strain*/
    Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1.0 / 3.0) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
    Real hydrostatic_pressure = (1.0 / 3.0) * stress_tensor_3D_[index_i].trace();
    elastic_strain_tensor_3D_[index_i] = deviatoric_stress / (2.0 * plastic_continuum_.getShearModulus(E_, nu_)) +
                                         hydrostatic_pressure * Mat3d::Identity() / (9.0 * plastic_continuum_.getBulkModulus(E_, nu_));
    Mat3d plastic_strain_tensor_3D = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D_[index_i];
    acc_deviatoric_plastic_strain_[index_i] = getDeviatoricPlasticStrain(plastic_strain_tensor_3D);
}
//=================================================================================================//
template <class RiemannSolverType>
PlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::
    PlasticIntegration2ndHalf(BaseContactRelation &wall_contact_relation)
    : BaseIntegrationWithWall(wall_contact_relation),
      riemann_solver_(plastic_continuum_, plastic_continuum_) {}
//=================================================================================================//
template <class RiemannSolverType>
void PlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    Vecd vel_i = vel_[index_i];
    Matd velocity_gradient = Matd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_[index_i];
            density_change_rate += (vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += mass_[index_i] * riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
            velocity_gradient -= (vel_i - vel_in_wall) * dW_ijV_j * e_ij.transpose();
        }
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation / rho_[index_i];
    velocity_gradient_[index_i] += velocity_gradient;
}
} // namespace continuum_dynamics
} // namespace SPH
