#include "non_newtonian_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
Oldroyd_BIntegration1stHalf<Inner<>>::
    Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalfInnerRiemann(inner_relation)
{
    particles_->registerVariable(tau_, "ElasticStress");
    particles_->registerVariable(dtau_dt_, "ElasticStressChangeRate");
    particles_->registerSortableVariable<Matd>("ElasticStress");
    particles_->addVariableToRestart<Matd>("ElasticStress");
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Inner<>>::initialization(size_t index_i, Real dt)
{
    Integration1stHalfInnerRiemann::initialization(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Inner<>>::interaction(size_t index_i, Real dt)
{
    Integration1stHalfInnerRiemann::interaction(index_i, dt);

    Vecd force = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        // elastic force
        force += mass_[index_i] * (tau_[index_i] + tau_[index_j]) * nablaW_ijV_j;
    }

    force_[index_i] += force / rho_[index_i];
}
//=================================================================================================//
Oldroyd_BIntegration1stHalf<Contact<Wall>>::
    Oldroyd_BIntegration1stHalf(BaseContactRelation &wall_contact_relation)
    : Integration1stHalfContactWallRiemann(wall_contact_relation),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")){};
//=================================================================================================//
void Oldroyd_BIntegration1stHalf<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Integration1stHalfContactWallRiemann::interaction(index_i, dt);

    Real rho_i = rho_[index_i];
    Matd tau_i = tau_[index_i];

    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            /** stress boundary condition. */
            force += mass_[index_i] * 2.0 * tau_i * nablaW_ijV_j / rho_i;
        }
    }

    force_[index_i] += force;
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalf<Inner<>>::
    Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : Integration2ndHalfInnerRiemann(inner_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf<Inner<>>::update(size_t index_i, Real dt)
{
    Integration2ndHalfInnerRiemann::update(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf<Inner<>>::interaction(size_t index_i, Real dt)
{
    Integration2ndHalfInnerRiemann::interaction(index_i, dt);

    Matd tau_i = tau_[index_i];
    Matd stress_rate = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        Matd velocity_gradient = -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
        stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
                       (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
    }

    dtau_dt_[index_i] = stress_rate;
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalf<Contact<Wall>>::
    Oldroyd_BIntegration2ndHalf(BaseContactRelation &wall_contact_relation)
    : Integration2ndHalfContactWallRiemann(wall_contact_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Integration2ndHalfContactWallRiemann::interaction(index_i, dt);

    Vecd vel_i = vel_[index_i];
    Matd tau_i = tau_[index_i];

    Matd stress_rate = Matd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];

            Matd velocity_gradient = -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
            stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
                           (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
        }
    }
    dtau_dt_[index_i] += stress_rate;
}

//=================================================================================================//
// Herschel_Bulkley
//=================================================================================================//
HerschelBulkleyAcceleration<Inner<>>::HerschelBulkleyAcceleration(BaseInnerRelation &inner_relation)
    : HerschelBulkleyAcceleration<FluidDataInner>(inner_relation),
      ForcePrior(&base_particles_, "ViscousForce"),
      mu_srd_(*particles_->getVariableByName<Real>("SRDViscosity"))
{
}

void HerschelBulkleyAcceleration<Inner<>>::interaction(size_t index_i, Real dt)
{
    //* Force Calculation *//

    Vecd force = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        //* Monaghan 2005 (Rep. Prog. Phys.) with averaged viscosity

        Real v_r_ij = (vel_[index_i] - vel_[index_j]).dot(r_ij * e_ij);
        Real avg_visc = 2.0 * (mu_srd_[index_i] * mu_srd_[index_j]) / (mu_srd_[index_i] + mu_srd_[index_j]);
        Real eta_ij = 8.0 * avg_visc * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
        force += eta_ij * mass_[index_i] * inner_neighborhood.dW_ijV_j_[n] * e_ij;
    }
    viscous_force_[index_i] = force / rho_[index_i];
}
//=================================================================================================//
HerschelBulkleyAcceleration<Contact<Wall>>::HerschelBulkleyAcceleration(BaseContactRelation &contact_relation)
    : BaseHerschelBulkleyAccelerationWithWall(contact_relation),
      ForcePrior(&base_particles_, "ViscousForce"),
      mu_srd_(*particles_->getVariableByName<Real>("SRDViscosity"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
    }
}

void HerschelBulkleyAcceleration<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    //* Force Calculation *//
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        //! whats the difference?
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        // StdLargeVec<Vecd> &vel_k = *(wall_vel_ave_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];

            Real v_r_ij = (vel_[index_i] - vel_k[index_j]).dot(r_ij * e_ij);
            Real eta_ij = 2.0 * mu_srd_[index_i] * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
            force += eta_ij * mass_[index_i] * wall_neighborhood.dW_ijV_j_[n] * e_ij;
        }
    }
    viscous_force_[index_i] += force / rho_[index_i];
}
//=================================================================================================//
VelocityGradientInner::VelocityGradientInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      FluidDataInner(inner_relation),
      vel_(particles_->vel_)
{
    particles_->registerVariable(combined_velocity_gradient_, "CombinedVelocityGradient");
    particles_->addVariableToWrite<Matd>("CombinedVelocityGradient");
}
//=================================================================================================//
void VelocityGradientInner::interaction(size_t index_i, Real dt)
{
    //? calculate velocity gradient
    combined_velocity_gradient_[index_i] = Matd::Zero();
    Matd velocity_gradient = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        velocity_gradient += -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
    }
    combined_velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
VelocityGradientContact::VelocityGradientContact(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      FluidContactData(contact_relation),
      vel_(particles_->vel_),
      combined_velocity_gradient_(*particles_->getVariableByName<Matd>("CombinedVelocityGradient"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
    }
}
//=================================================================================================//
void VelocityGradientContact::interaction(size_t index_i, Real dt)
{
    //? calculate velocity gradient
    Matd velocity_gradient = Matd::Zero();

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            velocity_gradient += -(vel_[index_i] - vel_k[index_j]) * nablaW_ijV_j.transpose();
        }
    }
    combined_velocity_gradient_[index_i] += velocity_gradient;
}
//=================================================================================================//
ShearRateDependentViscosity::ShearRateDependentViscosity(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      FluidDataInner(inner_relation),
      combined_velocity_gradient_(*particles_->getVariableByName<Matd>("CombinedVelocityGradient")),
      herschel_bulkley_fluid_(DynamicCast<HerschelBulkleyFluid>(this, this->particles_->getBaseMaterial()))
{
    particles_->registerVariable(mu_srd_, "SRDViscosity");
    particles_->addVariableToWrite<Real>("SRDViscosity");
    particles_->registerVariable(scalar_shear_rate_, "ScalarShearRate");
    particles_->addVariableToWrite<Real>("ScalarShearRate");
}

void ShearRateDependentViscosity::interaction(size_t index_i, Real dt)
{
    //* Viscosity Calculation *//
    Real tau_y = herschel_bulkley_fluid_.getYieldStress();
    Real m = herschel_bulkley_fluid_.getPowerIndex();
    Real K = herschel_bulkley_fluid_.getConsistencyIndex();
    Real min_shear_rate = herschel_bulkley_fluid_.getMinShearRate();
    Real max_shear_rate = herschel_bulkley_fluid_.getMaxShearRate();

    //? Rate of Strain Tensor
    Matd D_ij = 0.5 * (combined_velocity_gradient_[index_i] + combined_velocity_gradient_[index_i].transpose());

    //? D_ij D_ij
    // Real second_invariant = D_ij(0, 0) * D_ij(0, 0) + D_ij(1, 0) * D_ij(1, 0) + D_ij(2, 0) * D_ij(2, 0) + D_ij(0, 1) * D_ij(0, 1) + D_ij(1, 1) * D_ij(1, 1) + D_ij(2, 1) * D_ij(2, 1) + D_ij(0, 2) * D_ij(0, 2) + D_ij(1, 2) * D_ij(1, 2) + D_ij(2, 2) * D_ij(2, 2);
    D_ij = D_ij.cwiseProduct(D_ij);
    Real second_invariant = D_ij.sum();

    //? sqrt(2*D_ij D_ij)
    second_invariant = (Real)std::sqrt(2 * second_invariant);

    //? Herschel Blukley Model tau_y / gamma + K*gamma^(n-1) | Excess of shear rate prevention
    Real capped_shear_rate = std::max(second_invariant, min_shear_rate);
    capped_shear_rate = std::min(capped_shear_rate, max_shear_rate);
    scalar_shear_rate_[index_i] = capped_shear_rate;
    mu_srd_[index_i] = (tau_y + K * (std::pow(capped_shear_rate, m))) / capped_shear_rate;
}
} // namespace fluid_dynamics
} // namespace SPH
