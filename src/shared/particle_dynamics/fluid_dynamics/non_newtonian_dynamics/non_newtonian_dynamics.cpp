#include "non_newtonian_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
Oldroyd_BIntegration1stHalf ::Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalfInnerDissipativeRiemann(inner_relation)
{
    particles_->registerVariable(tau_, "ElasticStress");
    particles_->registerVariable(dtau_dt_, "ElasticStressChangeRate");
    particles_->registerSortableVariable<Matd>("ElasticStress");
    particles_->addVariableToRestart<Matd>("ElasticStress");
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    Integration1stHalfInnerDissipativeRiemann::initialization(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf::interaction(size_t index_i, Real dt)
{
    Integration1stHalfInnerDissipativeRiemann::interaction(index_i, dt);

    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        // elastic force
        acceleration += (tau_[index_i] + tau_[index_j]) * nablaW_ijV_j;
    }

    acc_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalf::
    Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : Integration2ndHalfInnerDissipativeRiemann(inner_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf::update(size_t index_i, Real dt)
{
    Integration2ndHalfInnerDissipativeRiemann::update(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf::
    interaction(size_t index_i, Real dt)
{
    Integration2ndHalfInnerDissipativeRiemann::interaction(index_i, dt);

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
Oldroyd_BMomentumWallBoundary::
    Oldroyd_BMomentumWallBoundary(BaseContactRelation &wall_contact_relation)
    : MomentumWallBoundaryDissipativeRiemann(wall_contact_relation),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")){};
//=================================================================================================//
void Oldroyd_BMomentumWallBoundary::interaction(size_t index_i, Real dt)
{
    MomentumWallBoundaryDissipativeRiemann::interaction(index_i, dt);

    Real rho_i = rho_[index_i];
    Matd tau_i = tau_[index_i];

    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            /** stress boundary condition. */
            acceleration += 2.0 * tau_i * nablaW_ijV_j / rho_i;
        }
    }

    acc_[index_i] += acceleration;
}
//=================================================================================================//
Oldroyd_BContinuityWallBoundary::
    Oldroyd_BContinuityWallBoundary(BaseContactRelation &wall_contact_relation)
    : ContinuityWallBoundaryDissipativeRiemann(wall_contact_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BContinuityWallBoundary::interaction(size_t index_i, Real dt)
{
    ContinuityWallBoundaryDissipativeRiemann::interaction(index_i, dt);

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
} // namespace fluid_dynamics
} // namespace SPH
