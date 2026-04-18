#include "eulerian_compressible_fluid_integration.h"
#include "muscl_reconstruction.hpp"

#include "adaptation.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BaseIntegrationInCompressible::BaseIntegrationInCompressible(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation),
      compressible_fluid_(CompressibleFluid(1.0, 1.4)),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      E_(particles_->registerStateVariableData<Real>("TotalEnergy")),
      dE_dt_(particles_->registerStateVariableData<Real>("TotalEnergyChangeRate")),
      dmass_dt_(particles_->registerStateVariableData<Real>("MassChangeRate")),
      mom_(particles_->registerStateVariableData<Vecd>("Momentum")),
      force_(particles_->registerStateVariableData<Vecd>("Force")),
      force_prior_(particles_->registerStateVariableData<Vecd>("ForcePrior")) {};
//=================================================================================================//
inline MUSCLHLLCBridgeConfig make_default_bridge_config(CompressibleFluid &fluid, MUSCLHLLCBridgeConfig cfg_in)
{
    cfg_in.muscl_cfg.gamma = fluid.HeatCapacityRatio();
    cfg_in.muscl_cfg.limiter = SlopeLimiter::MC;
    return cfg_in;
}
//-------------------------------------------------------------------------------------------------//
inline Vecd row_to_vecd(const Matd &m, size_t r)
{
#if SPH_NDIM == 2
    return Vecd(m(r, 0), m(r, 1));
#else
    return Vecd(m(r, 0), m(r, 1), m(r, 2));
#endif
}
//=================================================================================================//
EulerianCompressibleIntegration1stHalfMUSCL<Inner<>>::EulerianCompressibleIntegration1stHalfMUSCL(
    BaseInnerRelation &inner_relation, const MUSCLHLLCBridgeConfig &bridge_cfg)
    : BaseIntegrationInCompressible(inner_relation),
      bridge_cfg_(make_default_bridge_config(compressible_fluid_, bridge_cfg)),
      bridge_(compressible_fluid_, compressible_fluid_, bridge_cfg_),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")) {}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration1stHalfMUSCL<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    Real &rho_i = rho_[index_i];
    Vecd &vel_i = vel_[index_i];
    Real &p_i = p_[index_i];
    CompressibleFluidState state_i(rho_i, vel_i, p_i, energy_per_volume_i);
    Vecd momentum_change_rate = force_prior_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        Real &rho_j = rho_[index_j];
        Vecd &vel_j = vel_[index_j];
        Real &p_j = p_[index_j];
        CompressibleFluidState state_j(rho_j, vel_j, p_j, energy_per_volume_j);

        Vecd grad_rho_i = rho_grad_[index_i];
        Vecd grad_rho_j = rho_grad_[index_j];
        Matd vg_i = vel_grad_[index_i];
        Matd vg_j = vel_grad_[index_j];
        Vecd grad_u_i = row_to_vecd(vg_i, 0);
        Vecd grad_u_j = row_to_vecd(vg_j, 0);
        Vecd grad_v_i = row_to_vecd(vg_i, 1);
        Vecd grad_v_j = row_to_vecd(vg_j, 1);
        Vecd grad_p_i = p_grad_[index_i];
        Vecd grad_p_j = p_grad_[index_j];

        const Vecd &xi = pos_[index_i];
        const Vecd &xj = pos_[index_j];
        Vecd xf = 0.5 * (xi + xj);

        CompressibleFluidStarState interface_state = bridge_.getInterfaceState(
            state_i, state_j, xi, xj, xf, e_ij,
            grad_rho_i, grad_rho_j,
            grad_u_i, grad_u_j,
            grad_v_i, grad_v_j,
#if SPH_NDIM == 3
            row_to_vecd(vg_i, 2), row_to_vecd(vg_j, 2),
#endif
            grad_p_i, grad_p_j);

        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
        momentum_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j *
                                (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij;
    }
    force_[index_i] = momentum_change_rate;
}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration1stHalfMUSCL<Inner<>>::update(size_t index_i, Real dt)
{
    mom_[index_i] += force_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//-------------------------------------------------------------------------------------------------//
EulerianCompressibleIntegration2ndHalfMUSCL<Inner<>>::EulerianCompressibleIntegration2ndHalfMUSCL(
    BaseInnerRelation &inner_relation, const MUSCLHLLCBridgeConfig &bridge_cfg)
    : BaseIntegrationInCompressible(inner_relation),
      bridge_cfg_(make_default_bridge_config(compressible_fluid_, bridge_cfg)),
      bridge_(compressible_fluid_, compressible_fluid_, bridge_cfg_),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")) {}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration2ndHalfMUSCL<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    Real &rho_i = rho_[index_i];
    Vecd &vel_i = vel_[index_i];
    Real &p_i = p_[index_i];
    CompressibleFluidState state_i(rho_i, vel_i, p_i, energy_per_volume_i);
    Real mass_change_rate = 0.0;
    Real energy_change_rate = force_prior_[index_i].dot(vel_[index_i]);
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        Real &rho_j = rho_[index_j];
        Vecd &vel_j = vel_[index_j];
        Real &p_j = p_[index_j];
        CompressibleFluidState state_j(rho_j, vel_j, p_j, energy_per_volume_j);

        Vecd grad_rho_i = rho_grad_[index_i];
        Vecd grad_rho_j = rho_grad_[index_j];
        Matd vg_i = vel_grad_[index_i];
        Matd vg_j = vel_grad_[index_j];
        Vecd grad_u_i = row_to_vecd(vg_i, 0);
        Vecd grad_u_j = row_to_vecd(vg_j, 0);
        Vecd grad_v_i = row_to_vecd(vg_i, 1);
        Vecd grad_v_j = row_to_vecd(vg_j, 1);
        Vecd grad_p_i = p_grad_[index_i];
        Vecd grad_p_j = p_grad_[index_j];

        const Vecd &xi = pos_[index_i];
        const Vecd &xj = pos_[index_j];
        Vecd xf = 0.5 * (xi + xj);

        CompressibleFluidStarState interface_state = bridge_.getInterfaceState(
            state_i, state_j, xi, xj, xf, e_ij,
            grad_rho_i, grad_rho_j,
            grad_u_i, grad_u_j,
            grad_v_i, grad_v_j,
#if SPH_NDIM == 3
            row_to_vecd(vg_i, 2), row_to_vecd(vg_j, 2),
#endif
            grad_p_i, grad_p_j);

        mass_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j *
                            (interface_state.rho_ * interface_state.vel_).dot(e_ij);
        energy_change_rate -= 2.0 * Vol_[index_i] * dW_ijV_j *
                              ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(e_ij);
    }
    dmass_dt_[index_i] = mass_change_rate;
    dE_dt_[index_i] = energy_change_rate;
}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration2ndHalfMUSCL<Inner<>>::update(size_t index_i, Real dt)
{
    E_[index_i] += dE_dt_[index_i] * dt;
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];
    Real rho_e = E_[index_i] / Vol_[index_i] - 0.5 * (mom_[index_i] / mass_[index_i]).squaredNorm() * rho_[index_i];
    p_[index_i] = compressible_fluid_.getPressure(rho_[index_i], rho_e);
}
//-------------------------------------------------------------------------------------------------//
using MUSCLWallBase = InteractionWithWall<BaseIntegrationInCompressibleForWall>;
//-------------------------------------------------------------------------------------------------//
EulerianCompressibleIntegration1stHalfMUSCL<Contact<Wall>>::EulerianCompressibleIntegration1stHalfMUSCL(
    BaseContactRelation &contact_relation, const MUSCLHLLCBridgeConfig &bridge_cfg)
    : MUSCLWallBase(contact_relation),
      bridge_cfg_(make_default_bridge_config(this->compressible_fluid_, bridge_cfg)),
      bridge_(this->compressible_fluid_, this->compressible_fluid_, bridge_cfg_),
      rho_grad_(this->particles_->getVariableDataByName<Vecd>("DensityGradient")),
      vel_grad_(this->particles_->getVariableDataByName<Matd>("VelocityGradient")),
      p_grad_(this->particles_->getVariableDataByName<Vecd>("PressureGradient")) {}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration1stHalfMUSCL<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    CompressibleFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i], energy_per_volume_i);
    Vecd momentum_change_rate = this->force_prior_[index_i];
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            Real rho_g = this->rho_[index_i];
            Real p_g = this->p_[index_i];
            Real energy_per_volume_g = energy_per_volume_i;
            CompressibleFluidState state_g(rho_g, vel_reflect, p_g, energy_per_volume_g);

            Vecd grad_rho_i = rho_grad_[index_i];
            Vecd grad_rho_g = grad_rho_i;
            Matd vg_i = vel_grad_[index_i];
            Vecd grad_u_i = row_to_vecd(vg_i, 0);
            Vecd grad_v_i = row_to_vecd(vg_i, 1);
            Vecd grad_u_g = grad_u_i;
            Vecd grad_v_g = grad_v_i;
            Vecd grad_p_i = p_grad_[index_i];
            Vecd grad_p_g = grad_p_i;

            const Vecd &xi = this->pos_[index_i];
            Vecd xj = xi - contact_neighborhood.r_ij_[n] * e_ij;
            Vecd xf = 0.5 * (xi + xj);

            CompressibleFluidStarState interface_state = bridge_.getInterfaceState(
                state_i, state_g, xi, xj, xf, e_ij,
                grad_rho_i, grad_rho_g,
                grad_u_i, grad_u_g,
                grad_v_i, grad_v_g,
#if SPH_NDIM == 3
                row_to_vecd(vg_i, 2), row_to_vecd(vg_i, 2),
#endif
                grad_p_i, grad_p_g);

            Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * this->Vol_[index_i] * dW_ijV_j *
                                    (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij;
        }
    }
    this->force_[index_i] += momentum_change_rate;
}
//-------------------------------------------------------------------------------------------------//
EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>>::EulerianCompressibleIntegration2ndHalfMUSCL(
    BaseContactRelation &contact_relation, const MUSCLHLLCBridgeConfig &bridge_cfg)
    : MUSCLWallBase(contact_relation),
      bridge_cfg_(make_default_bridge_config(this->compressible_fluid_, bridge_cfg)),
      bridge_(this->compressible_fluid_, this->compressible_fluid_, bridge_cfg_),
      rho_grad_(this->particles_->getVariableDataByName<Vecd>("DensityGradient")),
      vel_grad_(this->particles_->getVariableDataByName<Matd>("VelocityGradient")),
      p_grad_(this->particles_->getVariableDataByName<Vecd>("PressureGradient")) {}
//-------------------------------------------------------------------------------------------------//
void EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    CompressibleFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i], energy_per_volume_i);
    Real mass_change_rate = 0.0;
    Real energy_change_rate = this->force_prior_[index_i].dot(this->vel_[index_i]);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            Real rho_g = this->rho_[index_i];
            Real p_g = this->p_[index_i];
            Real energy_per_volume_g = energy_per_volume_i;
            CompressibleFluidState state_g(rho_g, vel_reflect, p_g, energy_per_volume_g);

            Vecd grad_rho_i = rho_grad_[index_i];
            Vecd grad_rho_g = grad_rho_i;
            Matd vg_i = vel_grad_[index_i];
            Vecd grad_u_i = row_to_vecd(vg_i, 0);
            Vecd grad_v_i = row_to_vecd(vg_i, 1);
            Vecd grad_u_g = grad_u_i;
            Vecd grad_v_g = grad_v_i;
            Vecd grad_p_i = p_grad_[index_i];
            Vecd grad_p_g = grad_p_i;

            const Vecd &xi = this->pos_[index_i];
            Vecd xj = xi - contact_neighborhood.r_ij_[n] * e_ij;
            Vecd xf = 0.5 * (xi + xj);

            CompressibleFluidStarState interface_state = bridge_.getInterfaceState(
                state_i, state_g, xi, xj, xf, e_ij,
                grad_rho_i, grad_rho_g,
                grad_u_i, grad_u_g,
                grad_v_i, grad_v_g,
#if SPH_NDIM == 3
                row_to_vecd(vg_i, 2), row_to_vecd(vg_i, 2),
#endif
                grad_p_i, grad_p_g);

            mass_change_rate -= 2.0 * this->Vol_[index_i] * dW_ijV_j *
                                (interface_state.rho_ * interface_state.vel_).dot(e_ij);
            energy_change_rate -= 2.0 * this->Vol_[index_i] * dW_ijV_j *
                                  ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(e_ij);
        }
    }
    this->dmass_dt_[index_i] += mass_change_rate;
    this->dE_dt_[index_i] += energy_change_rate;
}
//=================================================================================================//
CompressibleFluidInitialCondition::CompressibleFluidInitialCondition(SPHBody &sph_body)
    : FluidInitialCondition(sph_body),
      mom_(particles_->getVariableDataByName<Vecd>("Momentum")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      E_(particles_->getVariableDataByName<Real>("TotalEnergy")) {}
//=================================================================================================//
EulerianCompressibleAcousticTimeStepSize::
    EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
    : AcousticTimeStep(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()),
      compressible_fluid_(CompressibleFluid(1.0, 1.4))
{
    acousticCFL_ = acousticCFL;
};
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return acousticCFL_ / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH