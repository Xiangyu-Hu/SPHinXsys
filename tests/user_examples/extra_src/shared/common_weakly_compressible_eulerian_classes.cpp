
#include "common_weakly_compressible_eulerian_classes.h"

namespace SPH
{
//=================================================================================================//
EulerianWCTimeStepInitialization::EulerianWCTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : TimeStepInitialization(sph_body, gravity_ptr), rho_(particles_->rho_), pos_(particles_->pos_), vel_(particles_->vel_),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
void EulerianWCTimeStepInitialization::update(size_t index_i, Real dt)
{
    dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
}
//=================================================================================================//
Real EulerianWCAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return acousticCFL_ / Dimensions * smoothing_length_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
WCEulerianViscousAccelerationInner::WCEulerianViscousAccelerationInner(BaseInnerRelation &inner_relation)
    : BaseViscousAccelerationInner(inner_relation),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
void WCEulerianViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Vecd vel_i = vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + TinyReal);
        acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
    }
    dmom_dt_prior_[index_i] += rho_i * acceleration;
}
//=================================================================================================//
FluidStarState NoRiemannSolverInWCEulerianMethod::getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;

    Real p_star = (rhol_cl * state_j.p_ + rhor_cr * state_i.p_) / (rhol_cl + rhor_cr);
    Real u_star = (rhol_cl * ul + rhor_cr * ur) / (rhol_cl + rhor_cr);
    Vecd vel_star = (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) - e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));

    FluidStarState interface_state(vel_star, p_star);
    interface_state.vel_ = vel_star;
    interface_state.p_ = p_star;

    return interface_state;
}
//=================================================================================================//
FluidStarState AcousticRiemannSolverInEulerianMethod::getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
    Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);

    Real p_star = (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ + rhol_cl * rhor_cr * (ul - ur) * SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1))) / (rhol_cl + rhor_cr);
    Real u_star = (rhol_cl * ul + rhor_cr * ur + (state_i.p_ - state_j.p_) * pow(SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1)), 2)) / (rhol_cl + rhor_cr);
    Vecd vel_star = (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) - e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));

    FluidStarState interface_state(vel_star, p_star);
    interface_state.vel_ = vel_star;
    interface_state.p_ = p_star;

    return interface_state;
}
//=================================================================================================//
EulerianBaseIntegration::EulerianBaseIntegration(BaseInnerRelation &inner_relation) : BaseIntegration(inner_relation),
                                                                                      Vol_(particles_->Vol_), mom_(*particles_->getVariableByName<Vecd>("Momentum")),
                                                                                      dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
                                                                                      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
NonReflectiveBoundaryVariableCorrection::NonReflectiveBoundaryVariableCorrection(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner<BaseParticles>(inner_relation),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), rho_farfield_(0.0), sound_speed_(0.0), vel_farfield_(Vecd::Zero()),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), Vol_(particles_->Vol_), vel_(particles_->vel_),
      mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
      surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator"))
{
    particles_->registerVariable(n_, "NormalDirection");
    particles_->registerVariable(inner_weight_summation_, "InnerWeightSummation");
    particles_->registerVariable(rho_average_, "DensityAverage");
    particles_->registerVariable(vel_normal_average_, "VelocityNormalAverage");
    particles_->registerVariable(vel_tangential_average_, "VelocityTangentialAverage");
    particles_->registerVariable(vel_average_, "VelocityAverage");
    particles_->registerVariable(surface_inner_particle_indicator_, "SurfaceInnerParticleIndicator");
};
//=================================================================================================//
void NonReflectiveBoundaryVariableCorrection::initialization(size_t index_i, Real dt)
{
    if (surface_indicator_[index_i] == 1)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            // if particle j could be searched by indicator 1, this partical j also is seemed as boundary condition particle
            size_t index_j = inner_neighborhood.j_[n];
            surface_inner_particle_indicator_[index_j] = 1;
        }
    };
}
//=================================================================================================//
void NonReflectiveBoundaryVariableCorrection::interaction(size_t index_i, Real dt)
{
    Shape &body_shape = *sph_body_.body_shape_;
    if (surface_indicator_[index_i] == 1 || surface_inner_particle_indicator_[index_i] == 1)
    {
        Vecd normal_direction = body_shape.findNormalDirection(pos_[index_i]);
        n_[index_i] = normal_direction;
        Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);
        // judge it is the inflow condition
        if (n_[index_i][0] <= 0.0 || fabs(n_[index_i][1]) > fabs(n_[index_i][0]))
        {
            // subsonic inflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                inner_weight_summation_[index_i] = 0.0;
                Real rho_summation = 0.0;
                Real vel_normal_summation(0.0);
                size_t total_inner_neighbor_particles = 0;
                const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (surface_indicator_[index_j] != 1)
                    {
                        Real W_ij = inner_neighborhood.W_ij_[n];
                        inner_weight_summation_[index_i] += W_ij * Vol_[index_j];
                        rho_summation += rho_[index_j];
                        vel_normal_summation += vel_[index_j].dot(n_[index_i]);
                        total_inner_neighbor_particles += 1;
                    }
                }
                rho_average_[index_i] = rho_summation / (total_inner_neighbor_particles + TinyReal);
                vel_normal_average_[index_i] = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
            }
        }
        // judge it is the outflow condition
        else
        {
            // supersonic outflow condition
            if (fabs(velocity_boundary_normal) >= sound_speed_)
            {
                Real rho_summation = 0.0;
                Vecd vel_summation = Vecd::Zero();
                size_t total_inner_neighbor_particles = 0;
                const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (surface_indicator_[index_j] != 1)
                    {
                        rho_summation += rho_[index_j];
                        vel_summation += vel_[index_j];
                        total_inner_neighbor_particles += 1;
                    }
                }
                rho_average_[index_i] = rho_summation / (total_inner_neighbor_particles + TinyReal);
                vel_average_[index_i] = vel_summation / (total_inner_neighbor_particles + TinyReal);
            }

            // subsonic outflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                inner_weight_summation_[index_i] = 0.0;
                Real rho_summation = 0.0;
                Real vel_normal_summation(0.0);
                Vecd vel_tangential_summation = Vecd::Zero();
                size_t total_inner_neighbor_particles = 0;
                const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (surface_indicator_[index_j] != 1)
                    {
                        Real W_ij = inner_neighborhood.W_ij_[n];
                        inner_weight_summation_[index_i] += W_ij * Vol_[index_j];
                        rho_summation += rho_[index_j];
                        vel_normal_summation += vel_[index_j].dot(n_[index_i]);
                        vel_tangential_summation += vel_[index_j] - vel_[index_j].dot(n_[index_i]) * n_[index_i];
                        total_inner_neighbor_particles += 1;
                    }
                }
                rho_average_[index_i] = rho_summation / (total_inner_neighbor_particles + TinyReal);
                vel_normal_average_[index_i] = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
                vel_tangential_average_[index_i] = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);
            }
        }
    }
}
//=================================================================================================//
void NonReflectiveBoundaryVariableCorrection::update(size_t index_i, Real dt)
{
    Shape &body_shape = *sph_body_.body_shape_;
    if (surface_indicator_[index_i] == 1 || surface_inner_particle_indicator_[index_i] == 1)
    {
        Vecd normal_direction = body_shape.findNormalDirection(pos_[index_i]);
        n_[index_i] = normal_direction;
        Real velocity_farfield_normal = vel_farfield_.dot(n_[index_i]);
        Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);

        // judge it is the inflow condition
        if (n_[index_i][0] <= 0.0 || fabs(n_[index_i][1]) > fabs(n_[index_i][0]))
        {
            // supersonic inflow condition
            if (fabs(velocity_boundary_normal) >= sound_speed_)
            {
                vel_[index_i] = vel_farfield_;
                rho_[index_i] = rho_farfield_;
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
            // subsonic inflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
        }
        // judge it is the outflow condition
        else
        {
            // supersonic outflow condition
            if (fabs(velocity_boundary_normal) >= sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] + TinyReal;
                vel_[index_i] = vel_average_[index_i];
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }

            // subsonic outflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average_[index_i];
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
        }
    }
    //=================================================================================================//
}
} // namespace SPH
  //=================================================================================================//