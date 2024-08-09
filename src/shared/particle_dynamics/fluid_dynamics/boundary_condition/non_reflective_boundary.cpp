#include "non_reflective_boundary.h"

#include "weakly_compressible_fluid.h"
namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
NonReflectiveBoundaryCorrection::NonReflectiveBoundaryCorrection(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      rho_farfield_(0.0), sound_speed_(0.0), vel_farfield_(Vecd::Zero()),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      mom_(particles_->getVariableDataByName<Vecd>("Momentum")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      inner_weight_summation_(particles_->registerStateVariable<Real>("InnerWeightSummation")),
      rho_average_(particles_->registerStateVariable<Real>("DensityAverage")),
      vel_normal_average_(particles_->registerStateVariable<Real>("VelocityNormalAverage")),
      vel_tangential_average_(particles_->registerStateVariable<Vecd>("VelocityTangentialAverage")),
      vel_average_(particles_->registerStateVariable<Vecd>("VelocityAverage")),
      indicator_(particles_->getVariableDataByName<int>("Indicator")),
      smeared_surface_(particles_->getVariableDataByName<int>("SmearedSurface")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")){};
//=================================================================================================//
void NonReflectiveBoundaryCorrection::interaction(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1 || smeared_surface_[index_i] == 1)
    {
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
                    if (indicator_[index_j] != 1)
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
                    if (indicator_[index_j] != 1)
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
                    if (indicator_[index_j] != 1)
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
void NonReflectiveBoundaryCorrection::update(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1 || smeared_surface_[index_i] == 1)
    {
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
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }
            // subsonic inflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
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
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }

            // subsonic outflow condition
            if (fabs(velocity_boundary_normal) < sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average_[index_i];
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }
        }
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
