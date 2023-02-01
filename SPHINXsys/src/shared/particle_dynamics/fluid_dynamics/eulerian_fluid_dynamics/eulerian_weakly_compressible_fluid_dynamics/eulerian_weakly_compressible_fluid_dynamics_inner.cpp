#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.hpp"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{	//=================================================================================================//
		WeaklyCompressibleFluidInitialCondition::
			WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body), EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			  pos_(particles_->pos_), vel_(particles_->vel_), mom_(particles_->mom_),
			  rho_(particles_->rho_), p_(particles_->p_) {}
		//=================================================================================================//
		EulerianFlowTimeStepInitialization::
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
			: BaseTimeStepInitialization(sph_body, gravity_ptr),
			  EulerianWeaklyCompressibleFluidDataSimple(sph_body), rho_(particles_->rho_),
			  pos_(particles_->pos_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		void EulerianFlowTimeStepInitialization::update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_),
			  EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			  vel_(particles_->vel_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			  mu_(particles_->fluid_.ReferenceViscosity()),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd &vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				// viscous force
				vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}

			dmom_dt_prior_[index_i] += rho_i * acceleration;
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			  EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			  fluid_(particles_->fluid_), rho_(particles_->rho_),
			  p_(particles_->p_), vel_(particles_->vel_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		BaseIntegration::BaseIntegration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_),
			  EulerianWeaklyCompressibleFluidDataInner(inner_relation), fluid_(particles_->fluid_),
			  Vol_(particles_->Vol_), mass_(particles_->mass_), rho_(particles_->rho_),
			  p_(particles_->p_), drho_dt_(particles_->drho_dt_), vel_(particles_->vel_), mom_(particles_->mom_),
			  dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		void NonReflectiveBoundaryVariableCorrection::interaction(size_t index_i, Real dt)
		{
			Shape &body_shape = *sph_body_.body_shape_;
			if (surface_indicator_[index_i] == 1)
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
						p_[index_i] = p_farfield_;
						rho_[index_i] = rho_farfield_;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}

					// subsonic inflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_[index_j];
								vel_normal_summation += vel_[index_j].dot(n_[index_i]);
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
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
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Vecd vel_summation = Vecd::Zero();
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								rho_summation += rho_[index_j];
								vel_summation += vel_[index_j];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_average = vel_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average;
						rho_[index_i] = rho_average;
						vel_[index_i] = vel_average;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}

					// subsonic outflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						Vecd vel_tangential_summation = Vecd::Zero();
						size_t total_inner_neighbor_particles = 0;
						Neighborhood &inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_[index_j];
								vel_normal_summation += vel_[index_j].dot(n_[index_i]);
								vel_tangential_summation += vel_[index_j] - (vel_[index_j].dot(n_[index_i])) * n_[index_i];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_tangential_average = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
						vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average;
						mom_[index_i] = rho_[index_i] * vel_[index_i];
					}
				}
			}
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//