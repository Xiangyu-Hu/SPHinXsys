/**
 * @file 	eulerian_weakly_compressible_fluid_dynamics_inner.cpp
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
 */

#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.hpp"

 //=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		//=================================================================================================//
		EulerianFlowTimeStepInitialization::EulerianFlowTimeStepInitialization(SPHBody &sph_body)
			: ParticleDynamicsSimple(sph_body), EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			rho_n_(particles_->rho_n_), pos_n_(particles_->pos_n_), mass_(particles_->mass_),
			vel_n_(particles_->vel_n_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			gravity_(gravity_ptr_keeper_.createPtr<Gravity>(Vecd(0))) {}
		//=================================================================================================//
		EulerianFlowTimeStepInitialization::
			EulerianFlowTimeStepInitialization(SPHBody &sph_body, Gravity &gravity)
			: ParticleDynamicsSimple(sph_body), EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			rho_n_(particles_->rho_n_), pos_n_(particles_->pos_n_), mass_(particles_->mass_),
			vel_n_(particles_->vel_n_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			gravity_(&gravity) {}
		//=================================================================================================//
		void EulerianFlowTimeStepInitialization::setupDynamics(Real dt)
		{
			particles_->total_ghost_particles_ = 0;
		}
		//=================================================================================================//
		void EulerianFlowTimeStepInitialization::Update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_n_[index_i] * gravity_->InducedAcceleration(pos_n_[index_i]);
		}
		//=================================================================================================//
		FreeSurfaceIndicationInner::
			FreeSurfaceIndicationInner(BaseBodyRelationInner &inner_relation, Real thereshold)
			: InteractionDynamicsWithUpdate(*inner_relation.sph_body_),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			thereshold_by_dimensions_(thereshold * (Real)Dimensions),
			Vol_(particles_->Vol_),
			surface_indicator_(particles_->surface_indicator_),
			smoothing_length_(inner_relation.sph_body_->sph_adaptation_->ReferenceSmoothingLength())
		{
			particles_->registerAVariable<Real>(pos_div_, "PositionDivergence");
		}
		//=================================================================================================//
		void FreeSurfaceIndicationInner::Interaction(size_t index_i, Real dt)
		{
			Real pos_div = 0.0;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				pos_div -= inner_neighborhood.dW_ij_[n] * inner_neighborhood.r_ij_[n] * Vol_[inner_neighborhood.j_[n]];
			}
			pos_div_[index_i] = pos_div;
		}
		//=================================================================================================//
		void FreeSurfaceIndicationInner::Update(size_t index_i, Real dt)
		{
			bool is_free_surface = pos_div_[index_i] < thereshold_by_dimensions_ ? true : false;

			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				/** Two layer particles.*/
				if (pos_div_[inner_neighborhood.j_[n]] < thereshold_by_dimensions_ &&
					inner_neighborhood.r_ij_[n] < smoothing_length_)
				{
					is_free_surface = true;
					break;
				}
			}
			surface_indicator_[index_i] = is_free_surface ? 1 : 0;
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
			vel_n_(particles_->vel_n_),
			dmom_dt_prior_(particles_->dmom_dt_prior_),
			mu_(material_->ReferenceViscosity()),
			smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			const Vecd &vel_i = vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//viscous force
				vel_derivative = (vel_i - vel_n_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}

			dmom_dt_prior_[index_i] += rho_i * acceleration;
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(EulerianFluidBody &fluid_body)
			: ParticleDynamicsReduce<Real, ReduceMax>(fluid_body),
			EulerianWeaklyCompressibleFluidDataSimple(fluid_body), rho_n_(particles_->rho_n_),
			p_(particles_->p_), vel_n_(particles_->vel_n_),
			smoothing_length_(sph_adaptation_->ReferenceSmoothingLength())
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return material_->getSoundSpeed(p_[index_i], rho_n_[index_i]) + vel_n_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			particles_->signal_speed_max_ = reduced_value;
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		VorticityInner::
			VorticityInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			Vol_(particles_->Vol_), vel_n_(particles_->vel_n_)
		{
			particles_->registerAVariable<AngularVecd>(vorticity_, "VorticityInner");
			particles_->addAVariableToWrite<AngularVecd>("VorticityInner");
		}
		//=================================================================================================//
		void VorticityInner::Interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_i = vel_n_[index_i];

			AngularVecd vorticity(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd r_ij = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

				Vecd vel_diff = vel_i - vel_n_[index_j];
				vorticity += SimTK::cross(vel_diff, r_ij) * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}

			vorticity_[index_i] = vorticity;
		}
		//=================================================================================================//
		BaseRelaxation::BaseRelaxation(BaseBodyRelationInner &inner_relation)
			: ParticleDynamics1Level(*inner_relation.sph_body_),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_), rho_n_(particles_->rho_n_),
			p_(particles_->p_), drho_dt_(particles_->drho_dt_), vel_n_(particles_->vel_n_), mom_(particles_->mom_),
			dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BasePressureRelaxation::Initialization(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = material_->getPressure(rho_n_[index_i]);
		}
		//=================================================================================================//
		void BasePressureRelaxation::Update(size_t index_i, Real dt)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_n_[index_i] = mom_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		BaseDensityAndEnergyRelaxation::
			BaseDensityAndEnergyRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BaseDensityAndEnergyRelaxation::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void NonReflectiveBoundaryVariableCorrection::Interaction(size_t index_i, Real dt)
		{
			ComplexShape &body_shape = sph_body_->body_shape_;
			if (surface_indicator_[index_i] == 1)
			{
				Vecd normal_direction = body_shape.findNormalDirection(pos_n_[index_i]);
				n_[index_i] = normal_direction;
				Real velocity_farfield_normal = dot(vel_farfield_, n_[index_i]);
				Real velocity_boundary_normal = dot(vel_n_[index_i], n_[index_i]);

				//judge it is the inflow condition
				if (n_[index_i][0] <= 0.0 | fabs(n_[index_i][1]) > fabs(n_[index_i][0]))
				{
					//supersonic inflow condition
					if (fabs(velocity_boundary_normal) >= sound_speed_)
					{
						vel_n_[index_i] = vel_farfield_;
						p_[index_i] = p_farfield_;
						rho_n_[index_i] = rho_farfield_;
						mom_[index_i] = rho_n_[index_i] * vel_n_[index_i];
					}

					//subsonic inflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						size_t total_inner_neighbor_particles = 0;
						Neighborhood& inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_n_[index_j];
								vel_normal_summation += dot(vel_n_[index_j], n_[index_i]);
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_n_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
						vel_n_[index_i] = vel_normal * n_[index_i] + (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
						mom_[index_i] = rho_n_[index_i] * vel_n_[index_i];
					}
				}
				//judge it is the outflow condition
				else
				{
					//supersonic outflow condition
					if (fabs(velocity_boundary_normal) >= sound_speed_)
					{
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Vecd vel_summation(0.0);
						size_t total_inner_neighbor_particles = 0;
						Neighborhood& inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								rho_summation += rho_n_[index_j];
								vel_summation += vel_n_[index_j];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_average = vel_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);


						p_[index_i] = p_average + TinyReal;
						rho_n_[index_i] = rho_average + TinyReal;
						vel_n_[index_i] = vel_average + TinyReal;
						mom_[index_i] = rho_n_[index_i] * vel_n_[index_i];
					}

					//subsonic outflow condition
					if (fabs(velocity_boundary_normal) < sound_speed_)
					{
						Real inner_weight_summation = 0.0;
						Real rho_summation = 0.0;
						Real p_summation = 0.0;
						Real vel_normal_summation(0.0);
						Vecd vel_tangential_summation(0.0);
						size_t total_inner_neighbor_particles = 0;
						Neighborhood& inner_neighborhood = inner_configuration_[index_i];
						for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
						{
							size_t index_j = inner_neighborhood.j_[n];
							if (surface_indicator_[index_j] != 1)
							{
								Real W_ij = inner_neighborhood.W_ij_[n];
								inner_weight_summation += W_ij * Vol_[index_j];
								rho_summation += rho_n_[index_j];
								vel_normal_summation += dot(vel_n_[index_j], n_[index_i]);
								vel_tangential_summation += vel_n_[index_j] - dot(vel_n_[index_j], n_[index_i])*n_[index_i];
								p_summation += p_[index_j];
								total_inner_neighbor_particles += 1;
							}
						}
						Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
						Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
						Vecd vel_tangential_average = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);
						Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);


						p_[index_i] = p_average * inner_weight_summation + p_farfield_ * (1.0 - inner_weight_summation);
						rho_n_[index_i] = rho_average * inner_weight_summation + rho_farfield_ * (1.0 - inner_weight_summation);
						Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
						vel_n_[index_i] = vel_normal * n_[index_i] + vel_tangential_average;
						mom_[index_i] = rho_n_[index_i] * vel_n_[index_i];
					}
				}
			}
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//