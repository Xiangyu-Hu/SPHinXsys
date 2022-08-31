/**
 * @file 	fluid_boundary.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_boundary.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseFlowBoundaryCondition::BaseFlowBoundaryCondition(BodyPartByCell &body_part)
			: LocalDynamics(body_part.getSPHBody()), FluidDataSimple(sph_body_),
			  rho_(particles_->rho_), p_(particles_->p_),
			  pos_(particles_->pos_), vel_(particles_->vel_){};
		//=================================================================================================//
		FlowVelocityBuffer::FlowVelocityBuffer(BodyPartByCell &body_part)
			: BaseFlowBoundaryCondition(body_part), relaxation_rate_(0.3){};
		//=================================================================================================//
		void FlowVelocityBuffer::update(size_t index_i, Real dt)
		{
			vel_[index_i] += relaxation_rate_ * (getTargetVelocity(pos_[index_i], vel_[index_i]) - vel_[index_i]);
		}
		//=================================================================================================//
		InflowVelocityCondition::InflowVelocityCondition(BodyAlignedBoxByCell &aligned_box_part)
			: BaseFlowBoundaryCondition(aligned_box_part),
			  transform_(aligned_box_part.aligned_box_.getTransform()),
			  halfsize_(aligned_box_part.aligned_box_.HalfSize()) {}
		//=================================================================================================//
		void InflowVelocityCondition::update(size_t index_i, Real dt)
		{
			Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
			Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
			Vecd prescribed_velocity =
				transform_.xformFrameVecToBase(getPrescribedVelocity(frame_position, frame_velocity));
			vel_[index_i] = prescribed_velocity;
		}
		//=================================================================================================//
		DampingBoundaryCondition::DampingBoundaryCondition(BodyRegionByCell &body_part)
			: BaseFlowBoundaryCondition(body_part), strength_(5.0),
			  damping_zone_bounds_(body_part.body_part_shape_.getBounds()){};
		//=================================================================================================//
		void DampingBoundaryCondition::update(size_t index_i, Real dt)
		{
			Real damping_factor = (pos_[index_i][0] - damping_zone_bounds_.first[0]) /
								  (damping_zone_bounds_.second[0] - damping_zone_bounds_.first[0]);
			vel_[index_i] *= (1.0 - dt * strength_ * damping_factor * damping_factor);
		}
		//=================================================================================================//
		EmitterInflowCondition::
			EmitterInflowCondition(BodyAlignedBoxByParticle &aligned_box_part)
			: LocalDynamics(aligned_box_part.getSPHBody()), FluidDataSimple(sph_body_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  rho_(particles_->rho_), p_(particles_->p_), inflow_pressure_(0),
			  rho0_(material_->ReferenceDensity()), aligned_box_(aligned_box_part.aligned_box_),
			  updated_transform_(aligned_box_.getTransform()),
			  old_transform_(updated_transform_) {}
		//=================================================================================================//
		void EmitterInflowCondition ::update(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			Vecd frame_position = old_transform_.shiftBaseStationToFrame(pos_[sorted_index_i]);
			Vecd frame_velocity = old_transform_.xformBaseVecToFrame(vel_[sorted_index_i]);
			pos_[sorted_index_i] = updated_transform_.shiftFrameStationToBase(frame_position);
			vel_[sorted_index_i] = updated_transform_.xformFrameVecToBase(getTargetVelocity(frame_position, frame_velocity));
			rho_[sorted_index_i] = rho0_;
			p_[sorted_index_i] = material_->getPressure(rho_[sorted_index_i]);
		}
		//=================================================================================================//
		EmitterInflowInjecting::EmitterInflowInjecting(BodyAlignedBoxByParticle &aligned_box_part,
													   size_t body_buffer_width, int axis, bool positive)
			: LocalDynamics(aligned_box_part.getSPHBody()), FluidDataSimple(sph_body_),
			  pos_(particles_->pos_), rho_(particles_->rho_), p_(particles_->p_),
			  axis_(axis), aligned_box_(aligned_box_part.aligned_box_)
		{
			size_t total_body_buffer_particles = aligned_box_part.body_part_particles_.size() * body_buffer_width_;
			particles_->addBufferParticles(total_body_buffer_particles);
			sph_body_.allocateConfigurationMemoriesForBufferParticles();

			checking_bound_ = positive ? std::bind(&EmitterInflowInjecting::checkUpperBound, this, _1, _2)
									   : std::bind(&EmitterInflowInjecting::checkLowerBound, this, _1, _2);
		}
		//=================================================================================================//
		void EmitterInflowInjecting::update(size_t unsorted_index_i, Real dt)
		{
			mutex_switch_to_buffer_.lock();
			checking_bound_(unsorted_index_i, dt);
			mutex_switch_to_buffer_.unlock();
		}
		//=================================================================================================//
		void EmitterInflowInjecting::checkUpperBound(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			if (aligned_box_.checkUpperBound(axis_, pos_[sorted_index_i]))
			{
				if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
				{
					std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
							  << "Not enough body buffer particles! Exit the code."
							  << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
				/** Realize the buffer particle by increasing the number of real particle in the body.  */
				particles_->total_real_particles_ += 1;
				/** Periodic bounding. */
				pos_[sorted_index_i] = aligned_box_.getUpperPeriodic(axis_, pos_[sorted_index_i]);
				rho_[sorted_index_i] = material_->ReferenceDensity();
				p_[sorted_index_i] = material_->getPressure(rho_[sorted_index_i]);
			}
		}
		//=================================================================================================//
		void EmitterInflowInjecting::checkLowerBound(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			if (aligned_box_.checkLowerBound(axis_, pos_[sorted_index_i]))
			{
				if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
				{
					std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
							  << "Not enough body buffer particles! Exit the code."
							  << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
				/** Realize the buffer particle by increasing the number of real particle in the body.  */
				particles_->total_real_particles_ += 1;
				pos_[sorted_index_i] = aligned_box_.getUpperPeriodic(axis_, pos_[sorted_index_i]);
			}
		}
		//=================================================================================================//
		StaticConfinementDensity::StaticConfinementDensity(NearShapeSurface &near_surface)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			  rho0_(particles_->rho0_), inv_sigma0_(1.0 / particles_->sigma0_),
			  mass_(particles_->mass_), rho_sum_(particles_->rho_sum_), pos_(particles_->pos_),
			  level_set_shape_(&near_surface.level_set_shape_) {}
		//=================================================================================================//
		void StaticConfinementDensity::update(size_t index_i, Real dt)
		{
			Real inv_Vol_0_i = rho0_ / mass_[index_i];
			rho_sum_[index_i] +=
				level_set_shape_->computeKernelIntegral(pos_[index_i]) * inv_Vol_0_i * rho0_ * inv_sigma0_;
		}
		//=================================================================================================//
		StaticConfinementPressureRelaxation::StaticConfinementPressureRelaxation(NearShapeSurface &near_surface)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			  rho_(particles_->rho_), p_(particles_->p_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  acc_(particles_->acc_),
			  level_set_shape_(&near_surface.level_set_shape_),
			  riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		void StaticConfinementPressureRelaxation::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			Vecd normal_to_fluid = -kernel_gradient / (kernel_gradient.norm() + TinyReal);

			FluidState state(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd vel_in_wall = -state.vel_;
			FluidState state_in_wall(rho_[index_i], vel_in_wall, p_[index_i]);

			// always solving one-side Riemann problem for wall boundaries
			Real p_star = riemann_solver_.getPStar(state, state_in_wall, normal_to_fluid);
			acc_[index_i] -= 2.0 * p_star * kernel_gradient / state.rho_;
		}
		//=================================================================================================//
		StaticConfinementDensityRelaxation::StaticConfinementDensityRelaxation(NearShapeSurface &near_surface)
			: LocalDynamics(near_surface.getSPHBody()), FluidDataSimple(sph_body_),
			  rho_(particles_->rho_), p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  level_set_shape_(&near_surface.level_set_shape_),
			  riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		void StaticConfinementDensityRelaxation::update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
			Vecd normal_to_fluid = -kernel_gradient / (kernel_gradient.norm() + TinyReal);

			FluidState state(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd vel_in_wall = -state.vel_;
			FluidState state_in_wall(rho_[index_i], vel_in_wall, p_[index_i]);

			// always solving one-side Riemann problem for wall boundaries
			Vecd vel_star = riemann_solver_.getVStar(state, state_in_wall, normal_to_fluid);
			drho_dt_[index_i] += 2.0 * state.rho_ * dot(state.vel_ - vel_star, kernel_gradient);
		}
		//=================================================================================================//
		StaticConfinement::StaticConfinement(NearShapeSurface &near_surface)
			: density_summation_(near_surface), pressure_relaxation_(near_surface),
			  density_relaxation_(near_surface) {}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//