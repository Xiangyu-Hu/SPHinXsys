/**
 * @file 	fluid_dynamics.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_inner.h"
#include "in_output.h"
#include "geometry_level_set.h"
//=================================================================================================//
using namespace std;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FluidInitialCondition::
			FluidInitialCondition(FluidBody* body)
			: ParticleDynamicsSimple(body), FluidDataSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_)
		{
		}
		//=================================================================================================//
		FreeSurfaceIndicationInner::
			FreeSurfaceIndicationInner(BaseInnerBodyRelation* inner_relation, Real thereshold) :
			InteractionDynamicsWithUpdate(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			thereshold_by_dimensions_(thereshold*(Real)Dimensions), 
			Vol_(particles_->Vol_), is_free_surface_(particles_->is_free_surface_)
		{
			particles_->registerAVariable<indexScalar, Real>(pos_div_, "PositionDivergence", true);
			smoothing_length_ = inner_relation->sph_body_->particle_adaptation_->ReferenceSmoothingLength();
		}
		//=================================================================================================//
		void FreeSurfaceIndicationInner::Interaction(size_t index_i, Real dt)
		{
			Real pos_div = 0.0;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				pos_div -= inner_neighborhood.dW_ij_[n] 
					* inner_neighborhood.r_ij_[n] * Vol_[inner_neighborhood.j_[n]];
			}
			pos_div_[index_i] = pos_div;
		}
		//=================================================================================================//
		void FreeSurfaceIndicationInner::Update(size_t index_i, Real dt)
		{
			bool is_free_surface = pos_div_[index_i] < thereshold_by_dimensions_ ? true : false;
	
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				/** Two layer particles.*/
				if (pos_div_[inner_neighborhood.j_[n]] < thereshold_by_dimensions_ && inner_neighborhood.r_ij_[n] < smoothing_length_)
				/** Three layer particles. */
				//if (pos_div_[inner_neighborhood.j_[n]] < thereshold_by_dimensions_)
				{
					is_free_surface = true;
					break;
				}
			is_free_surface_[index_i] = is_free_surface;
		}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseInnerBodyRelation* inner_relation) :
			InteractionDynamicsWithUpdate(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			rho_sum_(particles_->rho_sum_)
		{
			W0_ = particle_adaptation_->getKernel()->W0(Vecd(0));
			rho_0_ = particles_->rho_0_;
			inv_sigma_0_ = 1.0 / particles_->sigma_0_;
		}
		//=================================================================================================//
		void DensitySummationInner::Interaction(size_t index_i, Real dt)
		{
			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			rho_sum_[index_i] = sigma * rho_0_ * inv_sigma_0_;
		}
		//=================================================================================================//
		void DensitySummationInner::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] = ReinitializedDensity(rho_sum_[index_i], rho_0_, rho_n_[index_i]);
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseInnerBodyRelation* inner_relation) :
			InteractionDynamics(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_), 
			vel_n_(particles_->vel_n_), dvel_dt_others_(particles_->dvel_dt_others_)
		{
			mu_ = material_->ReferenceViscosity();
			smoothing_length_ = particle_adaptation_->ReferenceSmoothingLength();
		}		
		//=================================================================================================//
		void ViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//viscous force
				vel_derivative = (vel_i - vel_n_[index_j])
								/ (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative 
								* Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}

			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		void AngularConservativeViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			Vecd& vel_i = vel_n_[index_i];

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
				Real v_r_ij = dot(vel_i - vel_n_[index_j], r_ij * e_ij);
				Real eta_ij = 8.0 * mu_  * v_r_ij /	(r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * Vol_[index_j] / rho_i
								* inner_neighborhood.dW_ij_[n] * e_ij;
			}
	
			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseInnerBodyRelation* inner_relation) :
			InteractionDynamics(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), 
			pos_n_(particles_->pos_n_), 
			is_free_surface_(particles_->is_free_surface_), 
			p_background_(0) {}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::setupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->ReferenceDensity();
			p_background_ =  10.0 * density * speed_max * speed_max;
		}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];

			Vecd acceleration_trans(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_*Vol_[index_j] * nablaW_ij / rho_i;
			}

			/** correcting particle position */
			if(!is_free_surface_[index_i]) pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		TotalMechanicalEnergy::TotalMechanicalEnergy(FluidBody* body, Gravity* gravity)
			: ParticleDynamicsReduce<Real, ReduceSum<Real>>(body), 
			FluidDataSimple(body), mass_(particles_->mass_), 
			vel_n_(particles_->vel_n_), pos_n_(particles_->pos_n_), gravity_(gravity)
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real TotalMechanicalEnergy::ReduceFunction(size_t index_i, Real dt)
		{
			return 0.5 * mass_[index_i] * vel_n_[index_i].normSqr()
				+ mass_[index_i] * gravity_->getPotential(pos_n_[index_i]);
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(FluidBody* body)
			: ParticleDynamicsReduce<Real, ReduceMax>(body),
			FluidDataSimple(body), rho_n_(particles_->rho_n_),
			p_(particles_->p_), vel_n_(particles_->vel_n_)
		{
			smoothing_length_ = particle_adaptation_->ReferenceSmoothingLength();
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
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(FluidBody* body, Real U_max)
			: ParticleDynamicsReduce<Real, ReduceMax>(body),
			FluidDataSimple(body), vel_n_(particles_->vel_n_)
		{
			smoothing_length_ = particle_adaptation_->ReferenceSmoothingLength();
			Real rho_0 = material_->ReferenceDensity();
			Real mu = material_->ReferenceViscosity();
			Real viscous_speed = mu / rho_0 / smoothing_length_;
			Real u_max = SMAX(viscous_speed, U_max);
			initial_reference_ = u_max * u_max;
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return vel_n_[index_i].normSqr();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::OutputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			particles_->speed_max_ = speed_max;
			return 0.25 * smoothing_length_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSizeForImplicitViscosity::
			AdvectionTimeStepSizeForImplicitViscosity(FluidBody* body, Real U_max)
			: AdvectionTimeStepSize(body, U_max)
		{
			initial_reference_ = U_max * U_max;
		}
		//=================================================================================================//
		VorticityInner::
			VorticityInner(BaseInnerBodyRelation* body_inner_relation) : 
			InteractionDynamics(body_inner_relation->sph_body_),
			FluidDataInner(body_inner_relation),
			Vol_(particles_->Vol_), vel_n_(particles_->vel_n_)
		{
			SPHBody* sph_body = body_inner_relation->sph_body_;
			BaseParticles* base_particles = sph_body->base_particles_;
			//register particle variable defined in this class
			base_particles->registerAVariable<indexAngularVector, AngularVecd>(vorticity_, "VorticityInner", true);
		}
		//=================================================================================================//
		void VorticityInner::Interaction(size_t index_i, Real dt)
		{
			Vecd& vel_i = vel_n_[index_i];

			AngularVecd vorticity(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
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
		BaseRelaxation::BaseRelaxation(BaseInnerBodyRelation* inner_relation) :
			ParticleDynamics1Level(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_), rho_n_(particles_->rho_n_), 
			p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			dvel_dt_(particles_->dvel_dt_), dvel_dt_others_(particles_->dvel_dt_others_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseInnerBodyRelation* inner_relation) :
			BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BasePressureRelaxation::Initialization(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
			p_[index_i] = material_->getPressure(rho_n_[index_i]);
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BasePressureRelaxation::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		Vecd BasePressureRelaxation::computeNonConservativeAcceleration(size_t index_i)
		{
			Real rho_i = rho_n_[index_i];
			Real p_i = p_[index_i];
			Vecd acceleration = dvel_dt_others_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				Real rho_j = rho_n_[index_j];
				Real p_j = p_[index_j];

				Real p_star = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
				acceleration += (p_i - p_star) * Vol_[index_j] * dW_ij * e_ij / rho_i;
			}
			return acceleration;
		}
		//=================================================================================================//
		BaseDensityRelaxation::
			BaseDensityRelaxation(BaseInnerBodyRelation* inner_relation) :
			BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BaseDensityRelaxation::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BaseDensityRelaxation::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		PressureRelaxationRiemannInnerOldroyd_B ::
			PressureRelaxationRiemannInnerOldroyd_B(BaseInnerBodyRelation* inner_relation) : 
			PressureRelaxationRiemannInner(inner_relation),
			tau_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->tau_),
			dtau_dt_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->dtau_dt_) {}
		//=================================================================================================//
		void PressureRelaxationRiemannInnerOldroyd_B::Initialization(size_t index_i, Real dt)
		{
			PressureRelaxationRiemannInner::Initialization(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationRiemannInnerOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxationRiemannInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//elastic force
				acceleration += (tau_i + tau_[index_j]) * nablaW_ij * Vol_[index_j] / rho_i;
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		DensityRelaxationRiemannInnerOldroyd_B::
			DensityRelaxationRiemannInnerOldroyd_B(BaseInnerBodyRelation* inner_relation) :
			DensityRelaxationRiemannInner(inner_relation),
			tau_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->tau_),
			dtau_dt_(dynamic_cast<ViscoelasticFluidParticles*>(body_->base_particles_)->dtau_dt_) 
		{
			Oldroyd_B_Fluid *oldroy_b_fluid 
				= dynamic_cast<Oldroyd_B_Fluid*>(body_->base_particles_->base_material_);
			mu_p_ = oldroy_b_fluid->ReferencePolymericViscosity();
			lambda_ = oldroy_b_fluid->getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void DensityRelaxationRiemannInnerOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxationRiemannInner::Interaction(index_i, dt);
			
			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				Matd velocity_gradient = - SimTK::outer((vel_i - vel_n_[index_j]), nablaW_ij) * Vol_[index_j];
				stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient 
					- tau_i / lambda_ + (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}

			dtau_dt_[index_i] = stress_rate;
		}
		//=================================================================================================//
		void DensityRelaxationRiemannInnerOldroyd_B::Update(size_t index_i, Real dt)
		{
			DensityRelaxationRiemannInner::Update(index_i, dt);

			tau_[index_i] +=  dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		FlowRelaxationBuffer::
			FlowRelaxationBuffer(FluidBody* body, BodyPartByCell* body_part) :
			PartDynamicsByCell(body, body_part), FluidDataSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), relaxation_rate_(0.1)
		{
		};
		//=================================================================================================//
		void FlowRelaxationBuffer
			::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += 
				relaxation_rate_ * ( getTargetVelocity(pos_n_[index_i], vel_n_[index_i]) - vel_n_[index_i]);
		}
		//=================================================================================================//
		DampingBoundaryCondition::
			DampingBoundaryCondition(FluidBody* body, BodyPartByCell* body_part) :
			PartDynamicsByCell(body, body_part), FluidDataSimple(body),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), strength_(5.0)
		{
			damping_zone_bounds_ = body_part->BodyPartBounds();
		};
		//=================================================================================================//
		void DampingBoundaryCondition::Update(size_t index_i, Real dt)
		{
			Real damping_factor = (pos_n_[index_i][0] - damping_zone_bounds_.first[0]) /
								  (damping_zone_bounds_.second[0]- damping_zone_bounds_.first[0]);
			vel_n_[index_i] *=  (1.0 - dt * strength_ * damping_factor * damping_factor);
		}
		//=================================================================================================//
		StaticConfinementDensity::
			StaticConfinementDensity(FluidBody* body, NearShapeSurface* near_surface) :
			PartDynamicsByCell(body, near_surface), FluidDataSimple(body),
			rho_0_(particles_->rho_0_), inv_sigma_0_(1.0 / particles_->sigma_0_),
			mass_(particles_->mass_), rho_sum_(particles_->rho_sum_), pos_n_(particles_->pos_n_), 
			level_set_complex_shape_(near_surface->getLevelSetComplexShape()) {}
		//=================================================================================================//
		void StaticConfinementDensity::Update(size_t index_i, Real dt)
		{
			Real inv_Vol_0_i = rho_0_ / mass_[index_i];
			rho_sum_[index_i] += 
				level_set_complex_shape_->computeKernelIntegral(pos_n_[index_i]) * inv_Vol_0_i * rho_0_ * inv_sigma_0_ ;
		}
		//=================================================================================================//
		StaticConfinementPressureRelaxation::
			StaticConfinementPressureRelaxation(FluidBody* body, NearShapeSurface* near_surface) :
			PartDynamicsByCell(body, near_surface), FluidDataSimple(body),
			rho_n_(particles_->rho_n_), p_(particles_->p_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			dvel_dt_(particles_->dvel_dt_), 
			level_set_complex_shape_(near_surface->getLevelSetComplexShape()),
			riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		void StaticConfinementPressureRelaxation::Update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_complex_shape_->computeKernelGradientIntegral(pos_n_[index_i]);
			Vecd normal_to_fluid = -kernel_gradient / (kernel_gradient.norm() + TinyReal);

			FluidState state(rho_n_[index_i], vel_n_[index_i], p_[index_i]);
			Vecd vel_in_wall = -state.vel_;
			FluidState state_in_wall(rho_n_[index_i], vel_in_wall, p_[index_i]);

			//always solving one-side Riemann problem for wall boundaries
			Real p_star = riemann_solver_.getPStar(state, state_in_wall, normal_to_fluid);
			dvel_dt_[index_i] -= 2.0 * p_star * kernel_gradient / state.rho_;
		}
		//=================================================================================================//
		StaticConfinementDensityRelaxation::
			StaticConfinementDensityRelaxation(FluidBody* body, NearShapeSurface* near_surface) :
			PartDynamicsByCell(body, near_surface), FluidDataSimple(body),
			rho_n_(particles_->rho_n_), p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			level_set_complex_shape_(near_surface->getLevelSetComplexShape()),
			riemann_solver_(*material_, *material_) {}
		//=================================================================================================//
		void StaticConfinementDensityRelaxation::Update(size_t index_i, Real dt)
		{
			Vecd kernel_gradient = level_set_complex_shape_->computeKernelGradientIntegral(pos_n_[index_i]);
			Vecd normal_to_fluid = -kernel_gradient / (kernel_gradient.norm() + TinyReal);

			FluidState state(rho_n_[index_i], vel_n_[index_i], p_[index_i]);
			Vecd vel_in_wall = -state.vel_;
			FluidState state_in_wall(rho_n_[index_i], vel_in_wall, p_[index_i]);

			//always solving one-side Riemann problem for wall boundaries
			Vecd vel_star = riemann_solver_.getVStar(state, state_in_wall, normal_to_fluid);
			drho_dt_[index_i] += 2.0 * state.rho_ * dot(state.vel_ - vel_star, kernel_gradient);
		}
		//=================================================================================================//
		StaticConfinement::StaticConfinement(FluidBody* body, NearShapeSurface* near_surface) :
			density_summation_(body, near_surface),	pressure_relaxation_(body, near_surface),
			density_relaxation_(body, near_surface) {}
		//=================================================================================================//
		EmitterInflowCondition::
			EmitterInflowCondition(FluidBody* body, BodyPartByParticle* body_part) :
			PartSimpleDynamicsByParticle(body, body_part), FluidDataSimple(body),
			rho_n_(particles_->rho_n_), p_(particles_->p_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), inflow_pressure_(0)
		{
			rho_0_ = material_->ReferenceDensity();
		}
		//=================================================================================================//
		void EmitterInflowCondition
			::Update(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			vel_n_[sorted_index_i] = getTargetVelocity(pos_n_[sorted_index_i], vel_n_[sorted_index_i]);
			rho_n_[sorted_index_i] = rho_0_;
			p_[sorted_index_i] = material_->getPressure(rho_n_[sorted_index_i]);
		}
		//=================================================================================================//
		EmitterInflowInjecting
			::EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_width, int axis_direction, bool positive)
			: PartSimpleDynamicsByParticle(body, body_part), FluidDataSimple(body), pos_n_(particles_->pos_n_),
			axis_(axis_direction), periodic_translation_(0), body_buffer_width_(body_buffer_width) 
		{
			body_part_bounds_ = body_part->getBodyPartShape()->findBounds();
			periodic_translation_[axis_] = body_part_bounds_.second[axis_] - body_part_bounds_.first[axis_];
			size_t total_body_buffer_particles = body_part_particles_.size() * body_buffer_width_;
			for (size_t i = 0; i < total_body_buffer_particles; ++i)
			{
				particles_->addABufferParticle();
			}
			particles_->real_particles_bound_ += total_body_buffer_particles;
			body_->allocateConfigurationMemoriesForBufferParticles();

			checking_bound_ = positive ?
				std::bind(&EmitterInflowInjecting::checkUpperBound, this, _1, _2)
				: std::bind(&EmitterInflowInjecting::checkLowerBound, this, _1, _2);
		}
		//=================================================================================================//
		void EmitterInflowInjecting::checkUpperBound(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			if (pos_n_[sorted_index_i][axis_] > body_part_bounds_.second[axis_]) {
				if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				particles_->total_real_particles_ += 1;
				/** Periodic bounding. */
				pos_n_[sorted_index_i][axis_] -= periodic_translation_[axis_];

			}
		}
		//=================================================================================================//
		void EmitterInflowInjecting::checkLowerBound(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			if (pos_n_[sorted_index_i][axis_] < body_part_bounds_.first[axis_]) {
				if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
				{
					cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
				/** Realize the buffer particle by increasing the number of real particle in the body.  */
				particles_->total_real_particles_ += 1;
				pos_n_[sorted_index_i][axis_] += periodic_translation_[axis_];
			}
		}
		//=================================================================================================//
		SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerBodyRelation* inner_relation, Real gamma) 
		: InteractionDynamics(inner_relation->sph_body_), FluidDataInner(inner_relation), Vol_(particles_->Vol_), 
			mass_(particles_->mass_), dvel_dt_others_(particles_->dvel_dt_others_), gamma_(gamma)
		{
			color_grad_ = particles_->getVariableByName<indexVector, Vecd>("ColorGradient");
			surface_norm_ = particles_->getVariableByName<indexVector, Vecd>("SurfaceNorm");
			pos_div_ = particles_->getVariableByName<indexScalar, Real>("PositionDivergence");
		}
		//=================================================================================================//
		SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseInnerBodyRelation* inner_relation)
		: SurfaceTensionAccelerationInner(inner_relation, 1.0) {}
		//=================================================================================================//
		void SurfaceTensionAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd n_i = (*surface_norm_)[index_i];
			Real curvature(0.0);
			Real renormal_curvature(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd n_j = (*surface_norm_)[index_j];
				if(n_j.norm() > 1.0e-1)
				{
					Vecd n_ij = n_i - n_j;
					curvature -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * dot(n_ij, inner_neighborhood.e_ij_[n]);
				}
			}
			/** Normalize the curvature. */
			renormal_curvature = n_i.size() * curvature / ABS((*pos_div_)[index_i] + TinyReal);
			Vecd acceleration = gamma_ * renormal_curvature * ((*color_grad_)[index_i]) * Vol_[index_i];
			dvel_dt_others_[index_i] -= acceleration / mass_[index_i];
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//