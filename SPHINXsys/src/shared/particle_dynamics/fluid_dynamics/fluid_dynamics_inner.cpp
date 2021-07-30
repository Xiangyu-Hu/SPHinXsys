/**
 * @file 	fluid_dynamics.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_inner.h"
#include "in_output.h"
#include "geometry_level_set.h"

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
			FreeSurfaceIndicationInner(BaseBodyRelationInner* inner_relation, Real thereshold) :
			InteractionDynamicsWithUpdate(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			thereshold_by_dimensions_(thereshold*(Real)Dimensions), 
			Vol_(particles_->Vol_), 
			pos_div_(*particles_->createAVariable<indexScalar, Real>("PositionDivergence")),
			surface_indicator_(particles_->surface_indicator_)
		{
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
			{
				/** Two layer particles.*/
				if (pos_div_[inner_neighborhood.j_[n]] < thereshold_by_dimensions_ && inner_neighborhood.r_ij_[n] < smoothing_length_)
				{
					is_free_surface = true;
					break;
				}
			}
			surface_indicator_[index_i] = is_free_surface ? 1 : 0;
		}
		//=================================================================================================//
		MultilayeredSurfaceParticlesIdentification::
			MultilayeredSurfaceParticlesIdentification(BaseBodyRelationInner* inner_relation) :
			InteractionDynamicsWithUpdate(inner_relation->sph_body_),
			FluidDataInner(inner_relation),pos_n_(particles_->pos_n_),
			surface_indicator_(particles_->surface_indicator_),
			previous_surface_indicator_(*particles_->createAVariable<indexInteger, int>("PreviousSurfaceIndicator")) {}
		//=================================================================================================//
		void MultilayeredSurfaceParticlesIdentification::Interaction(size_t index_i, Real dt)
		{
			if (surface_indicator_[index_i] != 1)
			{
				Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					if (surface_indicator_[inner_neighborhood.j_[n]] == 1)
					{
						surface_indicator_[index_i] = 2;
						break;
					}
				}
			}
		}
		//=================================================================================================//
		void MultilayeredSurfaceParticlesIdentification::Update(size_t index_i, Real dt)
		{
			if (surface_indicator_[index_i] == 2)surface_indicator_[index_i] = 1;
			previous_surface_indicator_[index_i] = surface_indicator_[index_i];
		}
		//=================================================================================================//
		FreeStreamInletOutletSurfaceParticleIdentification  ::
			FreeStreamInletOutletSurfaceParticleIdentification(BaseBodyRelationInner* inner_relation,
				int axis_direction) :ParticleDynamicsSimple(inner_relation->sph_body_),
			FluidDataInner(inner_relation), pos_n_(particles_->pos_n_),
			fluid_body_domain_bounds_(body_->getBodyDomainBounds()), axis_(axis_direction),
			particle_spacing_(body_->particle_adaptation_->ReferenceSpacing()),
			surface_indicator_(particles_->surface_indicator_),
			previous_surface_indicator_(*particles_->createAVariable<indexInteger, int>("PreviousSurfaceIndicator")) {}
		//=================================================================================================//
		void FreeStreamInletOutletSurfaceParticleIdentification::Update(size_t index_i, Real dt)
		{
			if (surface_indicator_[index_i] == 1)
			{
				Real distance_to_inlet_bound = pos_n_[index_i][axis_] - fluid_body_domain_bounds_.first[axis_];
				Real distance_to_outlet_bound = fluid_body_domain_bounds_.second[axis_] - pos_n_[index_i][axis_];

				if (distance_to_inlet_bound < 3.5 * particle_spacing_ || distance_to_outlet_bound < 3.5 * particle_spacing_)
					surface_indicator_[index_i] = 2;
			}
			previous_surface_indicator_[index_i] = surface_indicator_[index_i];
		}
		//=================================================================================================//
		TransportVelocityCorrectionComplex::
			TransportVelocityCorrectionComplex(BaseBodyRelationInner* inner_relation,
				BaseBodyRelationContact* conatct_relation) :
			ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>(
				inner_relation, conatct_relation)
		{
			prepareContactData();
		}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseBodyRelationInner* inner_relation) :
			InteractionDynamicsWithUpdate(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			rho_sum_(particles_->rho_sum_)
		{
			W0_ = particle_adaptation_->getKernel()->W0(Vecd(0));
			rho0_ = particles_->rho0_;
			inv_sigma0_ = 1.0 / particles_->sigma0_;
		}
		//=================================================================================================//
		void DensitySummationInner::Interaction(size_t index_i, Real dt)
		{
			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
		}
		//=================================================================================================//
		void DensitySummationInner::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] = ReinitializedDensity(rho_sum_[index_i], rho0_, rho_n_[index_i]);
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		void DensitySummationFreeStreamInner::Update(size_t index_i, Real dt)
		{
			if (surface_indicator_[index_i] == 1 || surface_indicator_[index_i] == 2)
				rho_n_[index_i] = ReinitializedDensity(rho_sum_[index_i], rho0_, rho_n_[index_i]);
			else
				rho_n_[index_i] = rho_sum_[index_i];

			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseBodyRelationInner* inner_relation) :
			InteractionDynamics(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_), 
			vel_n_(particles_->vel_n_), dvel_dt_prior_(particles_->dvel_dt_prior_)
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

			dvel_dt_prior_[index_i] += acceleration;
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
	
			dvel_dt_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseBodyRelationInner* inner_relation) :
			InteractionDynamics(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), 
			pos_n_(particles_->pos_n_), 		
			surface_indicator_(particles_->surface_indicator_), p_background_(0){}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::setupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->ReferenceDensity();
			p_background_ =  7.0 * density * speed_max * speed_max;
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

			if (surface_indicator_[index_i] == 0) pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;			
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
			VorticityInner(BaseBodyRelationInner* body_inner_relation) : 
			InteractionDynamics(body_inner_relation->sph_body_),
			FluidDataInner(body_inner_relation),
			Vol_(particles_->Vol_), vel_n_(particles_->vel_n_),
			vorticity_(*particles_->createAVariable<indexAngularVector, AngularVecd>("VorticityInner")) 
		{
			particles_->addAVariableToWrite<indexAngularVector, AngularVecd>("VorticityInner");
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
		BaseRelaxation::BaseRelaxation(BaseBodyRelationInner* inner_relation) :
			ParticleDynamics1Level(inner_relation->sph_body_),
			FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_), rho_n_(particles_->rho_n_), 
			p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			dvel_dt_(particles_->dvel_dt_), dvel_dt_prior_(particles_->dvel_dt_prior_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseBodyRelationInner* inner_relation) :
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
			Vecd acceleration = dvel_dt_prior_[index_i];
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
			BaseDensityRelaxation(BaseBodyRelationInner* inner_relation) :
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
		void FreeStreamBoundaryVelocityCorrection::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
			dvel_dt_[index_i] = Vecd(0.0, 0.0);
			
			if (surface_indicator_[index_i] == 1)
			{
				Real run_time_ = GlobalStaticVariables::physical_time_;
				Real u_ave_ = run_time_ < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time_ / t_ref_)) : u_ref_;
				vel_n_[index_i][0] = u_ave_ + SMIN(rho_sum[index_i], rho_ref_) * (vel_n_[index_i][0] - u_ave_) / rho_ref_;
			}							
		}
		//=================================================================================================//
		PressureRelaxationRiemannInnerOldroyd_B ::
			PressureRelaxationRiemannInnerOldroyd_B(BaseBodyRelationInner* inner_relation) : 
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
			DensityRelaxationRiemannInnerOldroyd_B(BaseBodyRelationInner* inner_relation) :
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
			pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_), relaxation_rate_(0.3)
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
			rho0_(particles_->rho0_), inv_sigma0_(1.0 / particles_->sigma0_),
			mass_(particles_->mass_), rho_sum_(particles_->rho_sum_), pos_n_(particles_->pos_n_), 
			level_set_complex_shape_(near_surface->getLevelSetComplexShape()) {}
		//=================================================================================================//
		void StaticConfinementDensity::Update(size_t index_i, Real dt)
		{
			Real inv_Vol_0_i = rho0_ / mass_[index_i];
			rho_sum_[index_i] += 
				level_set_complex_shape_->computeKernelIntegral(pos_n_[index_i]) * inv_Vol_0_i * rho0_ * inv_sigma0_ ;
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
			rho0_ = material_->ReferenceDensity();
		}
		//=================================================================================================//
		void EmitterInflowCondition
			::Update(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			vel_n_[sorted_index_i] = getTargetVelocity(pos_n_[sorted_index_i], vel_n_[sorted_index_i]);
			rho_n_[sorted_index_i] = rho0_;
			p_[sorted_index_i] = material_->getPressure(rho_n_[sorted_index_i]);
		}
		//=================================================================================================//
		InletOutletInflowCondition::
			InletOutletInflowCondition(FluidBody* body, BodyPartByParticle* body_part) :
			EmitterInflowCondition(body, body_part){}
		//=================================================================================================//
		void InletOutletInflowCondition
			::Update(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			vel_n_[sorted_index_i] = getTargetVelocity(pos_n_[sorted_index_i], vel_n_[sorted_index_i]);
		}
		//=================================================================================================//
		EmitterInflowInjecting
			::EmitterInflowInjecting(FluidBody* body, BodyPartByParticle* body_part,
				size_t body_buffer_width, int axis_direction, bool positive)
			: PartSimpleDynamicsByParticle(body, body_part), FluidDataSimple(body), 
			pos_n_(particles_->pos_n_),rho_n_(particles_->rho_n_), p_(particles_->p_),
			axis_(axis_direction), periodic_translation_(0), body_buffer_width_(body_buffer_width)			
		{
			body_part_bounds_ = body_part->getBodyPartShape()->findBounds();
			periodic_translation_[axis_] = body_part_bounds_.second[axis_] - body_part_bounds_.first[axis_];

			size_t total_body_buffer_particles = body_part_particles_.size() * body_buffer_width_;
			particles_->addBufferParticles(total_body_buffer_particles);
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
					std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
						<< "Not enough body buffer particles! Exit the code." << "\n";
					exit(0);
				}
				/** Buffer Particle state copied from real particle. */
				particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
				/** Realize the buffer particle by increas�ng the number of real particle in the body.  */
				particles_->total_real_particles_ += 1;
				/** Periodic bounding. */
				pos_n_[sorted_index_i][axis_] -= periodic_translation_[axis_];
				rho_n_[sorted_index_i] = material_->ReferenceDensity();
				p_[sorted_index_i] = material_->getPressure(rho_n_[sorted_index_i]);
			}
		}
		//=================================================================================================//
		void EmitterInflowInjecting::checkLowerBound(size_t unsorted_index_i, Real dt)
		{
			size_t sorted_index_i = sorted_id_[unsorted_index_i];
			if (pos_n_[sorted_index_i][axis_] < body_part_bounds_.first[axis_]) {
				if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
				{
					std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
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
		ColorFunctionGradientInner::ColorFunctionGradientInner(BaseBodyRelationInner* inner_relation)
		: InteractionDynamics(inner_relation->sph_body_), FluidDataInner(inner_relation),
			Vol_(particles_->Vol_), 
			surface_indicator_(particles_->surface_indicator_),
			color_grad_(*particles_->createAVariable<indexVector, Vecd>("ColorGradient")),
			surface_norm_(*particles_->createAVariable<indexVector, Vecd>("SurfaceNormal")),
			pos_div_(*particles_->getVariableByName<indexScalar, Real>("PositionDivergence"))
		{
			//register particle variable defined in this class
			thereshold_by_dimensions_ = (0.75 * (Real)Dimensions);
		}
		//=================================================================================================//
		void ColorFunctionGradientInner::Interaction(size_t index_i, Real dt)
		{
			Vecd gradient(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			if(pos_div_[index_i] < thereshold_by_dimensions_)
			{
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					size_t index_j = inner_neighborhood.j_[n];
					gradient -= inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
				}

			}
			color_grad_[index_i] = gradient;
			surface_norm_[index_i] = gradient / (gradient.norm() + TinyReal);
		}
		//=================================================================================================//
		ColorFunctionGradientInterplationInner::ColorFunctionGradientInterplationInner(BaseBodyRelationInner* inner_relation) 
		: InteractionDynamics(inner_relation->sph_body_), FluidDataInner(inner_relation), Vol_(particles_->Vol_), 
			surface_indicator_(particles_->surface_indicator_),
			color_grad_(*particles_->getVariableByName<indexVector, Vecd>("ColorGradient")),
			surface_norm_(*particles_->getVariableByName<indexVector, Vecd>("SurfaceNormal")),
			pos_div_(*particles_->getVariableByName<indexScalar, Real>("PositionDivergence"))
		{
			thereshold_by_dimensions_ = (0.75 * (Real)Dimensions);
			particles_->addAVariableToWrite<indexVector, Vecd>("SurfaceNormal");
			particles_->addAVariableToWrite<indexVector, Vecd>("ColorGradient");
		}
		//=================================================================================================//
		void ColorFunctionGradientInterplationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd grad(0);
			Real weight(0);
			Real total_weight(0);
			if (surface_indicator_[index_i] == 1 && pos_div_[index_i] > thereshold_by_dimensions_)
			{
				Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					size_t index_j = inner_neighborhood.j_[n];
					if (surface_indicator_[index_j] == 1 && pos_div_[index_j] < thereshold_by_dimensions_)
					{
						weight = inner_neighborhood.W_ij_[n] * Vol_[index_j];
						grad += weight * color_grad_[index_j];
						total_weight += weight;
					}
				}
				Vecd grad_norm =  grad / (total_weight + TinyReal);
				color_grad_[index_i] = grad_norm;
				surface_norm_[index_i] =  grad_norm / (grad_norm.norm() + TinyReal);
			}
		}
		//=================================================================================================//
		SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseBodyRelationInner* inner_relation, Real gamma) 
		: InteractionDynamics(inner_relation->sph_body_), FluidDataInner(inner_relation),
			gamma_(gamma),  Vol_(particles_->Vol_), 
			mass_(particles_->mass_), dvel_dt_prior_(particles_->dvel_dt_prior_), surface_indicator_(particles_->surface_indicator_),
			color_grad_(*particles_->getVariableByName<indexVector, Vecd>("ColorGradient")),
			surface_norm_(*particles_->getVariableByName<indexVector, Vecd>("SurfaceNormal")){}
		//=================================================================================================//
		SurfaceTensionAccelerationInner::SurfaceTensionAccelerationInner(BaseBodyRelationInner* inner_relation)
		: SurfaceTensionAccelerationInner(inner_relation, 1.0) {}
		//=================================================================================================//
		void SurfaceTensionAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd n_i = surface_norm_[index_i];
			Real curvature(0.0);
			Real renormal_curvature(0);
			Real pos_div(0);
			if(surface_indicator_[index_i] == 1)
			{
				Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					size_t index_j = inner_neighborhood.j_[n];
					if(surface_indicator_[index_j] == 1)
					{
						Vecd n_j = surface_norm_[index_j];
						Vecd n_ij = n_i - n_j;
						curvature -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * dot(n_ij, inner_neighborhood.e_ij_[n]);
						pos_div -= inner_neighborhood.dW_ij_[n] * inner_neighborhood.r_ij_[n] * Vol_[index_j];
					}
				}
			}
			/**
			 Adami et al. 2010 is wrong in equation.
			 (dv / dt)_s = (1.0 / rho) (-sigma * k * n * delta) 
			 			 = (1/rho) * curvature * color_grad 
						 = (1/m) * curvature * color_grad * vol
			 */
			renormal_curvature = (Real)Dimensions * curvature / ABS(pos_div + TinyReal);
			Vecd acceleration = gamma_ * renormal_curvature* color_grad_[index_i] * Vol_[index_i];
			dvel_dt_prior_[index_i] -= acceleration / mass_[index_i];
		}
		//=================================================================================================//
	}		
//=================================================================================================//
}
//=================================================================================================//