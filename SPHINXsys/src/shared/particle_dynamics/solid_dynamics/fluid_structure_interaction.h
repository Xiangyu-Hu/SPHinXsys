/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	fluid_structure_interaction.h
* @brief 	Here, we define the algorithm classes for fluid structure interaction.   
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef FLUID_STRUCTURE_INTERACTION_H
#define FLUID_STRUCTURE_INTERACTION_H

#include "all_particle_dynamics.h"
#include "base_material.h"
#include "fluid_dynamics_complex.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace solid_dynamics
	{
		typedef DataDelegateSimple<SolidBody, SolidParticles, Solid> SolidDataSimple;
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, FluidBody, FluidParticles, Fluid> FSIContactData;
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, EulerianFluidBody,
			FluidParticles, Fluid> EFSIContactData; //EFSIContactData=Eulerian Fluid contact Data

		/**
		* @class FluidViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidViscousForceOnSolid : public InteractionDynamics, public FSIContactData
		{
		public:
			explicit FluidViscousForceOnSolid(BaseBodyRelationContact &contact_relation);
			virtual ~FluidViscousForceOnSolid(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_ave_;
			StdVec<StdLargeVec<Real> *> contact_Vol_, contact_rho_n_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_;
			StdVec<Real> mu_;
			StdVec<Real> smoothing_length_;
			StdLargeVec<Vecd> viscous_force_from_fluid_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidViscousForceOnSolidInEuler
		* @brief Computing the viscous force from the fluid in eulerian framework
		*/
		class FluidViscousForceOnSolidInEuler : public InteractionDynamics, public EFSIContactData
		{
		public:
			explicit FluidViscousForceOnSolidInEuler(BaseBodyRelationContact &contact_relation);
			virtual ~FluidViscousForceOnSolidInEuler(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_ave_;
			StdVec<StdLargeVec<Real> *> contact_Vol_, contact_rho_n_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_;
			StdVec<Real> mu_;
			StdVec<Real> smoothing_length_;
			StdLargeVec<Vecd> viscous_force_from_fluid_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidAngularConservativeViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidAngularConservativeViscousForceOnSolid : public FluidViscousForceOnSolid
		{
		public:
			explicit FluidAngularConservativeViscousForceOnSolid(BaseBodyRelationContact &contact_relation)
				: FluidViscousForceOnSolid(contact_relation){};
			virtual ~FluidAngularConservativeViscousForceOnSolid(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class BaseFluidPressureForce
		* @brief Template class fro computing the pressure force from the fluid with different Riemann solvers.
		* The pressrue force is added on the viscous force of the latter is computed.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		template <class RiemannSolverType>
		class BaseFluidPressureForceOnSolid : public InteractionDynamics, public FSIContactData
		{
		public:
			explicit BaseFluidPressureForceOnSolid(BaseBodyRelationContact &contact_relation)
				: InteractionDynamics(*contact_relation.sph_body_),
				  FSIContactData(contact_relation),
				  Vol_(particles_->Vol_), vel_ave_(particles_->vel_ave_),
				  force_from_fluid_(particles_->force_from_fluid_),
				  dvel_dt_ave_(particles_->dvel_dt_ave_), n_(particles_->n_)
			{
				for (size_t k = 0; k != contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
					contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
					contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
					contact_p_.push_back(&(contact_particles_[k]->p_));
					contact_dvel_dt_prior_.push_back(&(contact_particles_[k]->dvel_dt_prior_));
					riemann_solvers_.push_back(RiemannSolverType(*contact_material_[k], *contact_material_[k]));
				}
			};
			virtual ~BaseFluidPressureForceOnSolid(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_ave_, &force_from_fluid_, &dvel_dt_ave_, &n_;
			StdVec<StdLargeVec<Real> *> contact_Vol_, contact_rho_n_, contact_p_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_, contact_dvel_dt_prior_;
			StdVec<RiemannSolverType> riemann_solvers_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override
			{
				const Vecd &dvel_dt_ave_i = dvel_dt_ave_[index_i];
				Real Vol_i = Vol_[index_i];
				const Vecd &vel_ave_i = vel_ave_[index_i];
				const Vecd &n_i = n_[index_i];

				Vecd force(0);
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
					StdLargeVec<Real> &p_k = *(contact_p_[k]);
					StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
					StdLargeVec<Vecd> &dvel_dt_prior_k = *(contact_dvel_dt_prior_[k]);
					Fluid *fluid_k = contact_material_[k];
					RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Vecd e_ij = contact_neighborhood.e_ij_[n];
						Real r_ij = contact_neighborhood.r_ij_[n];
						Real face_wall_external_acceleration = dot((dvel_dt_prior_k[index_j] - dvel_dt_ave_i), e_ij);
						Real p_in_wall = p_k[index_j] + rho_n_k[index_j] * r_ij * SMAX(0.0, face_wall_external_acceleration);
						Real rho_in_wall = fluid_k->DensityFromPressure(p_in_wall);
						Vecd vel_in_wall = 2.0 * vel_ave_i - vel_n_k[index_j];

						FluidState state_l(rho_n_k[index_j], vel_n_k[index_j], p_k[index_j]);
						FluidState state_r(rho_in_wall, vel_in_wall, p_in_wall);
						Real p_star = riemann_solver_k.getPStar(state_l, state_r, n_i);
						force -= 2.0 * p_star * e_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
					}
				}
				force_from_fluid_[index_i] = force;
			};
		};
		using FluidPressureForceOnSolid = BaseFluidPressureForceOnSolid<NoRiemannSolver>;
		using FluidPressureForceOnSolidRiemann = BaseFluidPressureForceOnSolid<AcousticRiemannSolver>;

		/**
		* @class BaseFluidPressureForceOnSolidInEuler
		* @brief Template class fro computing the pressure force from the fluid with different Riemann solvers.
		* The pressrue force is added on the viscous force of the latter is computed.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		template <class RiemannSolverType>
		class BaseFluidPressureForceOnSolidInEuler : public InteractionDynamics, public EFSIContactData
		{
		public:
			explicit BaseFluidPressureForceOnSolidInEuler(BaseBodyRelationContact &contact_relation)
				: InteractionDynamics(*contact_relation.sph_body_),
				EFSIContactData(contact_relation),
				Vol_(particles_->Vol_), vel_ave_(particles_->vel_ave_),
				force_from_fluid_(particles_->force_from_fluid_), n_(particles_->n_)
			{
				for (size_t k = 0; k != contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
					contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
					contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
					contact_p_.push_back(&(contact_particles_[k]->p_));
					riemann_solvers_.push_back(RiemannSolverType(*contact_material_[k], *contact_material_[k]));
				}
			};
			virtual ~BaseFluidPressureForceOnSolidInEuler() {};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_ave_, &force_from_fluid_, &n_;
			StdVec<StdLargeVec<Real> *> contact_Vol_, contact_rho_n_, contact_p_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_;
			StdVec<RiemannSolverType> riemann_solvers_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override
			{
				Real Vol_i = Vol_[index_i];
				const Vecd &vel_ave_i = vel_ave_[index_i];
				const Vecd &n_i = n_[index_i];

				Vecd force(0);
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
					StdLargeVec<Real> &p_k = *(contact_p_[k]);
					StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
					Fluid *fluid_k = contact_material_[k];
					RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Vecd e_ij = contact_neighborhood.e_ij_[n];
						Real r_ij = contact_neighborhood.r_ij_[n];
						Real p_in_wall = p_k[index_j];
						Real rho_in_wall = fluid_k->DensityFromPressure(p_in_wall);
						Vecd vel_in_wall = - vel_n_k[index_j];

						FluidState state_l(rho_n_k[index_j], vel_n_k[index_j], p_k[index_j]);
						FluidState state_r(rho_in_wall, vel_in_wall, p_in_wall);
						FluidState interface_state = riemann_solver_k.getInterfaceState(state_l, state_r, n_i);
						Real p_star = interface_state.p_;
						force -= 2.0 * p_star * e_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
					}
				}
				force_from_fluid_[index_i] = force;
			};
		};
		using FluidPressureForceOnSolidInEuler = BaseFluidPressureForceOnSolidInEuler<NoRiemannSolver>;
		using FluidPressureForceOnSolidAcousticRiemannInEuler = BaseFluidPressureForceOnSolidInEuler<AcousticRiemannSolver>;
		using FluidPressureForceOnSolidHLLCRiemannInEuler = BaseFluidPressureForceOnSolidInEuler<HLLCRiemannSolverInWeaklyCompressibleFluid>;
		using FluidPressureForceOnSolidHLLCWithLimiterRiemannInEuler = BaseFluidPressureForceOnSolidInEuler<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		/**
		* @class BaseFluidForceOnSolidUpdate
		* @brief template class for computing force from fluid with updated viscous force
		*/
		template <class PressureForceType, class ViscousForceType>
		class BaseFluidForceOnSolidUpdate : public PressureForceType
		{
		public:
			explicit BaseFluidForceOnSolidUpdate(BaseBodyRelationContact &contact_relation)
				: PressureForceType(contact_relation), viscous_force_(contact_relation),
				  viscous_force_from_fluid_(*this->particles_->template getVariableByName<Vecd>("ViscousForceFromFluid")){};
			virtual ~BaseFluidForceOnSolidUpdate(){};

			ViscousForceType viscous_force_;

		protected:
			StdLargeVec<Vecd> &viscous_force_from_fluid_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override
			{
				PressureForceType::Interaction(index_i, dt);
				this->force_from_fluid_[index_i] += viscous_force_from_fluid_[index_i];
			};
		};
		using FluidForceOnSolidUpdate =
			BaseFluidForceOnSolidUpdate<FluidPressureForceOnSolid, FluidViscousForceOnSolid>;
		using FluidForceOnSolidUpdateRiemann =
			BaseFluidForceOnSolidUpdate<FluidPressureForceOnSolidRiemann, FluidViscousForceOnSolid>;
		using FluidForceOnSolidUpdateInEuler =
			BaseFluidForceOnSolidUpdate<FluidPressureForceOnSolidHLLCRiemannInEuler, FluidViscousForceOnSolidInEuler>;
		using FluidForceOnSolidUpdateRiemannWithLimiterInEuler =
			BaseFluidForceOnSolidUpdate<FluidPressureForceOnSolidHLLCWithLimiterRiemannInEuler, FluidViscousForceOnSolidInEuler>;


		/**
		* @class TotalViscousForceOnSolid
		* @brief Computing the total viscous force from fluid
		*/
		class TotalViscousForceOnSolid : public ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>, public SolidDataSimple
		{
		public:
			explicit TotalViscousForceOnSolid(SolidBody &solid_body);
			virtual ~TotalViscousForceOnSolid(){};

		protected:
			StdLargeVec<Vecd> &viscous_force_from_fluid_;
			Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TotalForceOnSolid
		 * @brief Computing total force from fluid.
		 */
		class TotalForceOnSolid : public ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>, public SolidDataSimple
		{
		public:
			explicit TotalForceOnSolid(SolidBody &solid_body);
			virtual ~TotalForceOnSolid(){};

		protected:
			StdLargeVec<Vecd> &force_from_fluid_;
			Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class InitializeDisplacement
		* @brief initialize the displacement for computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class InitializeDisplacement : public ParticleDynamicsSimple, public SolidDataSimple
		{
		public:
			explicit InitializeDisplacement(SolidBody &solid_body, StdLargeVec<Vecd> &pos_temp);
			virtual ~InitializeDisplacement(){};

		protected:
			StdLargeVec<Vecd> &pos_temp_, &pos_n_, &vel_ave_, &dvel_dt_ave_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class UpdateAverageVelocityAndAcceleration
		* @brief Computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class UpdateAverageVelocityAndAcceleration : public InitializeDisplacement
		{
		public:
			explicit UpdateAverageVelocityAndAcceleration(SolidBody &solid_body, StdLargeVec<Vecd> &pos_temp)
				: InitializeDisplacement(solid_body, pos_temp){};
			virtual ~UpdateAverageVelocityAndAcceleration(){};

		protected:
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AverageVelocityAndAcceleration
		* @brief Impose force matching between fluid and solid dynamics.
		* Note that the fluid time step should be larger than that of solid time step.
		* Otherwise numerical instability may occur. 
		*/
		class AverageVelocityAndAcceleration
		{
		protected:
			StdLargeVec<Vecd> pos_temp_;

		public:
			InitializeDisplacement initialize_displacement_;
			UpdateAverageVelocityAndAcceleration update_averages_;

			explicit AverageVelocityAndAcceleration(SolidBody &solid_body);
			~AverageVelocityAndAcceleration(){};
		};
	}
}
#endif //FLUID_STRUCTURE_INTERACTION_H