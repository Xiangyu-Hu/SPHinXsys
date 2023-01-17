/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	shell_fluid_interaction.h
 * @brief 	Here, we define the algorithm classes for fluid force on shell structure. 
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SHELL_FLUID_INTERACTION_H
#define SHELL_FLUID_INTERACTION_H

#include "all_particle_dynamics.h"
#include "base_material.h"
#include "fluid_dynamics_complex.h"
#include "elastic_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
	namespace solid_dynamics
	{
        typedef DataDelegateContact<ShellParticles, FluidParticles> FluidShellContactData;
		/**
		 * @class FluidViscousForceOnShell
		 * @brief Computing the viscous force on shell structure from the fluid
		 */
		class FluidViscousForceOnShell : public LocalDynamics, public FluidShellContactData
		{
		public:
			explicit FluidViscousForceOnShell(BaseContactRelation &contact_relation);
			virtual ~FluidViscousForceOnShell(){};

			void interaction(size_t index_i, Real dt = 0.0);
			StdLargeVec<Vecd> &getViscousForceFromFluid() { return viscous_force_from_fluid_; };

		protected:
            StdLargeVec<Real> &Vol_;
            StdLargeVec<Vecd> &vel_ave_;
			StdLargeVec<Vecd> viscous_force_from_fluid_;

			StdVec<Fluid *> contact_fluids_;
			StdVec<Real> mu_;
			StdVec<Real> smoothing_length_;

			StdVec<StdLargeVec<Real> *> contact_rho_n_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_;
		};

		/**
		 * @class FluidAngularConservativeViscousForceOnShell
		 * @brief Computing the viscous force on shell structure from the fluid
		 * TODO: new test for this.
		 */
		class FluidAngularConservativeViscousForceOnShell : public FluidViscousForceOnShell
		{
		public:
			explicit FluidAngularConservativeViscousForceOnShell(BaseContactRelation &contact_relation)
				: FluidViscousForceOnShell(contact_relation)
            {};
			virtual ~FluidAngularConservativeViscousForceOnShell(){};

		protected:
			void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class   BaseFluidPressureForceOnShell
		 * @brief   Template class fro computing the pressure force from the fluid with different Riemann solvers.
		 *          The pressure force is added on the viscous force of the latter is computed.
		 *          This class is for fluid-shell-interaction applications to achieve smaller solid dynamics
		 *          time step size compared to the fluid dynamics
		 */
		template <class RiemannSolverType>
		class BaseFluidPressureForceOnShell : public LocalDynamics, public FluidShellContactData
		{
		public:
			explicit BaseFluidPressureForceOnShell(BaseContactRelation &contact_relation)
				: LocalDynamics(contact_relation.sph_body_)
                , FluidShellContactData(contact_relation)
                , vel_ave_(*particles_->AverageVelocity())
                , acc_prior_(particles_->acc_prior_)
                , acc_ave_(*particles_->AverageAcceleration())
                , n_(particles_->n_)
			{
				particles_->registerVariable(force_from_fluid_, "ForceFromFluid");

				contact_fluids_.reserve(contact_particles_.size());
				riemann_solvers_.reserve(contact_particles_.size());
				contact_rho_n_.reserve(contact_particles_.size());
				contact_p_.reserve(contact_particles_.size());
				contact_vel_n_.reserve(contact_particles_.size());
				contact_acc_prior_.reserve(contact_particles_.size());

				for (const auto& cp: contact_particles_)
				{
					contact_fluids_.push_back(&cp->fluid_);
                    riemann_solvers_.push_back(RiemannSolverType(cp->fluid_, cp->fluid_));

					contact_rho_n_.push_back(&(cp->rho_));
                    contact_p_.push_back(&(cp->p_));
					contact_vel_n_.push_back(&(cp->vel_));
					contact_acc_prior_.push_back(&(cp->acc_prior_));
				}

			};
			virtual ~BaseFluidPressureForceOnShell(){};

			void interaction(size_t index_i, Real dt = 0.0) 
			{
				const Vecd &acc_ave_i = acc_ave_[index_i];
				const Vecd &vel_ave_i = vel_ave_[index_i];
				const Vecd &n_i = n_[index_i];

				Vecd force = Vecd::Zero();
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
					StdLargeVec<Real> &p_k = *(contact_p_[k]);
					StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
					StdLargeVec<Vecd> &acc_prior_k = *(contact_acc_prior_[k]);
					Fluid *fluid_k = contact_fluids_[k];
					RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Vecd e_ij = contact_neighborhood.e_ij_[n];
						Real r_ij = contact_neighborhood.r_ij_[n];
						Vecd correct_n = -SGN( e_ij.dot(n_i) ) * n_i;

						Real face_wall_external_acceleration = (acc_prior_k[index_j] - acc_ave_i).dot(e_ij);
						Real p_in_wall = p_k[index_j] + rho_n_k[index_j] * r_ij * SMAX(0.0, face_wall_external_acceleration);
						
						Vecd vel_in_shell = 2.0 * vel_ave_i - vel_n_k[index_j];
						Real u_jump = (vel_n_k[index_j] - vel_in_shell).dot(correct_n);

						force -= (p_in_wall + p_k[index_j] - riemann_solver_k.DissipativePJump(u_jump)) * e_ij 
								 * contact_neighborhood.dW_ijV_j_[n] * particles_->ParticleVolume(index_i);
					}
				}
				force_from_fluid_[index_i] = force;
				acc_prior_[index_i] = force / particles_->ParticleMass(index_i);
			};

		protected:
			StdLargeVec<Vecd> &vel_ave_, &acc_prior_, &acc_ave_, &n_;
            StdLargeVec<Vecd> force_from_fluid_;

			StdVec<Fluid *> contact_fluids_;
            StdVec<RiemannSolverType> riemann_solvers_;
			StdVec<StdLargeVec<Real> *> contact_rho_n_, contact_p_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_, contact_acc_prior_;
		};
		using FluidPressureForceOnShell = BaseFluidPressureForceOnShell<NoRiemannSolver>;
		using FluidPressureForceOnShellRiemann = BaseFluidPressureForceOnShell<AcousticRiemannSolver>;

		/**
		 * @class BaseFluidForceOnShellUpdate
		 * @brief template class for computing force on shell structure from fluid with updated viscous force
		 */
		template <class PressureForceType>
		class BaseFluidForceOnShellUpdate : public PressureForceType
		{
		public:
			template <class ViscousForceOnShellType>
			BaseFluidForceOnShellUpdate(BaseContactRelation &contact_relation, ViscousForceOnShellType &viscous_force_on_solid)
				: PressureForceType(contact_relation)
                , viscous_force_from_fluid_(viscous_force_on_solid.getViscousForceFromFluid())
            {};
			virtual ~BaseFluidForceOnShellUpdate(){};

			void interaction(size_t index_i, Real dt = 0.0)
			{
				PressureForceType::interaction(index_i, dt);
				this->force_from_fluid_[index_i] += viscous_force_from_fluid_[index_i];
				this->acc_prior_[index_i] += viscous_force_from_fluid_[index_i] / this->particles_->ParticleMass(index_i);
			};

		protected:
			StdLargeVec<Vecd> &viscous_force_from_fluid_;
		};

		using FluidForceOnShellUpdate = BaseFluidForceOnShellUpdate<FluidPressureForceOnShell>;
		using FluidForceOnShellUpdateRiemann = BaseFluidForceOnShellUpdate<FluidPressureForceOnShellRiemann>;
	}
}
#endif // FLUID_STRUCTURE_INTERACTION_H