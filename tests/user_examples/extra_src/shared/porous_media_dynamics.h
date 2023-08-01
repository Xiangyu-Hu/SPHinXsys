/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	POROUS_ELASTIC_DYNAMICS_H.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef POROUS_ELASTIC_DYNAMICS_H
#define POROUS_ELASTIC_DYNAMICS_H

#include "elastic_dynamics.h"
#include "porous_media_solid.h"
#include "porous_solid_particles.h"
#include "all_particle_dynamics.h"
namespace SPH
{
namespace multi_species_continuum
{
//----------------------------------------------------------------------
//		for porous elastic solid dynamics
//----------------------------------------------------------------------
    typedef DataDelegateSimple<SPH::multi_species_continuum::PorousMediaParticles> PorousMediaSolidDataSimple;
	typedef DataDelegateInner<SPH::multi_species_continuum::PorousMediaParticles> PorousMediaSolidDataInner;
 
		/** 
	 	* @class GetSaturationTimeStepSize
		* @brief Computing the time step size based on diffusion coefficient and particle smoothing length
		*/
		class GetSaturationTimeStepSize
			: public LocalDynamicsReduce<Real, ReduceMin>,
			public PorousMediaSolidDataSimple
		{
		protected:
			Real saturation_time_step_;
			Real smoothing_length_;
		public:
			explicit GetSaturationTimeStepSize(SPHBody &sph_body):
				LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)), 
				PorousMediaSolidDataSimple(sph_body)
			{ 
				smoothing_length_ = sph_body.sph_adaptation_->ReferenceSmoothingLength();
			};
			virtual ~GetSaturationTimeStepSize() {};

			Real reduce(size_t index_i, Real dt = 0.0) 
			{
				return 0.5 * smoothing_length_ * smoothing_length_ 
						/ particles_->porous_solid_.getDiffusivityConstant() / (Real)Dimensions;
			};		
	  };
     

		/**@class MomentumConstraint
		 * @brief MomentumConstraint with zero momentum.
		 */
		class MomentumConstraint : public  BaseLocalDynamics<BodyPartByParticle>, public PorousMediaSolidDataSimple
		{
		public:
			explicit MomentumConstraint(BodyPartByParticle &body_part);
			virtual ~MomentumConstraint() {};
           
		   void update(size_t index_i, Real dt = 0.0) { total_momentum_[index_i] = Vecd::Zero(); }; 

		protected:
			StdLargeVec<Vecd> &total_momentum_;
		};
 
   	/**
		 * @class BasePorousStressRelaxation
		 * @brief base class for porous media relaxation
		 */
		class BasePorousMediaRelaxation : public LocalDynamics, public PorousMediaSolidDataInner
		{
		public:
			explicit BasePorousMediaRelaxation(BaseInnerRelation &inner_relation);
			virtual ~BasePorousMediaRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &pos_, &vel_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_;
			Real rho0_, inv_rho0_;
			Real smoothing_length_;
		 
		};

    /**
		 * @class PorousMediaStressRelaxationFirstHalf
		 * @brief computing Porous Media stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class PorousMediaStressRelaxationFirstHalf
			: public BasePorousMediaRelaxation
		{
		public:
			PorousMediaStressRelaxationFirstHalf(BaseInnerRelation &body_inner_relation);
			virtual ~PorousMediaStressRelaxationFirstHalf() {};
		protected:
		
			StdLargeVec<Real> &Vol_update_, &fluid_saturation_,  &total_mass_, &fluid_mass_,  &dfluid_mass_dt_ ;
			StdLargeVec<Vecd> &total_momentum_, &dtotal_momentum_dt_, &fluid_velocity_, &relative_fluid_flux_;  
			StdLargeVec<Matd> &outer_fluid_velocity_relative_fluid_flux_, &Stress_;
		
			Real diffusivity_constant, fluid_initial_density, 
				 water_pressure_constant;

			const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;
			Real numerical_dissipation_factor_;
			Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(Vecd(0));
		

			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};

    	/**
		 * @class PorousMediaStressRelaxationSecondHalf
		 * @brief computing Porous Media stress relaxation process by verlet time stepping
		 * This is the second step
		 */
		class PorousMediaStressRelaxationSecondHalf
			: public PorousMediaStressRelaxationFirstHalf
		{
		public:
			PorousMediaStressRelaxationSecondHalf(BaseInnerRelation &body_inner_relation) :
				PorousMediaStressRelaxationFirstHalf(body_inner_relation) {};
			virtual ~PorousMediaStressRelaxationSecondHalf() {};
		protected:
			void initialization(size_t index_i, Real dt = 0.0) ;
		 	void interaction(size_t index_i, Real dt = 0.0) ;
			void update(size_t index_i, Real dt = 0.0) ;
		};
   

   	/**
		* @class PorousMediaSaturationDynamicsInitialCondition
		* @brief Set initial condition  in porous media.
		*/
		class PorousMediaSaturationDynamicsInitialCondition: public  BaseLocalDynamics<BodyPartByParticle>, public PorousMediaSolidDataSimple
		{
		public:
			  PorousMediaSaturationDynamicsInitialCondition(BodyPartByParticle &body_part):
				 BaseLocalDynamics<BodyPartByParticle>(body_part),PorousMediaSolidDataSimple(body_part.getSPHBody()),
			    fluid_mass_(particles_->fluid_mass_), fluid_saturation_(particles_->fluid_saturation_), 
				total_mass_(particles_->total_mass_), rho_n_(particles_->rho_),
				Vol_update_(particles_->Vol_update_), pos_(particles_->pos_) {};
	
			virtual ~PorousMediaSaturationDynamicsInitialCondition() {};

		protected:
			StdLargeVec<Real> &fluid_mass_, &fluid_saturation_,&total_mass_, &rho_n_, &Vol_update_;
			StdLargeVec<Vecd> &pos_;
			
		};

	  /**
		 * @class SaturationRelaxationInPorousMedia
		 * @brief computing saturation relaxation process in porous media 
		 */
		class SaturationRelaxationInPorousMedia
			: public PorousMediaStressRelaxationFirstHalf
		{
		public:
			SaturationRelaxationInPorousMedia(BaseInnerRelation &body_inner_relation)
				: PorousMediaStressRelaxationFirstHalf(body_inner_relation){};
				virtual ~SaturationRelaxationInPorousMedia() {};
		protected:
			void initialization(size_t index_i, Real Dt = 0.0) ;
		 	void interaction(size_t index_i, Real Dt = 0.0) ;
			void update(size_t index_i, Real Dt = 0.0);
		};

} // namespace multi_species_continuum
} // namespace SPH
#endif // POROUS_ELASTIC_DYNAMICS_H
