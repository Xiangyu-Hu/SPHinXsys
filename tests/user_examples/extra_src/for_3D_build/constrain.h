/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
 /**
  * @file 	constraint_dynamics.h
  * @brief 	Here, we define the algorithm classes for solid dynamics.
  * @details 	We consider here a weakly compressible solids.
  * @author	Luhui Han, Chi Zhang and Xiangyu Hu
  */

#ifndef CONSTRAIN_3d_H
#define CONSTRAIN_3d_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
//#include "elastic_solid.h"
#include "constraint_dynamics.h"
namespace SPH
{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
		//typedef DataDelegateInner<SolidParticles> SolidDataInner;


		/**
 * @class QuantityMoment
 * @brief Compute the moment of a body
 */
		//template <typename VariableType>
		class QuantityMomentOfMomentum : public QuantitySummation<Vecd>
		{
		protected:
			StdLargeVec<Real> &mass_;
			StdLargeVec<Vecd> &pos_;
			Vecd mass_center_;

		public:
			explicit QuantityMomentOfMomentum(SPHBody &sph_body, Vecd mass_center)
				: QuantitySummation<Vecd>(sph_body, "Velocity"),
				mass_center_(mass_center), mass_(this->particles_->mass_), pos_(this->particles_->pos_)
			{
				this->quantity_name_ = "Moment of Momentum";
			};
			virtual ~QuantityMomentOfMomentum() {};

			Vecd reduce(size_t index_i, Real dt = 0.0)
			{
				return (pos_[index_i] - mass_center_).cross(this->variable_[index_i]) * mass_[index_i];
			};
		};

		/**
	* @class QuantityMoment
	* @brief Compute the moment of a body
	*/
		//template <typename VariableType>
		class QuantityMomentOfInertia : public QuantitySummation<Real>
		{
		protected:
			StdLargeVec<Vecd> &pos_;
			Vecd mass_center_;
			Real p_1_;
			Real p_2_;

		public:
			explicit QuantityMomentOfInertia(SPHBody &sph_body, Vecd mass_center, Real position_1, Real position_2)
				: QuantitySummation<Real>(sph_body, "MassiveMeasure"),
				pos_(this->particles_->pos_), mass_center_(mass_center), p_1_(position_1), p_2_(position_2)
			{
				this->quantity_name_ = "Moment of Inertia";
			};
			virtual ~QuantityMomentOfInertia() {};

			Real reduce(size_t index_i, Real dt = 0.0)
			{
				if (p_1_ == p_2_)
				{
					return  ((pos_[index_i] - mass_center_).norm() * (pos_[index_i] - mass_center_).norm() -(pos_[index_i][p_1_] - mass_center_[p_1_]) * (pos_[index_i][p_2_] - mass_center_[p_2_]))* this->variable_[index_i];
				}
				else
				{
					return -(pos_[index_i][p_1_] - mass_center_[p_1_]) * (pos_[index_i][p_2_] - mass_center_[p_2_]) *this->variable_[index_i];
				}
			};
		};
		/**
	* @class QuantityMoment
	* @brief Compute the moment of a body
	*/
	//template <typename VariableType>
		class QuantityMassPosition : public QuantitySummation<Vecd>
		{
		protected:
			StdLargeVec<Real> &mass_;


		public:
			explicit QuantityMassPosition(SPHBody &sph_body)
				: QuantitySummation<Vecd>(sph_body, "Position"),
				mass_(this->particles_->mass_)
			{
				this->quantity_name_ = "Mass*Position";
			};
			virtual ~QuantityMassPosition() {};

			Vecd reduce(size_t index_i, Real dt = 0.0)
			{
				return this->variable_[index_i] * mass_[index_i];
			};
		};

		class Constrain3DSolidBodyRotation : public LocalDynamics, public SolidDataSimple
		{
		private:
			Vecd mass_center_;
			Matd moment_of_inertia_;
			Vecd angular_velocity_;
			Vecd linear_velocity_;
			ReduceDynamics<QuantityMomentOfMomentum> compute_total_moment_of_momentum_;
			StdLargeVec<Vecd> &vel_;
			StdLargeVec<Vecd> &pos_;

		protected:
			virtual void setupDynamics(Real dt = 0.0) override;

		public:
			explicit Constrain3DSolidBodyRotation(SPHBody &sph_body, Vecd mass_center, Matd inertia_tensor);
			virtual ~Constrain3DSolidBodyRotation() {};

			void update(size_t index_i, Real dt = 0.0);
		};
}
#endif //CONSTRAIN_H