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
 * @file 	fluid_surface_complex.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details 	We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_SURFACE_COMPLEX_H
#define FLUID_SURFACE_COMPLEX_H

#include "fluid_surface_inner.hpp"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
		 * @class FreeSurfaceIndicationComplex
		 * @brief indicate the particles near the free fluid surface.
		 */
		class FreeSurfaceIndicationComplex : public FreeSurfaceIndicationInner, public FluidContactData
		{
		public:
			FreeSurfaceIndicationComplex(BaseBodyRelationInner &inner_relation,
										 BaseBodyRelationContact &contact_relation, Real threshold = 0.75);
			explicit FreeSurfaceIndicationComplex(ComplexBodyRelation &complex_relation, Real threshold = 0.75);
			virtual ~FreeSurfaceIndicationComplex(){};

		protected:
			StdVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using SpatialTemporalFreeSurfaceIdentificationComplex =
			SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationComplex>;

		/** the case with free surface */
		using DensitySummationFreeSurfaceComplex = DensitySummation<DensitySummationFreeSurfaceInner>;
		/** the case with free stream */
		using DensitySummationFreeStreamComplex = DensitySummation<DensitySummationFreeStreamInner>;

		/**
		 * @class ColorFunctionGradientComplex
		 * @brief indicate the particles near the free fluid surface.
		 */
		class ColorFunctionGradientComplex : public ColorFunctionGradientInner, public FluidContactData
		{
		public:
			ColorFunctionGradientComplex(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation);
			ColorFunctionGradientComplex(ComplexBodyRelation &complex_relation);
			virtual ~ColorFunctionGradientComplex(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdVec<StdLargeVec<Real> *> contact_Vol_;
		};

		/**
		 * @class 	SurfaceNormWithWall
		 * @brief  Modify surface norm when contact with wall
		 */
		class SurfaceNormWithWall : public LocalDynamics, public FSIContactData
		{
		public:
			SurfaceNormWithWall(BaseBodyRelationContact &contact_relation, Real contact_angle);
			virtual ~SurfaceNormWithWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real contact_angle_;
			Real smoothing_length_;
			Real particle_spacing_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &surface_norm_;
			StdLargeVec<Real> &pos_div_;
			StdVec<StdLargeVec<Vecd> *> wall_n_;
		};
	}
}
#endif // FLUID_SURFACE_COMPLEX_H