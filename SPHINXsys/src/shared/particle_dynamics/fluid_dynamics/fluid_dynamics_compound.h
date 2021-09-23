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
* @file fluid_dynamics_compound.h
* @brief Here, we define the algorithm classes for compound fluid dynamics, 
* which is involving with fluid dynamcis inner and fluid dynamics complex.   
* @author	Chi ZHang and Xiangyu Hu
*/


#ifndef FLUID_DYNAMCIS_COMPOUND_H
#define FLUID_DYNAMCIS_COMPOUND_H

#include "fluid_dynamics_complex.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
        * @class SurfaceParticlesIndicator
		* @brief this compound class contains SpatialTemporalFreeSurfaceIdentificationComplex method,
		* @brief MultilayeredSurfaceParticlesIdentification method and FreeStreamInletOutletSurfaceParticleIdentification method.
		* @brief SpatialTemporalFreeSurfaceIdentificationComplex method detect first layer suface particles.
		* @brief MultilayeredSurfaceParticlesIdentification method detect 2nd and 3rd layer surface particles.
		* @brief FreeStreamInletOutletSurfaceParticleIdentification method can detect inlet and outlet surface particles seperately.
        * @brief you can also combine these three subclasses according to your detailed case.      
        */
		class SurfaceParticlesIndicator
		{
		protected:
			/**
			* @class FirstyLayerSurfaceParticleIndicator
			* @brief detect first layer suface particles.
			*/
			class FirstyLayerSurfaceParticleIndicator : public SpatialTemporalFreeSurfaceIdentificationComplex
			{
			public:
				FirstyLayerSurfaceParticleIndicator(BaseBodyRelationInner* inner_relation,
					BaseBodyRelationContact* contact_relation, Real thereshold = 0.75) :
					SpatialTemporalFreeSurfaceIdentificationComplex(inner_relation, contact_relation, thereshold) {};
				FirstyLayerSurfaceParticleIndicator(ComplexBodyRelation* body_complex_relation, Real thereshold = 0.75) :
					FirstyLayerSurfaceParticleIndicator(body_complex_relation->inner_relation_,
						body_complex_relation->contact_relation_, thereshold) {};
				virtual ~FirstyLayerSurfaceParticleIndicator() {};
			};

			/**
			* @class SecondAndThirdLayersSurfaceParticlesIndicator
			* @brief detect 2nd and 3rd layer surface particles at inlet and outlet.
			*/
			class SecondAndThirdLayersSurfaceParticlesIndicator : public MultilayeredSurfaceParticlesIdentification
			{
			public:
				SecondAndThirdLayersSurfaceParticlesIndicator(BaseBodyRelationInner* inner_relation) :
					MultilayeredSurfaceParticlesIdentification(inner_relation) {};
				virtual ~SecondAndThirdLayersSurfaceParticlesIndicator() {};
			};

			/**
            * @class InletOutletSurfaceParticleIndicator
            * @brief detect inlet and outlet surface particles independantly.
            */
			class InletOutletSurfaceParticleIndicator : public FreeStreamInletOutletSurfaceParticleIdentification
			{
			public:
				InletOutletSurfaceParticleIndicator(BaseBodyRelationInner* inner_relation, int axis_direction) :
					FreeStreamInletOutletSurfaceParticleIdentification(inner_relation, axis_direction) {};
				virtual ~InletOutletSurfaceParticleIndicator() {};
			};

		public:
			SurfaceParticlesIndicator(ComplexBodyRelation* body_complex_relation, BaseBodyRelationInner* inner_relation, int axis_direction = 0) :
				first_layer(body_complex_relation), second_and_third_layers(inner_relation), inlet_outlet_layers(inner_relation, axis_direction) {};
			virtual ~SurfaceParticlesIndicator() {};

			FirstyLayerSurfaceParticleIndicator first_layer;
			SecondAndThirdLayersSurfaceParticlesIndicator second_and_third_layers;
			InletOutletSurfaceParticleIndicator inlet_outlet_layers;
		};
	}
}
#endif //FLUID_DYNAMCIS_COMPOUND_H