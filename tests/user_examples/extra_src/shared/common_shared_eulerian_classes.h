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
 *  HU1527/12-1 and HU1527/12-4													*
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
 * @file 	common_shared_eulerian_classes.h
 * @brief 	Here, we define the common shared eulerian classes for compressible and weakly compressible fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef COMMON_SHARED_EULERIAN_CLASSES_H
#define COMMON_SHARED_EULERIAN_CLASSES_H

#include "fluid_body.h"
#include "general_dynamics.h"
#include "fluid_dynamics_inner.h"
#include "compressible_fluid.h"
namespace SPH
{
	/**
	* @class EulerianFluidBody
	* @brief Eulerian Fluid body uses smoothing length to particle spacing 1.3
	*/
	class EulerianFluidBody : public FluidBody
	{
	public:
		explicit EulerianFluidBody(SPHSystem& system, SharedPtr<Shape> shape_ptr) : FluidBody(system, shape_ptr)
		{
			defineAdaptation<SPHAdaptation>(1.3);
		};
		virtual ~EulerianFluidBody() {};
		virtual EulerianFluidBody* ThisObjectPtr() override { return this; };
	};

	/**
	* @class KernalGredientWithCorrectionInner
	* @brief obtain the corrected initial configuration in strong form and correct kernel gredient
	*/
    class KernalGredientWithCorrectionInner : public LocalDynamics, public GeneralDataDelegateInner
    {
        public:
        KernalGredientWithCorrectionInner(BaseInnerRelation &inner_relation);
        virtual ~KernalGredientWithCorrectionInner(){};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

        protected:
        StdLargeVec<Matd> B_, local_configuration_inner_;
    };

	/**
	* @class KernalGredientWithCorrectionComplex
	* @brief obtain the corrected initial configuration in strong form and correct kernel gredient in complex topology
	*/
    class KernalGredientWithCorrectionComplex : public BaseInteractionComplex<KernalGredientWithCorrectionInner, GeneralDataDelegateContact>
    {
      public:
        template <typename... Args>
        KernalGredientWithCorrectionComplex(Args &&...args)
            : BaseInteractionComplex<KernalGredientWithCorrectionInner, GeneralDataDelegateContact>(std::forward<Args>(args)...){};
        virtual ~KernalGredientWithCorrectionComplex(){};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);
    };
}
#endif // COMMON_SHARED_EULERIAN_CLASSES_H