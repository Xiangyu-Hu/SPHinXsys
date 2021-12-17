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
 * @file    fluid_body.h
 * @brief 	This is the class for bodies used for fluid.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_BODY_H
#define FLUID_BODY_H



#include "base_body.h"

namespace SPH {
	class SPHSystem;
	/**
	 * @class FluidBody
	 * @brief Fluid body uses smoothing length to particle spacing 1.3 
	 * and carry out particle sorting every 100 iterations.
	 */
	class FluidBody : public RealBody
	{
	public:
		explicit FluidBody(SPHSystem &system, const std::string &body_name,
			SharedPtr<SPHAdaptation> sph_adaptation_ptr =  makeShared<SPHAdaptation>());
		virtual ~FluidBody() {};

		/** Update cell linked list with particle sorting. */
		virtual void updateCellLinkedList() override;
		virtual FluidBody* ThisObjectPtr() override {return this;};
	protected:
		size_t iteration_count_;
	};

	/**
	 * @class EulerianFluidBody
	 * @brief Eulerian Fluid body uses smoothing length to particle spacing 1.3 
	 */
	class EulerianFluidBody : public RealBody
	{
	public:
		explicit EulerianFluidBody(SPHSystem &system, const std::string &body_name,
			SharedPtr<SPHAdaptation> sph_adaptation_ptr =  makeShared<SPHAdaptation>());
		virtual ~EulerianFluidBody() {};

		virtual EulerianFluidBody* ThisObjectPtr() override { return this; };
	};
}
#endif //FLUID_BODY_H