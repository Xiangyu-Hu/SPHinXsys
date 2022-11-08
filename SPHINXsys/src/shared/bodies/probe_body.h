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
 * @file 	probe_body.h
 * @brief 	This is the base classes of SPH bodies. The real body is for
 *			that with cell linked list and the fictitious one does not.
 * 			Before the definition of the SPH bodies, the shapes with complex
 *			geometries, i.e. those are produced by advanced binary operation,
 * 			such as intersection, should be produced first.
 * 			Then, all shapes used in body definition should be either contain
 * 			or not contain each other.
 *			Partial overlap between them are not permitted.
 * @author	Chi ZHang and Xiangyu Hu
* @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 * @note	This body is changed from ObserverBody
 *			-- Chi ZHANG
 */

#ifndef PROBE_BODY_H
#define PROBE_BODY_H

#include "base_body.h"

namespace SPH
{
	class ProbeBody : public SPHBody
	{
	public:
		ProbeBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr);
		ProbeBody(SPHSystem &sph_system, const std::string &name);
		virtual ~ProbeBody(){};
	};
}
#endif // PROBE_BODY_H