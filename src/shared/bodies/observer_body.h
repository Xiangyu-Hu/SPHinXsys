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
 * @file 	observer_body.h
 * @brief 	This is the base classes of SPH bodies. The real body is for
 *			that with cell linked list and the fictitious one does not.
 * 			Before the definition of the SPH bodies, the shapes with complex
 *			geometries, i.e. those are produced by advanced binary operation,
 * 			such as intersection, should be produced first.
 * 			Then, all shapes used in body definition should be either contain
 * 			or not contain each other.
 *			Partial overlap between them are not permitted.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef OBSERVER_BODY_H
#define OBSERVER_BODY_H

#include "base_body.h"

namespace SPH
{
class ObserverBody : public SPHBody
{
  public:
    ObserverBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr);
    ObserverBody(SPHSystem &sph_system, const std::string &name);
    virtual ~ObserverBody(){};
};
} // namespace SPH
#endif // OBSERVER_BODY_H