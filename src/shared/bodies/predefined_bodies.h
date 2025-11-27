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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    predefined_bodies.h
 * @brief 	This is the class for several predefined bodies.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PREDEFINED_BODIES_H
#define PREDEFINED_BODIES_H

#include "base_body.h"

namespace SPH
{
class FluidBody : public RealBody
{
  public:
    template <typename... Args>
    FluidBody(Args &&...args) : RealBody(std::forward<Args>(args)...){};
    virtual ~FluidBody() {};
};

class SolidBody : public RealBody
{
    void addSolidBodyToSPHSystem();

  public:
    template <typename... Args>
    SolidBody(Args &&...args)
        : RealBody(std::forward<Args>(args)...)
    {
        addSolidBodyToSPHSystem();
        defineAdaptation<SPHAdaptation>(1.15);
    };
    virtual ~SolidBody() {};
};

class ObserverBody : public SPHBody
{
    void addObserverBodyToSPHSystem();

  public:
    template <typename... Args>
    ObserverBody(Args &&...args) : SPHBody(std::forward<Args>(args)...)
    {
        addObserverBodyToSPHSystem();
    };
    virtual ~ObserverBody() {};
};
} // namespace SPH
#endif // PREDEFINED_BODIES_H