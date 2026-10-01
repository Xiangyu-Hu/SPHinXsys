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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	simbody_system.h
 * @brief 	tbd.
 * @author	Xiangyu Hu
 */
#ifndef SIMBODY_SYSTEM_H
#define SIMBODY_SYSTEM_H

#include "base_data_type.h"
#include "ownership.h"
#include "vector_functions.h"

namespace SPH
{
class RealBody;
class Shape;

class SimbodySystem
{
  public:
    SimbodySystem();
    virtual ~SimbodySystem();
    void writeStateToXml(UnsignedInt iteration_step);
    void readStateFromXml(UnsignedInt iteration_step);
    std::string createRigidBody(RealBody &real_body, Shape &shape);
    std::string createFirstMobilizedPlanar(const std::string &body_name);
    void setUForMobilizedPlanar(
        const std::string &body_name, const Vec2d &velocity, Real angular_velocity);
    void realizeState();
    void initializeStateForIntegrator();
    void checkSimbodyState(const std::string &body_name);

  protected:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace SPH
#endif // SIMBODY_SYSTEM_H
