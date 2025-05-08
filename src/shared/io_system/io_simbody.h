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
 * @file 	io_simbody.h
 * @brief 	Classes for simbody relevant files.
 * @author	Chi Zhang, Shuoguo Zhang and Xiangyu Hu
 */

#ifndef IO_SIMBODY_H
#define IO_SIMBODY_H

#include "io_base.h"

namespace SPH
{
/**
 * @class WriteSimBodyStates
 * @brief base class for write SimBody states.
 */
template <class MobilizedBodyType>
class WriteSimBodyStates : public BaseIO
{
  public:
    WriteSimBodyStates(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                       MobilizedBodyType &mobody)
        : BaseIO(sph_system), integ_(integ), mobody_(mobody) {};
    virtual ~WriteSimBodyStates() {};

  protected:
    SimTK::RungeKuttaMersonIntegrator &integ_;
    MobilizedBodyType &mobody_;
};

/**
 * @class WriteSimBodyPinData
 * @brief Write total force acting a solid body.
 */
class WriteSimBodyPinData : public WriteSimBodyStates<SimTK::MobilizedBody::Pin>
{
  protected:
    std::string filefullpath_;

  public:
    WriteSimBodyPinData(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                        SimTK::MobilizedBody::Pin &pinbody);
    virtual ~WriteSimBodyPinData() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};

/**
 * @class WriteSimBodyCableData
 * @brief Write total force acting a single cable element.
 */
class WriteSimBodyCableData : public WriteSimBodyStates<SimTK::CableSpring>
{
  protected:
    std::string filefullpath_;

  public:
    WriteSimBodyCableData(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                          SimTK::CableSpring &cable1, std::string cable_inf);
    virtual ~WriteSimBodyCableData() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};

/**
 * @class WriteSimBodyPlanarData
 * @brief Write displacement and rotation of planar solid body.
 */
class WriteSimBodyPlanarData : public WriteSimBodyStates<SimTK::MobilizedBody::Planar>
{
  protected:
    std::string filefullpath_;

  public:
    WriteSimBodyPlanarData(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                           SimTK::MobilizedBody::Planar &planar_body);
    virtual ~WriteSimBodyPlanarData() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};

/**
 * @class WriteSimBodyFreeRotationMatrix
 * @brief Write displacement and rotation of planar solid body.
 */
class WriteSimBodyFreeRotationMatrix : public WriteSimBodyStates<SimTK::MobilizedBody::Free>
{
  protected:
    std::string filefullpath_;

  public:
    WriteSimBodyFreeRotationMatrix(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                                   SimTK::MobilizedBody::Free &free_body);
    virtual ~WriteSimBodyFreeRotationMatrix() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};

/**
 * @class WriteSimBodyVelocity
 * @brief Write displacement and rotation of planar solid body.
 */
class WriteSimBodyVelocity : public WriteSimBodyStates<SimTK::MobilizedBody::Free>
{
  protected:
    std::string filefullpath_;

  public:
    WriteSimBodyVelocity(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                         SimTK::MobilizedBody::Free &free_body);
    virtual ~WriteSimBodyVelocity() {};
    virtual void writeToFile(size_t iteration_step = 0) override;
};
} // namespace SPH
#endif // IO_SIMBODY_H
