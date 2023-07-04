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
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once

#include "io_base.h"

namespace SPH
{
/**
 * @class SimBodyStatesIO
 * @brief base class for write and read SimBody states.
 */
template <class MobilizedBodyType>
class SimBodyStatesIO
{
  protected:
    IOEnvironment &io_environment_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    MobilizedBodyType &mobody_;

  public:
    SimBodyStatesIO(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, MobilizedBodyType &mobody)
        : io_environment_(io_environment), integ_(integ), mobody_(mobody){};
    virtual ~SimBodyStatesIO(){};
};

/**
 * @class WriteSimBodyStates
 * @brief base class for write SimBody states.
 */
template <class MobilizedBodyType>
class WriteSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
{
  public:
    WriteSimBodyStates(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, MobilizedBodyType &mobody)
        : SimBodyStatesIO<MobilizedBodyType>(io_environment, integ, mobody){};
    virtual ~WriteSimBodyStates(){};

    virtual void writeToFile(size_t iteration_step) = 0;
};

/**
 * @class ReadSimBodyStates
 * @brief base class for read SimBody states.
 */
template <class MobilizedBodyType>
class ReadSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
{
  public:
    ReadSimBodyStates(IOEnvironment &io_environment, MobilizedBodyType *mobody)
        : SimBodyStatesIO<MobilizedBodyType>(io_environment, mobody){};
    ReadSimBodyStates(IOEnvironment &io_environment, StdVec<MobilizedBodyType *> mobodies)
        : SimBodyStatesIO<MobilizedBodyType>(io_environment, mobodies){};
    virtual ~ReadSimBodyStates(){};

    virtual void readFromFile(size_t iteration_step) = 0;
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
    WriteSimBodyPinData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody);
    virtual ~WriteSimBodyPinData(){};
    virtual void writeToFile(size_t iteration_step = 0) override;
};
} // namespace SPH
