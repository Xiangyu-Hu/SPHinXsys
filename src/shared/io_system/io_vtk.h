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
 * @file 	vtk_io.h
 * @brief 	Classes for input and output with vtk (Paraview) files.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once

#include "io_base.h"

using VtuStringData = std::map<std::string, std::string>;

namespace SPH
{
/**
 * @class BodyStatesRecordingToVtp
 * @brief  Write files for bodies
 * the output file is VTK XML format can visualized by ParaView the data type vtkPolyData
 */
class BodyStatesRecordingToVtp : public BodyStatesRecording
{
  public:
    BodyStatesRecordingToVtp(IOEnvironment &io_environment, SPHBody &body)
        : BodyStatesRecording(io_environment, body){};
    BodyStatesRecordingToVtp(IOEnvironment &io_environment, SPHBodyVector bodies)
        : BodyStatesRecording(io_environment, bodies){};
    virtual ~BodyStatesRecordingToVtp(){};

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
};

/**
 * @class BodyStatesRecordingToVtpString
 * @brief  Write strings for bodies
 * the output is map of strings with VTK XML format can visualized by ParaView
 * the data type vtkUnstructedGrid
 */
class BodyStatesRecordingToVtpString : public BodyStatesRecording
{
  public:
    BodyStatesRecordingToVtpString(IOEnvironment &io_environment, SPHBodyVector bodies)
        : BodyStatesRecording(io_environment, bodies){};
    virtual ~BodyStatesRecordingToVtpString() = default;

    const VtuStringData &GetVtuData() const;
    void clear()
    {
        _vtuData.clear();
    }

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    virtual void writeVtu(std::ostream &stream, SPHBody *body) const;

  private:
    VtuStringData _vtuData;
};

/**
 * @class WriteToVtpIfVelocityOutOfBound
 * @brief  output body sates if particle velocity is
 * out of a bound
 */
class WriteToVtpIfVelocityOutOfBound
    : public BodyStatesRecordingToVtp
{
  private:
    UniquePtrsKeeper<ReduceDynamics<VelocityBoundCheck>> check_bodies_ptr_keeper_;

  protected:
    bool out_of_bound_;
    StdVec<ReduceDynamics<VelocityBoundCheck> *> check_bodies_;
    virtual void writeWithFileName(const std::string &sequence) override;

  public:
    WriteToVtpIfVelocityOutOfBound(IOEnvironment &io_environment,
                                   SPHBodyVector bodies, Real velocity_bound);
    virtual ~WriteToVtpIfVelocityOutOfBound(){};
};
} // namespace SPH
