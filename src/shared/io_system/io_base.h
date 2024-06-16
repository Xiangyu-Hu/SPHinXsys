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
 * @file 	io_base.h
 * @brief 	base classes for io functions.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#ifndef IO_BASE_H
#define IO_BASE_H

#include "all_physical_dynamics.h"
#include "base_data_package.h"
#include "parameterization.h"
#include "sph_data_containers.h"
#include "xml_engine.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
namespace fs = std::filesystem;

namespace SPH
{
class SPHSystem;

/**
 * @class BaseIO
 * @brief base class for write and read.
 */
class BaseIO
{
  public:
    explicit BaseIO(SPHSystem &sph_system);
    virtual ~BaseIO(){};

    /** write with filename indicated by iteration step */
    virtual void writeToFile(size_t iteration_step) = 0;

  protected:
    SPHSystem &sph_system_;
    IOEnvironment &io_environment_;

    std::string convertPhysicalTimeToString(Real physical_time);

    template <typename T>
    std::string padValueWithZeros(T &&value, size_t max_string_width = 10)
    {
        std::ostringstream s_time;
        s_time << std::setw(max_string_width) << std::setfill('0') << value;
        return s_time.str();
    }
};

/**
 * @class BodyStatesRecording
 * @brief base class for write body states.
 */
class BodyStatesRecording : public BaseIO
{

  public:
    BodyStatesRecording(SPHSystem &sph_system);
    BodyStatesRecording(SPHBody &body);
    virtual ~BodyStatesRecording(){};
    /** write with filename indicated by physical time */
    void writeToFile();
    virtual void writeToFile(size_t iteration_step) override;

  protected:
    SPHBodyVector bodies_;
    bool state_recording_;

    virtual void writeWithFileName(const std::string &sequence) = 0;
};

/**
 * @class RestartIO
 * @brief Write and read the restart files in XML format.
 */
class RestartIO : public BaseIO
{
  protected:
    SPHBodyVector bodies_;
    std::string overall_file_path_;
    StdVec<std::string> file_names_;

    Real readRestartTime(size_t restart_step);

  public:
    RestartIO(SPHSystem &sph_system);
    virtual ~RestartIO(){};

    virtual void writeToFile(size_t iteration_step = 0) override;
    virtual void readFromFile(size_t iteration_step = 0);

    virtual Real readRestartFiles(size_t restart_step)
    {
        readFromFile(restart_step);
        return readRestartTime(restart_step);
    };
};

/**
 * @class ReloadParticleIO
 * @brief Write and read the particle-reloading files in XML format.
 */
class ReloadParticleIO : public BaseIO
{
  protected:
    SPHBodyVector bodies_;
    StdVec<std::string> file_names_;

  public:
    ReloadParticleIO(SPHBodyVector bodies);
    ReloadParticleIO(SPHBody &sph_body);
    ReloadParticleIO(SPHBody &sph_body, const std::string &given_body_name);
    virtual ~ReloadParticleIO(){};

    virtual void writeToFile(size_t iteration_step = 0) override;
    virtual void readFromFile(size_t iteration_step = 0);
};
} // namespace SPH
#endif // IO_BASE_H