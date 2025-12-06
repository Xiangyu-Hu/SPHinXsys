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
 * @file 	io_base.h
 * @brief 	base classes for io functions.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#ifndef IO_BASE_H
#define IO_BASE_H

#include "base_body.h"
#include "base_data_type_package.h"
#include "base_particle_dynamics.h"
#include "parameterization.h"
#include "sphinxsys_containers.h"
#include "xml_engine.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
namespace fs = std::filesystem;

namespace SPH
{
class SPHSystem;
class IOEnvironment;

/**
 * @class BaseIO
 * @brief base class for write and read.
 */
class BaseIO
{
  public:
    explicit BaseIO(SPHSystem &sph_system);
    virtual ~BaseIO() {};

    /** write with filename indicated by iteration step */
    virtual void writeToFile(size_t iteration_step) = 0;

  protected:
    SPHSystem &sph_system_;
    IOEnvironment &io_environment_;
    SingularVariable<Real> *sv_physical_time_;

    std::string convertPhysicalTimeToString(Real physical_time);

    template <typename T>
    std::string padValueWithZeros(T &&value, size_t max_string_width = 10)
    {
        std::ostringstream s_time;
        s_time << std::setw(max_string_width) << std::setfill('0') << value;
        return s_time.str();
    }

    bool isBodyIncluded(const SPHBodyVector &bodies, SPHBody *sph_body);
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
    virtual ~BodyStatesRecording() {};
    /** write with filename indicated by physical time */
    virtual void writeToFile();
    virtual void writeToFile(size_t iteration_step) override;

    template <typename DataType>
    BodyStatesRecording &addToWrite(SPHBody &sph_body, const std::string &name)
    {
        if (isBodyIncluded(bodies_, &sph_body))
        {
            sph_body.getBaseParticles().addVariableToWrite<DataType>(name);
        }
        else
        {
            std::cout << "\n Error: the body:" << sph_body.getName()
                      << " is not in the recording list" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        return *this;
    };

    template <typename DerivedVariableMethod,
              typename DynamicsIdentifier, typename... Args>
    BodyStatesRecording &addDerivedVariableRecording(DynamicsIdentifier &identifier, Args &&...args)
    {
        SPHBody &sph_body = identifier.getSPHBody();
        if (isBodyIncluded(bodies_, &sph_body))
        {
            derived_variables_.push_back(
                derived_variables_keeper_.createPtr<DerivedVariableMethod>(
                    identifier, std::forward<Args>(args)...));
        }
        else
        {
            std::cout << "\n Error: the body:" << sph_body.getName()
                      << " is not in the recording body list" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        return *this;
    };

  protected:
    SPHBodyVector bodies_;
    StdVec<BaseDynamics<void> *> derived_variables_;
    bool state_recording_;
    virtual void writeWithFileName(const std::string &sequence) = 0;

  private:
    UniquePtrsKeeper<BaseDynamics<void>> derived_variables_keeper_;
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
    virtual ~RestartIO() {};

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
    ReloadParticleIO(SPHSystem &sph_system);
    ReloadParticleIO(SPHBodyVector bodies);
    ReloadParticleIO(SPHBody &sph_body);
    ReloadParticleIO(SPHBody &sph_body, const std::string &given_body_name);
    virtual ~ReloadParticleIO() {};

    template <typename DataType>
    void addToReload(SPHBody &sph_body, const std::string &name)
    {
        if (isBodyIncluded(bodies_, &sph_body))
        {
            sph_body.getBaseParticles().addEvolvingVariable<DataType>(name);
        }
        else
        {
            std::cout << "\n Error: the body:" << sph_body.getName()
                      << " is not in the recording list" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    };

    virtual void writeToFile(size_t iteration_step = 0) override;
};

class ParticleGenerationRecording : public BaseIO
{

  public:
    ParticleGenerationRecording(SPHBody &body);
    virtual void writeToFile(size_t iteration_step) override;

  protected:
    SPHBody &sph_body_;
    bool state_recording_;
    virtual void writeWithFileName(const std::string &sequence) = 0;
};
} // namespace SPH
#endif // IO_BASE_H
