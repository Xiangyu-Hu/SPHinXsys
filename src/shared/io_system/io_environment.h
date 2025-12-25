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
 * @file 	io_environment.h
 * @brief 	environment for IO.
 * @author	Xiangyu Hu
 */

#ifndef IO_ENVIRONMENT_H
#define IO_ENVIRONMENT_H

#include "ownership.h"
#include "parameterization.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
namespace fs = std::filesystem;

namespace SPH
{
class SPHSystem;

/**
 * @class IOEnvironment
 * @brief The base class which defines folders for output,
 * restart and particle reload folders.
 */
class IOEnvironment
{
  private:
    UniquePtrKeeper<ParameterizationIO> parameterization_io_keeper_;

  public:
    explicit IOEnvironment(SPHSystem &sph_system);
    virtual ~IOEnvironment() {};
    void resetForRestart();
    ParameterizationIO *defineParameterizationIO();
    void appendOutputFolder(const std::string &append_name);
    void resetOutputFolder(const std::string &new_name);
    void reinitializeReloadFolder();
    std::string InputFolder() const { return input_folder_; }
    std::string OutputFolder() const { return output_folder_; }
    std::string RestartFolder() const { return restart_folder_; }
    std::string ReloadFolder() const { return reload_folder_; }

  protected:
    SPHSystem &sph_system_;
    std::string input_folder_;
    std::string output_folder_;
    std::string restart_folder_;
    std::string reload_folder_;
};
} // namespace SPH
#endif // IO_ENVIRONMENT_H
