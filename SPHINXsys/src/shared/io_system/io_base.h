/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	io_base.h
 * @brief 	base classes for io functions.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "all_physical_dynamics.h"
#include "parameterization.h"
#include "xml_engine.h"

#include <sstream>
#include <iomanip>
#include <fstream>
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH
{
	/**
	 * @class IOEnvironment
	 * @brief The base class which defines folders for output,
	 * restart and particle reload folders.
	 */
	class IOEnvironment
	{
	private:
		UniquePtrKeeper<ParameterizationIO> parameterization_io_ptr_keeper_;

	public:
		explicit IOEnvironment(SPHSystem &sph_system, bool delete_output = true);
		virtual ~IOEnvironment(){};

		SPHSystem &sph_system_;
		std::string input_folder_;
		std::string output_folder_;
		std::string restart_folder_;
		std::string reload_folder_;
		std::string restart_step_;

		ParameterizationIO &defineParameterizationIO();
	};

	/**
	 * @class BodyStatesIO
	 * @brief base class for write and read body states.
	 */
	class BodyStatesIO
	{
	protected:
		IOEnvironment &io_environment_;
		SPHBodyVector bodies_;

		std::string convertPhysicalTimeToString(Real physical_time);

	public:
		BodyStatesIO(IOEnvironment &io_environment, SPHBody &body)
			: io_environment_(io_environment), bodies_({&body}){};
		BodyStatesIO(IOEnvironment &io_environment, SPHBodyVector bodies)
			: io_environment_(io_environment), bodies_(bodies){};
		virtual ~BodyStatesIO(){};
	};

	/**
	 * @class BodyStatesRecording
	 * @brief base class for write body states.
	 */
	class BodyStatesRecording : public BodyStatesIO
	{
	public:
		BodyStatesRecording(IOEnvironment &io_environment, SPHBody &body)
			: BodyStatesIO(io_environment, body){};
		BodyStatesRecording(IOEnvironment &io_environment, SPHBodyVector bodies)
			: BodyStatesIO(io_environment, bodies){};
		virtual ~BodyStatesRecording(){};

		/** write with filename indicated by physical time */
		void writeToFile()
		{
			writeWithFileName(convertPhysicalTimeToString(GlobalStaticVariables::physical_time_));
		};

		/** write with filename indicated by iteration step */
		virtual void writeToFile(size_t iteration_step);

	protected:
		virtual void writeWithFileName(const std::string &sequence) = 0;
	};

	/**
	 * @class ReloadParticleIO
	 * @brief Write the reload particles file in XML format.
	 */
	class ReloadParticleIO : public BodyStatesIO
	{
	protected:
		StdVec<std::string> file_paths_;

	public:
		ReloadParticleIO(IOEnvironment &io_environment, SPHBodyVector bodies);
		ReloadParticleIO(IOEnvironment &io_environment, SPHBodyVector bodies, const StdVec<std::string> &given_body_names);
		virtual ~ReloadParticleIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
	};

	/**
	 * @class RestartIO
	 * @brief Write the restart file in XML format.
	 */
	class RestartIO : public BodyStatesIO
	{
	protected:
		std::string overall_file_path_;
		StdVec<std::string> file_paths_;

		Real readRestartTime(size_t restart_step);

	public:
		RestartIO(IOEnvironment &io_environment, SPHBodyVector bodies);
		virtual ~RestartIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
		virtual Real readRestartFiles(size_t restart_step)
		{
			readFromFile(restart_step);
			return readRestartTime(restart_step);
		};
	};

	/**
	 * @class ReloadMaterialParameterIO
	 * @brief For write  and read material property.
	 */
	class ReloadMaterialParameterIO
	{
	protected:
		IOEnvironment &io_environment_;
		BaseMaterial &base_material_;
		std::string file_path_;

	public:
		ReloadMaterialParameterIO(IOEnvironment &io_environment, SPHBody &sph_body);
		ReloadMaterialParameterIO(IOEnvironment &io_environment, SPHBody &sph_body,
								  const std::string &given_parameters_name);
		virtual ~ReloadMaterialParameterIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
	};
}
