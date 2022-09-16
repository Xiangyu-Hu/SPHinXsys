
/**
 * @file 	io_base.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "io_base.h"

#include "sph_system.h"

namespace
{
	/**
	 * @brief Convert any input to string and pad the output with zeros
	 * @todo Use external library for general string formatting, e.g. abseil, fmt library, or std::format
	 */
	template <typename T>
	std::string padValueWithZeros(T &&value, size_t max_string_width = 10)
	{
		std::ostringstream s_time;
		s_time << std::setw(max_string_width) << std::setfill('0') << value;
		return s_time.str();
	}
}
namespace SPH
{
	//=============================================================================================//
	IOEnvironment::IOEnvironment(SPHSystem &sph_system, bool delete_output)
		: sph_system_(sph_system),
		  input_folder_("./input"), output_folder_("./output"),
		  restart_folder_("./restart"), reload_folder_("./reload")
	{
		if (!fs::exists(input_folder_))
		{
			fs::create_directory(input_folder_);
		}

		if (!fs::exists(output_folder_))
		{
			fs::create_directory(output_folder_);
		}

		if (!fs::exists(restart_folder_))
		{
			fs::create_directory(restart_folder_);
		}

		if (sph_system.restart_step_ == 0)
		{
			fs::remove_all(restart_folder_);
			fs::create_directory(restart_folder_);
			if (delete_output == true)
			{
				fs::remove_all(output_folder_);
				fs::create_directory(output_folder_);
			}
		}

		restart_step_ = padValueWithZeros(sph_system.restart_step_);

		sph_system.io_environment_ = this;
	}
	//=============================================================================================//
	ParameterizationIO &IOEnvironment::defineParameterizationIO()
	{
		return parameterization_io_ptr_keeper_.createRef<ParameterizationIO>(input_folder_);
	}
	//=============================================================================================//
	std::string BodyStatesIO::convertPhysicalTimeToString(Real convertPhysicalTimeToStream)
	{
		int i_time = int(GlobalStaticVariables::physical_time_ * 1.0e6);
		return padValueWithZeros(i_time);
	}
	//=============================================================================================//
	void BodyStatesRecording::writeToFile(size_t iteration_step)
	{
		writeWithFileName(padValueWithZeros(iteration_step));
	};
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(IOEnvironment &io_environment, SPHBodyVector bodies) : BodyStatesIO(io_environment, bodies)
	{
		if (!fs::exists(io_environment.reload_folder_))
		{
			fs::create_directory(io_environment.reload_folder_);
		}

		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return io_environment.reload_folder_ + "/" + body->getName() + "_rld.xml"; });
	};
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(IOEnvironment &io_environment, SPHBodyVector bodies,
									   const StdVec<std::string> &given_body_names)
		: ReloadParticleIO(io_environment, bodies)
	{
		std::transform(given_body_names.begin(), given_body_names.end(), file_paths_.begin(),
					   [&](const std::string &body_name) -> std::string
					   { return io_environment.reload_folder_ + "/" + body_name + "_rld.xml"; });
	}
	//=============================================================================================//
	void ReloadParticleIO::writeToFile(size_t iteration_step)
	{
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i];

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			bodies_[i]->writeToXmlForReloadParticle(filefullpath);
		}
	}
	//=============================================================================================//
	void ReloadParticleIO::readFromFile(size_t restart_step)
	{
		std::cout << "\n Reloading particles from files." << std::endl;
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i];

			if (!fs::exists(filefullpath))
			{
				std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			bodies_[i]->readFromXmlForReloadParticle(filefullpath);
		}
	}
	//=============================================================================================//
	RestartIO::RestartIO(IOEnvironment &io_environment, SPHBodyVector bodies)
		: BodyStatesIO(io_environment, bodies), overall_file_path_(io_environment.restart_folder_ + "/Restart_time_")
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return io_environment.restart_folder_ + "/SPHBody_" + body->getName() + "_rst_"; });
	}
	//=============================================================================================//
	void RestartIO::writeToFile(size_t iteration_step)
	{
		std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(iteration_step) + ".dat";
		if (fs::exists(overall_filefullpath))
		{
			fs::remove(overall_filefullpath);
		}
		std::ofstream out_file(overall_filefullpath.c_str(), std::ios::app);
		out_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "   \n";
		out_file.close();

		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + padValueWithZeros(iteration_step) + ".xml";

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			bodies_[i]->writeParticlesToXmlForRestart(filefullpath);
		}
	}
	//=============================================================================================//
	Real RestartIO::readRestartTime(size_t restart_step)
	{
		std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
		std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".dat";
		if (!fs::exists(overall_filefullpath))
		{
			std::cout << "\n Error: the input file:" << overall_filefullpath << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		Real restart_time;
		std::ifstream in_file(overall_filefullpath.c_str());
		in_file >> restart_time;
		in_file.close();

		return restart_time;
	}
	//=============================================================================================//
	void RestartIO::readFromFile(size_t restart_step)
	{
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + padValueWithZeros(restart_step) + ".xml";

			if (!fs::exists(filefullpath))
			{
				std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			bodies_[i]->readParticlesFromXmlForRestart(filefullpath);
		}
	}
	//=================================================================================================//
	ReloadMaterialParameterIO::ReloadMaterialParameterIO(IOEnvironment &io_environment, BaseMaterial *base_material)
		: io_environment_(io_environment), base_material_(base_material),
		  file_path_(io_environment.reload_folder_ + "/Material_" + base_material->LocalParametersName() + "_rld.xml") {}
	//=================================================================================================//
	ReloadMaterialParameterIO::
		ReloadMaterialParameterIO(IOEnvironment &io_environment, BaseMaterial *base_material, const std::string &given_parameters_name)
		: io_environment_(io_environment), base_material_(base_material),
		  file_path_(io_environment.reload_folder_ + "/Material_" + given_parameters_name + "_rld.xml") {}
	//=================================================================================================//
	void ReloadMaterialParameterIO::writeToFile(size_t iteration_step)
	{
		std::string reload_material_folder = io_environment_.reload_folder_;
		if (!fs::exists(reload_material_folder))
		{
			fs::create_directory(reload_material_folder);
		}

		if (fs::exists(file_path_))
		{
			fs::remove(file_path_);
		}
		base_material_->writeToXmlForReloadLocalParameters(file_path_);
	}
	//=================================================================================================//
	void ReloadMaterialParameterIO::readFromFile(size_t restart_step)
	{
		if (!fs::exists(file_path_))
		{
			std::cout << "\n Error: the reloading material property file:" << file_path_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		base_material_->readFromXmlForLocalParameters(file_path_);
	}
	//=================================================================================================//
}
