
/**
 * @file 	in_output.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "in_output.h"
#include "all_bodies.h"
#include "level_set.h"
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
	void PltEngine::
		writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name)
	{
		out_file << "\"" << quantity_name << "\""
				 << "   ";
	}
	//=============================================================================================//
	void PltEngine::
		writeAQuantityHeader(std::ofstream &out_file, const Vecd &quantity, const std::string &quantity_name)
	{
		for (int i = 0; i != Dimensions; ++i)
			out_file << "\"" << quantity_name << "[" << i << "]\""
					 << "   ";
	}
	//=============================================================================================//
	void PltEngine::writeAQuantity(std::ofstream &out_file, const Real &quantity)
	{
		out_file << std::fixed << std::setprecision(9) << quantity << "   ";
	}
	//=============================================================================================//
	void PltEngine::writeAQuantity(std::ofstream &out_file, const Vecd &quantity)
	{
		for (int i = 0; i < Dimensions; ++i)
			out_file << std::fixed << std::setprecision(9) << quantity[i] << "   ";
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
	void BodyStatesRecordingToVtp::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				// TODO: we can short the file name by without using SPHBody
				std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtp";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
				// begin of the XML file
				out_file << "<?xml version=\"1.0\"?>\n";
				out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				out_file << " <PolyData>\n";

				BaseParticles *base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				out_file << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

				body->writeParticlesToVtpFile(out_file);

				out_file << "   </PointData>\n";

				// write empty cells
				out_file << "   <Verts>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i + 1 << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "   </Verts>\n";

				out_file << "  </Piece>\n";

				out_file << " </PolyData>\n";
				out_file << "</VTKFile>\n";

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtpString::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				const auto &vtuName = body->getName() + "_" + sequence + ".vtu";
				std::stringstream sstream;
				// begin of the XML file
				writeVtu(sstream, body);
				_vtuData[vtuName] = sstream.str();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtpString::writeVtu(std::ostream &stream, SPHBody *body) const
	{
		stream << "<?xml version=\"1.0\"?>\n";
		stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		stream << " <UnstructuredGrid>\n";

		BaseParticles *base_particles = body->base_particles_;
		size_t total_real_particles = base_particles->total_real_particles_;
		stream << "  <Piece Name =\"" << body->getName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

		body->writeParticlesToVtuFile(stream);

		stream << "   </PointData>\n";

		// write empty cells
		stream << "   <Cells>\n";
		stream << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "   </Cells>\n";

		stream << "  </Piece>\n";

		stream << " </UnstructuredGrid>\n";
		stream << "</VTKFile>\n";
	}
	//=============================================================================================//
	const VtuStringData &BodyStatesRecordingToVtpString::GetVtuData() const
	{
		return _vtuData;
	}
	//=============================================================================================//
	void BodyStatesRecordingToPlt::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".plt";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

				// begin of the plt file writing

				body->writeParticlesToPltFile(out_file);

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	WriteToVtpIfVelocityOutOfBound::
		WriteToVtpIfVelocityOutOfBound(IOEnvironment &io_environment, SPHBodyVector bodies, Real velocity_bound)
		: BodyStatesRecordingToVtp(io_environment, bodies), out_of_bound_(false)
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(check_bodies_),
					   [&](SPHBody *body) -> ReduceDynamics<VelocityBoundCheck>
					   { return ReduceDynamics<VelocityBoundCheck>(*body, velocity_bound); });
	}
	//=============================================================================================//
	void WriteToVtpIfVelocityOutOfBound::writeWithFileName(const std::string &sequence)
	{
		for (auto &check_body : check_bodies_)
		{
			out_of_bound_ = out_of_bound_ || check_body.parallel_exec();
		}

		if (out_of_bound_)
		{
			BodyStatesRecordingToVtp::writeWithFileName(sequence);
			std::cout << "\n Velocity is out of bound at iteration step " << sequence
					  << "\n The body states have been outputted and the simulation terminates here. \n";
		}
	}
	//=============================================================================================//
	MeshRecordingToPlt ::MeshRecordingToPlt(IOEnvironment &io_environment, SPHBody &body, BaseMeshField *mesh_field)
		: BodyStatesRecording(io_environment, body), mesh_field_(mesh_field)
	{
		filefullpath_ = io_environment_.output_folder_ + "/" + body.getName() + "_" + mesh_field_->Name() + ".dat";
	}
	//=============================================================================================//
	void MeshRecordingToPlt::writeWithFileName(const std::string &sequence)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		mesh_field_->writeMeshFieldToPlt(out_file);
		out_file.close();
	}
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
	//=============================================================================================//
	WriteSimBodyPinData::
		WriteSimBodyPinData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Pin>(io_environment, integ, pinbody),
		  filefullpath_(io_environment_.output_folder_ + "/mb_pinbody_data.dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "angles"
				 << " ";
		out_file << "  "
				 << "angle_rates"
				 << " ";
		out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyPinData::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		const SimTK::State &state = integ_.getState();

		out_file << "  " << mobody_.getAngle(state) << "  " << mobody_.getRate(state) << "  ";

		out_file << "\n";
		out_file.close();
	};
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
