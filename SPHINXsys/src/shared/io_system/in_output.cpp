
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
    template<typename T>
    std::string padValueWithZeros(T&& value, size_t max_string_width = 10)
	{
		std::ostringstream s_time;
		s_time << std::setw(max_string_width) << std::setfill('0') << value;
		return s_time.str();
	}
}
namespace SPH
{
	//=============================================================================================//
	In_Output::In_Output(SPHSystem& sph_system, bool delete_output)
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
			if(delete_output == true)
			{
				fs::remove_all(output_folder_);
				fs::create_directory(output_folder_);
			}
		}

		restart_step_ = padValueWithZeros(sph_system.restart_step_);

		sph_system.in_output_ = this;
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
				//TODO: we can short the file name by without using SPHBody
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".vtp";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
				//begin of the XML file
				out_file << "<?xml version=\"1.0\"?>\n";
				out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				out_file << " <PolyData>\n";

				BaseParticles *base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				out_file << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

				body->writeParticlesToVtpFile(out_file);

				out_file << "   </PointData>\n";

				//write empty cells
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
	void BodyStatesRecordingToVtpString::writeWithFileName(const std::string& sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				const auto& vtuName = body->getBodyName() + "_" + sequence + ".vtu";
				std::stringstream sstream;
				//begin of the XML file
				writeVtu(sstream, body);
				_vtuData[vtuName] = sstream.str();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtpString::writeVtu(std::ostream& stream, SPHBody* body) const
	{
		stream << "<?xml version=\"1.0\"?>\n";
		stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		stream << " <UnstructuredGrid>\n";

		BaseParticles* base_particles = body->base_particles_;
		size_t total_real_particles = base_particles->total_real_particles_;
		stream << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

		body->writeParticlesToVtuFile(stream);

		stream << "   </PointData>\n";

		//write empty cells
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
	const VtuStringData& BodyStatesRecordingToVtpString::GetVtuData() const
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
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".plt";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

				//begin of the plt file writing

				body->writeParticlesToPltFile(out_file);

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	WriteToVtpIfVelocityOutOfBound::
		WriteToVtpIfVelocityOutOfBound(In_Output &in_output, SPHBodyVector bodies, Real velocity_bound)
		: BodyStatesRecordingToVtp(in_output, bodies), out_of_bound_(false)
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(check_bodies_),
					   [&](SPHBody *body) -> VelocityBoundCheck
					   { return VelocityBoundCheck(*body, velocity_bound); });
	}
	//=============================================================================================//
	void WriteToVtpIfVelocityOutOfBound::writeWithFileName(const std::string &sequence)
	{
		for (auto check_body : check_bodies_)
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
	MeshRecordingToPlt ::MeshRecordingToPlt(In_Output &in_output, SPHBody &body, BaseMeshField *mesh_field)
		: BodyStatesRecording(in_output, body), mesh_field_(mesh_field)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body.getBodyName() + "_" + mesh_field_->Name() + ".dat";
	}
	//=============================================================================================//
	void MeshRecordingToPlt::writeWithFileName(const std::string &sequence)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		mesh_field_->writeMeshFieldToPlt(out_file);
		out_file.close();
	}
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output &in_output, SPHBodyVector bodies) : BodyStatesIO(in_output, bodies)
	{
		if (!fs::exists(in_output.reload_folder_))
		{
			fs::create_directory(in_output.reload_folder_);
		}

		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return in_output.reload_folder_ + "/" + body->getBodyName() + "_rld.xml"; });
	};
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output &in_output, SPHBodyVector bodies,
									   StdVec<std::string> given_body_names) : ReloadParticleIO(in_output, bodies)
	{
		std::transform(given_body_names.begin(), given_body_names.end(), file_paths_.begin(),
					   [&](const std::string &body_name) -> std::string
					   { return in_output.reload_folder_ + "/" + body_name + "_rld.xml"; });
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
	RestartIO::RestartIO(In_Output &in_output, SPHBodyVector bodies)
		: BodyStatesIO(in_output, bodies), overall_file_path_(in_output.restart_folder_ + "/Restart_time_")
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return in_output.restart_folder_ + "/SPHBody_" + body->getBodyName() + "_rst_"; });
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
		WriteSimBodyPinData(In_Output &in_output, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Pin>(in_output, integ, pinbody),
		filefullpath_(in_output_.output_folder_ + "/mb_pinbody_data.dat")
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
	ReloadMaterialParameterIO::ReloadMaterialParameterIO(In_Output &in_output, SharedPtr<BaseMaterial> material)
		: in_output_(in_output), material_(material.get()),
		  file_path_(in_output.reload_folder_ + "/Material_" + material->LocalParametersName() + "_rld.xml") {}
	//=================================================================================================//
	ReloadMaterialParameterIO::
		ReloadMaterialParameterIO(In_Output &in_output, SharedPtr<BaseMaterial> material, const std::string &given_parameters_name)
		: in_output_(in_output), material_(material.get()),
		  file_path_(in_output.reload_folder_ + "/Material_" + given_parameters_name + "_rld.xml") {}
	//=================================================================================================//
	void ReloadMaterialParameterIO::writeToFile(size_t iteration_step)
	{
		std::string reload_material_folder = in_output_.reload_folder_;
		if (!fs::exists(reload_material_folder))
		{
			fs::create_directory(reload_material_folder);
		}

		if (fs::exists(file_path_))
		{
			fs::remove(file_path_);
		}
		material_->writeToXmlForReloadLocalParameters(file_path_);
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
		material_->readFromXmlForLocalParameters(file_path_);
	}
	//=================================================================================================//
}
