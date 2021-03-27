/**
 * @file 	in_output.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "in_output.h"
#include "all_bodies.h"
#include "level_set.h"
#include "sph_system.h"

namespace SPH
{
	//=============================================================================================//
	In_Output::In_Output(SPHSystem& sph_system)
		: sph_system_(sph_system)
	{
		output_folder_ = sph_system.output_folder_;
		restart_folder_ = sph_system.restart_folder_;
		reload_folder_ = sph_system.reload_folder_;
		if (sph_system.restart_step_ == 0)
		{
			fs::remove_all(restart_folder_);
			fs::create_directory(restart_folder_);

			fs::remove_all(output_folder_);
			fs::create_directory(output_folder_);
		}
		restart_step_ = std::to_string(sph_system.restart_step_);
	}
	//=============================================================================================//
	void WriteBodyStatesToVtu::WriteToFile(Real time)
	{
		int Itime = int(time * 1.0e4);

		for (SPHBody* body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + std::to_string(Itime) + ".vtu";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), ios::trunc);
				//begin of the XML file
				out_file << "<?xml version=\"1.0\"?>\n";
				out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				out_file << " <UnstructuredGrid>\n";

				BaseParticles* base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				out_file << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

				body->writeParticlesToVtuFile(out_file);

				out_file << "   </PointData>\n";

				//write empty cells
				out_file << "   <Cells>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
				out_file << "    </DataArray>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
				out_file << "    </DataArray>\n";
				out_file << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
				out_file << "    </DataArray>\n";
				out_file << "   </Cells>\n";

				out_file << "  </Piece>\n";

				out_file << " </UnstructuredGrid>\n";
				out_file << "</VTKFile>\n";

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void WriteBodyStatesToPlt::WriteToFile(Real time)
	{
		for (SPHBody* body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + std::to_string(time) + ".plt";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), ios::trunc);

				//begin of the plt file writing

				body->writeParticlesToPltFile(out_file);

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	WriteToVtuIfVelocityOutOfBound
		::WriteToVtuIfVelocityOutOfBound(In_Output& in_output,
			SPHBodyVector bodies, Real velocity_bound)
		: WriteBodyStatesToVtu(in_output, bodies), out_of_bound_(false)
	{
		for (SPHBody* body : bodies_)
		{
			check_bodies_.push_back(new VelocityBoundCheck(body, velocity_bound));
		}
	}
	//=============================================================================================//
	void WriteToVtuIfVelocityOutOfBound::WriteToFile(Real time)
	{
		for (auto check_body : check_bodies_)
		{
			out_of_bound_ = out_of_bound_ || check_body->parallel_exec();
		}

		if (out_of_bound_) {
			WriteBodyStatesToVtu::WriteToFile(time);
			cout << "\n Velocity is out of bound at physical time " << GlobalStaticVariables::physical_time_
				<< "\n The body states have been outputted and the simulation terminates here. \n";
		}
	}
	//=============================================================================================//
	WriteMeshToPlt
		::WriteMeshToPlt(In_Output& in_output, SPHBody* body, Mesh* mesh)
		: WriteBodyStates(in_output, body), mesh_(mesh)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body->getBodyName() + "_" + mesh_->Name() + ".dat";
	}
	//=============================================================================================//
	void WriteMeshToPlt::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		mesh_->writeMeshToPltFile(out_file);
		out_file.close();
	}
	//=================================================================================================//
	WriteTotalMechanicalEnergy
		::WriteTotalMechanicalEnergy(In_Output& in_output, FluidBody* water_block, Gravity* gravity)
		: WriteBodyStates(in_output, water_block), TotalMechanicalEnergy(water_block, gravity)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + water_block->getBodyName()
			+ "_water_mechanical_energy_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << water_block->getBodyName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	void WriteTotalMechanicalEnergy::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
	WriteMaximumSpeed
		::WriteMaximumSpeed(In_Output& in_output, SPHBody* sph_body)
		: WriteBodyStates(in_output, sph_body), MaximumSpeed(sph_body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + sph_body->getBodyName()
			+ "_maximum_speed_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << sph_body->getBodyName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	void WriteMaximumSpeed::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	WriteTotalViscousForceOnSolid
		::WriteTotalViscousForceOnSolid(In_Output& in_output, SolidBody* solid_body)
		: WriteBodyStates(in_output, solid_body), TotalViscousForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = in_output_.output_folder_ + "/total_viscous_force_on_" + solid_body->getBodyName() + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		for (int i = 0; i != dimension_; ++i)
			out_file << "\"total_force[" << i << "]\"" << "   ";
		out_file << "\n";
		out_file.close();
	}
	//=============================================================================================//
	void WriteTotalViscousForceOnSolid::WriteToFile(Real time)
	{
		Vecd total_force = parallel_exec();

		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (int i = 0; i < dimension_; ++i)
			out_file << total_force[i] << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	WriteTotalForceOnSolid
		::WriteTotalForceOnSolid(In_Output& in_output, SolidBody* solid_body)
		: WriteBodyStates(in_output, solid_body), TotalForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = in_output_.output_folder_ + "/total_force_on_" + solid_body->getBodyName()
			+ "_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		for (int i = 0; i < dimension_; ++i)
			out_file << "\"total_force[" << i << "]\"" << "   ";
		out_file << "\n";
		out_file.close();
	}
	//=============================================================================================//
	void WriteTotalForceOnSolid::WriteToFile(Real time)
	{
		Vecd total_force = parallel_exec();

		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (int i = 0; i < dimension_; ++i)
			out_file << total_force[i] << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	WriteUpperFrontInXDirection
		::WriteUpperFrontInXDirection(In_Output& in_output, SPHBody* body)
		:WriteBodyStates(in_output, body), UpperFrontInXDirection(body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body->getBodyName()
			+ "_upper_bound_in_x_direction_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << body->getBodyName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	void WriteUpperFrontInXDirection::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies)
	{
		for (SPHBody* body : bodies)
		{
			file_paths_.push_back(in_output.reload_folder_ + "/SPHBody_" + body->getBodyName() + "_rld.xml");
		}
	}
	//=============================================================================================//
	WriteReloadParticle::WriteReloadParticle(In_Output& in_output, SPHBodyVector bodies)
		: ReloadParticleIO(in_output, bodies), WriteBodyStates(in_output, bodies)
	{
		if (!fs::exists(in_output.reload_folder_))
		{
			fs::create_directory(in_output.reload_folder_);
		}
	};
	WriteReloadParticle::WriteReloadParticle(In_Output& in_output, SPHBodyVector bodies,
		StdVec<string> given_body_names) : WriteReloadParticle(in_output, bodies)
	{
		for (size_t i = 0; i != bodies.size(); ++i)
		{
			file_paths_[i] = in_output.reload_folder_ + "/SPHBody_" + given_body_names[i] + "_rld.xml";
		}
	}
	//=============================================================================================//
	void WriteReloadParticle::WriteToFile(Real time)
	{
		std::string reload_particle_folder = in_output_.reload_folder_;
		if (!fs::exists(reload_particle_folder))
		{
			fs::create_directory(reload_particle_folder);
		}

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
	ReadReloadParticle::ReadReloadParticle(In_Output& in_output, SPHBodyVector bodies, StdVec<std::string> reload_body_names)
		: ReloadParticleIO(in_output, bodies), ReadBodyStates(in_output, bodies)
	{
		if (!fs::exists(in_output.reload_folder_))
		{
			std::cout << "\n Error: the particle reload folder:" << in_output.reload_folder_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		if (bodies.size() != reload_body_names.size())
		{
			std::cout << "\n Error: reload bodies boes not match" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=============================================================================================//
	void ReadReloadParticle::ReadFromFile(size_t restart_step)
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
	RestartIO::RestartIO(In_Output& in_output, SPHBodyVector bodies)
	{
		overall_file_path_ = in_output.restart_folder_ + "/Restart_time_";
		for (SPHBody* body : bodies)
		{
			file_paths_.push_back(in_output.restart_folder_ + "/SPHBody_" + body->getBodyName() + "_rst_");
		}
	}
	//=============================================================================================//
	void WriteRestart::WriteToFile(Real time)
	{
		int Itime = int(time);
		std::string overall_filefullpath = overall_file_path_ + std::to_string(Itime) + ".dat";
		if (fs::exists(overall_filefullpath))
		{
			fs::remove(overall_filefullpath);
		}
		std::ofstream out_file(overall_filefullpath.c_str(), ios::app);
		out_file << fixed << setprecision(9) << GlobalStaticVariables::physical_time_ << "   \n";
		out_file.close();

		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + std::to_string(Itime) + ".xml";

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			bodies_[i]->writeParticlesToXmlForRestart(filefullpath);
		}
	}
	//=============================================================================================//
	Real ReadRestart::ReadRestartTime(size_t restart_step)
	{
		std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
		std::string overall_filefullpath = overall_file_path_ + std::to_string(restart_step) + ".dat";
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
	void ReadRestart::ReadFromFile(size_t restart_step)
	{
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + std::to_string(restart_step) + ".xml";

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
	WriteSimBodyPinData::WriteSimBodyPinData(In_Output& in_output, SimTK::RungeKuttaMersonIntegrator& integ, SimTK::MobilizedBody::Pin& pinbody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Pin>(in_output, integ, pinbody)
	{
		filefullpath_ = in_output_.output_folder_ + "/mb_pinbody_data.dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);

		out_file << "\"time\"" << "   ";
		out_file << "  " << "angles" << " ";
		out_file << "  " << "angle_rates" << " ";
		out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyPinData::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		const SimTK::State& state = integ_.getState();

		out_file << "  " << mobody_.getAngle(state) << "  " << mobody_.getRate(state) << "  ";

		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
	ReloadMaterialPropertyIO::ReloadMaterialPropertyIO(In_Output& in_output, BaseMaterial* material)
	{

		file_path_ = in_output.reload_folder_ + "/Material_" + material->MaterialName() + "_rld.xml";
	}
	//=================================================================================================//
	void WriteReloadMaterialProperty::WriteToFile(Real time)
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
		material_->writeToXmlForReloadMaterialProperty(file_path_);
	}
	//=================================================================================================//
	ReadReloadMaterialProperty::ReadReloadMaterialProperty(In_Output& in_output, BaseMaterial* material)
		: ReloadMaterialPropertyIO(in_output, material), ReadMaterialProperty(in_output, material)
	{
		std::cout << "\n Reloading material property from files." << std::endl;

		if (!fs::exists(file_path_))
		{
			std::cout << "\n Error: the input file:" << file_path_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void ReadReloadMaterialProperty::ReadFromFile(size_t restart_step)
	{
		material_->readFromXmlForMaterialProperty(file_path_);
	}
	//=================================================================================================//
	WriteFreeSurfaceElevation
		::WriteFreeSurfaceElevation(In_Output& in_output, FluidBody* water_block, BodyPartByCell* body_part)
		: WriteBodyStates(in_output, water_block), FreeSurfaceProbeOnFluidBody(water_block, body_part)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body_part->BodyPartName() + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << body_part->BodyPartName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	void WriteFreeSurfaceElevation::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
}
