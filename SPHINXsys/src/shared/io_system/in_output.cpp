#include "in_output.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
#include "sph_system.h"

namespace SPH {

	//===============================================================//
	In_Output::In_Output(SPHSystem &sph_system)
		: sph_system_(sph_system)
	{
		output_folder_ = "./output";

		if (fs::exists(output_folder_) && sph_system.restart_step_ == 0)
		{
			fs::remove_all(output_folder_);
		}
		fs::create_directory(output_folder_);

		restart_folder_ = "./restart";
		if (!fs::exists(restart_folder_))
		{	
			fs::create_directory(restart_folder_);
		}
		restart_step_ = std::to_string(sph_system.restart_step_);

		particle_reload_folder_ = "./particle_reload";
		if (!fs::exists(particle_reload_folder_) && sph_system.reload_particle_)
		{
			std::cout << "\n Error: the input file:" << particle_reload_folder_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//===============================================================//
	void WriteBodyStatesToVtu::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->GetBodyName() + "_" + std::to_string(Itime) + ".vtu";
			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			std::ofstream out_file(filefullpath.c_str(), ios::trunc);
			//begin of the XML file
			out_file << "<?xml version=\"1.0\"?>\n";
			out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
			out_file << " <UnstructuredGrid>\n";

			body->WriteParticlesToVtuFile(out_file);
			out_file << " </UnstructuredGrid>\n";
			out_file << "</VTKFile>\n";

			out_file.close();
		}
	}
	//===============================================================//
	void WriteBodyStatesToPlt::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->GetBodyName() + "_" + std::to_string(Itime) + ".plt";
			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			std::ofstream out_file(filefullpath.c_str(), ios::trunc);

			//begin of the plt file writing

			body->WriteParticlesToPltFile(out_file);

			out_file.close();
		}
	}
	//===============================================================//
	WriteToPltIfVelocityOutOfBound
		::WriteToPltIfVelocityOutOfBound(In_Output& in_output,
			SPHBodyVector bodies, Real velocity_bound)
		: WriteBodyStatesToPlt(in_output, bodies)
	{
		for (SPHBody* body : bodies_)
		{
			check_bodies_.push_back(new VelocityBoundCheck(body, velocity_bound));
		}
	}
	//===============================================================//
	void WriteToPltIfVelocityOutOfBound::WriteToFile(Real time)
	{
		bool out_of_bound = false;

		for (auto check_body : check_bodies_)
		{
			out_of_bound = out_of_bound || check_body->parallel_exec();
		}

		if (out_of_bound) {
			WriteBodyStatesToPlt::WriteToFile(time);
			cout << "\n Velocity is out of bound at physical time " << GlobalStaticVariables::physical_time_
				 <<"\n The body states have been outputted and the simulation terminates here. \n";
			exit(0);
		}
	}
	//===============================================================//
	WriteRelaxBodyMeshToPlt
		::WriteRelaxBodyMeshToPlt(In_Output& in_output, RelaxBody * relax_body)
		: WriteBodyStates(in_output, relax_body), relax_body_(relax_body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + relax_body->GetBodyName() + "_background_mesh.dat";
	}
	//===============================================================//
	void WriteRelaxBodyMeshToPlt::WriteToFile(Real time)
	{
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		relax_body_->mesh_background_->WriteMeshToPltFile(out_file);
		out_file.close();
	}
	//===============================================================//
	WriteObservedFluidPressure
		::WriteObservedFluidPressure(In_Output& in_output,
			ObserverBody* observer, FluidBody *interacting_body)
		: WriteBodyStates(in_output, observer), observer_(observer),
		ObserveFluidPressure(observer, interacting_body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + observer->GetBodyName() 
					  + "_fluid_pressure_" + in_output_.restart_step_ +".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "run_time" << "   ";
		for (size_t i = 0; i != observer->number_of_particles_; ++i)
		{
			out_file << "  " << "pressures_[" << i << "]" << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	void WriteObservedFluidPressure::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		parallel_exec();
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != observer_->number_of_particles_; ++i)
		{
			out_file << "  " << fluid_quantities_[i] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	WriteObservedFluidVelocity
		::WriteObservedFluidVelocity(In_Output &in_output,
			ObserverBody* observer, FluidBody *interacting_body)
		: WriteBodyStates(in_output, observer), observer_(observer),
		ObserveFluidVelocity(observer, interacting_body)
	{
		dimension_ = Vecd(0).size();
		filefullpath_ = in_output_.output_folder_ + "/" + observer->GetBodyName() 
					  + "_fluid_velocity_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "run_time" << "   ";
		for (size_t i = 0; i != observer->number_of_particles_; ++i)
		{
			for (size_t j = 0; j != dimension_; ++j)
				out_file << "  " << "velocities_[" << i << "][" << j << "]" << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	void WriteObservedFluidVelocity::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		parallel_exec();
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != observer_->number_of_particles_; ++i)
		{
			for (size_t j = 0; j != dimension_; ++j)
				out_file << "  " << fluid_quantities_[i][j] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	WriteObservedElasticDisplacement
		::WriteObservedElasticDisplacement(In_Output &in_output,
			ObserverBody* observer, SolidBody *interacting_body)
		: WriteBodyStates(in_output, observer), observer_(observer),
		ObserveElasticDisplacement(observer, interacting_body)
	{
		dimension_ = Vecd(0).size();
		filefullpath_ = in_output_.output_folder_ + "/" + observer->GetBodyName() 
			+ "_elastic_displacement_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "run_time" << "   ";
		for (size_t i = 0; i != observer->number_of_particles_; ++i)
		{
			for (size_t j = 0; j != dimension_; ++j)
				out_file << "  " << "displacements_[" << i <<"]["<< j <<"]" << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	void WriteObservedElasticDisplacement::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		parallel_exec();
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != observer_->number_of_particles_; ++i)
		{
			for (size_t j = 0; j != dimension_; ++j)
				out_file << "  " << elastic_body_quantities_[i][j] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//=================================================================================================//
	WriteObservedVoltage::WriteObservedVoltage(In_Output &in_output,
			ObserverBody* observer, SolidBody *interacting_body)
		: WriteBodyStates(in_output, observer), observer_(observer),
		ObserveMuscleVoltage(observer, interacting_body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + observer->GetBodyName() + "_voltage.dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "run_time" << "   ";
		for (size_t i = 0; i != observer->number_of_particles_; ++i)
		{
			out_file << "  " << "voltage_[" << i << "]" << " ";
		}
		out_file << "\n";
		out_file.close();
	}
//=================================================================================================//
	void WriteObservedVoltage::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		parallel_exec();
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != observer_->number_of_particles_; ++i)
		{
			out_file << "  " << muscle_quantities_[i] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
//=================================================================================================//
	WriteWaterMechanicalEnergy
		::WriteWaterMechanicalEnergy(In_Output &in_output, FluidBody* water_block, ExternalForce* external_force)
		: WriteBodyStates(in_output, water_block), TotalMechanicalEnergy(water_block, external_force)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + water_block->GetBodyName() 
					  + "_water_mechnical_energy_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\""<<"   ";
		out_file << water_block->GetBodyName()<<"   ";
		out_file <<"\n";
		out_file.close();
	};
	//===============================================================//
	void WriteWaterMechanicalEnergy::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time<<"   ";
		out_file << parallel_exec()<<"   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	WriteTotalViscousForceOnSolid
		::WriteTotalViscousForceOnSolid(In_Output &in_output, SolidBody *solid_body)
		: WriteBodyStates(in_output, solid_body), TotalViscousForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = in_output_.output_folder_ + "/total_viscous_force_on_" + solid_body->GetBodyName() + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		for(size_t i=0; i!= dimension_; ++i)
		 out_file << "\"total_force["<<i<<"]\"" << "   ";
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	void WriteTotalViscousForceOnSolid::WriteToFile(Real time)
	{
		Vecd total_force = parallel_exec();

		int Itime = int(time*1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != dimension_; ++i)
			out_file << total_force[i] << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	WriteTotalForceOnSolid
		::WriteTotalForceOnSolid(In_Output& in_output, SolidBody *solid_body)
		: WriteBodyStates(in_output, solid_body), TotalForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = in_output_.output_folder_ + "/total_force_on_" + solid_body->GetBodyName() 
			+ "_"+ in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		for(size_t i=0; i!= dimension_; ++i)
		 out_file << "\"total_force["<<i<<"]\"" << "   ";
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	void WriteTotalForceOnSolid::WriteToFile(Real time)
	{
		Vecd total_force = parallel_exec();

		int Itime = int(time*1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		for (size_t i = 0; i != dimension_; ++i)
			out_file << total_force[i] << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	WriteUpperBoundInXDirection
		::WriteUpperBoundInXDirection(In_Output& in_output, SPHBody* body)
		:WriteBodyStates(in_output, body), UpperBoundInXDirection(body)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body->GetBodyName() 
			+ "_upper_bound_in_x_direction_" + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << body->GetBodyName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	void WriteUpperBoundInXDirection::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies)
	{
		for (SPHBody* body : bodies)
		{
			file_paths_.push_back(in_output.particle_reload_folder_ + "/SPHBody_" + body->GetBodyName() + "_rld.xml");
		}
	}
	//===============================================================//
	void WriteReloadParticle::WriteToFile(Real time)
	{
		std::string reload_particle_folder = in_output_.particle_reload_folder_;
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
			bodies_[i]->WriteToXmlForReloadParticle(filefullpath);
		}
	}
	//===============================================================//
	ReadReloadParticle::ReadReloadParticle(In_Output& in_output, SPHBodyVector bodies, StdVec<std::string> reload_body_names) 
	: ReloadParticleIO(in_output, bodies), ReadBodyStates(in_output, bodies) 
	{
		if (bodies.size() != reload_body_names.size())
		{
			std::cout << "\n Error: reload bodies boes not match" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			file_paths_.push_back(in_output.particle_reload_folder_ + "/SPHBody_" + reload_body_names[i] + "_rld.xml");
		}
	}
	//===============================================================//
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

			bodies_[i]->ReadFromXmlForReloadParticle(filefullpath);
		}
	}
	//===============================================================//
	RestartIO::RestartIO(In_Output& in_output, SPHBodyVector bodies)
	{
		overall_file_path_ = in_output.restart_folder_ + "/Restart_time_";
		for (SPHBody* body : bodies)
		{
			file_paths_.push_back(in_output.restart_folder_ + "/SPHBody_" + body->GetBodyName() + "_rst_");
		}
	}
	//===============================================================//
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
			bodies_[i]->WriteParticlesToXmlForRestart(filefullpath);
		}
	}
	//===============================================================//
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
	//===============================================================//
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

			bodies_[i]->ReadParticlesFromXmlForRestart(filefullpath);
		}
	}
	//===============================================================//
	WriteSimBodyPinAngleAndAngleRate
		::WriteSimBodyPinAngleAndAngleRate(In_Output& in_output, StdVec<SimTK::MobilizedBody::Pin *> mobodies, SimTK::RungeKuttaMersonIntegrator &integ)
		: WriteSimBodyStates<SimTK::MobilizedBody::Pin>(in_output, mobodies), integ_(integ)
	{
		filefullpath_ = in_output_.output_folder_ + "/simbody_pin_angles_and_anglerate" 
					  + in_output_.restart_step_ + ".dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		for (size_t i = 0; i != mobodies_.size(); ++i)
		{
			out_file << "  " << "simbody_pin_angles_[" << i << "]" << " ";
			out_file << "  " << "simbody_pin_angle_rates_[" << i << "]" << " ";
		}
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	void WriteSimBodyPinAngleAndAngleRate::WriteToFile(Real time)
	{
		int Itime = int(time * 1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		const SimTK::State *simbody_state = &integ_.getState();
		visulizer->report(*simbody_state);
		for (size_t i = 0; i != mobodies_.size(); ++i)
		{
			out_file << "  " << mobodies_[i]->getAngle(*simbody_state) <<"  ";
			out_file << "  " << mobodies_[i]->getRate(*simbody_state) <<"  ";
		}
		out_file << "\n";
		out_file.close();
	};
}
