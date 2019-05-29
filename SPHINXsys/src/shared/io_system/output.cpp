#include "output.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
#include "sph_system.h"

namespace SPH {

	//===============================================================//
	Output::Output(SPHSystem &sph_system)
		: sph_system_(sph_system)
	{
		output_folder_ = "./output";

		if (fs::exists(output_folder_))
		{
			fs::remove_all(output_folder_);
		}

		fs::create_directory(output_folder_);

		restart_folder_ = "./rstfile";
		if (!fs::exists(restart_folder_))
		{	
			fs::create_directory(restart_folder_);
		}
	}
	//===============================================================//
	void Output::WriteCaseSetup(Real start_time, Real D_Time,  Real end_Time)
	{
		std::string filefullpath = output_folder_ + "/CaseSetup" + ".dat";
		std::ofstream out_file(filefullpath.c_str(), ios::trunc);
		out_file << "Sph_System_Particle_Spacing_Ref :" << "   " << sph_system_.particle_spacing_ref_ << "\n";
		out_file << "Sph_System_Lower_Bound :" << "   " << sph_system_.lower_bound_ << "\n";
		out_file << "Sph_System_Upper_Bound :" << "   " << sph_system_.upper_bound_ << "\n";
		out_file << "StartTime :" << "   " << start_time << "\n";
		out_file << "EndTime :" << "   " << end_Time << "\n";
		out_file << "OutputPeriod :" << "   " << D_Time << "\n";

		for (SPHBody* body : sph_system_.bodies_)
		{
			body->GlobalBasicParameters(out_file);
			out_file << "\n";
		}
		out_file.close();
	}
	//===============================================================//
	WriteBodyStates
		::WriteBodyStates(Output &output, SPHBody *body)
		: output_(output), body_(body)
	{
		output_folder_ = output_.output_folder_;
		restart_folder_ = output_.restart_folder_;
	}
	//===============================================================//
	WriteBodyStates
		::WriteBodyStates(Output &output, SPHBodyVector bodies)
		: output_(output), bodies_(bodies)
	{
		output_folder_ = output_.output_folder_;
		restart_folder_ = output_.restart_folder_;
	}
	//===============================================================//
	WriteBodyStatesToVtu
		::WriteBodyStatesToVtu(Output &output, SPHBodyVector bodies)
		: WriteBodyStates(output, bodies)
	{

	}
	//===============================================================//
	void WriteBodyStatesToVtu::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = output_folder_ + "/SPHBody_" + body->GetBodyName() + "_" + std::to_string(Itime) + ".vtu";
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
	WriteBodyStatesToPlt
		::WriteBodyStatesToPlt(Output &output, SPHBodyVector bodies)
		: WriteBodyStates(output, bodies)
	{

	}
	//===============================================================//
	void WriteBodyStatesToPlt::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = output_folder_ + "/SPHBody_" + body->GetBodyName() + "_" + std::to_string(Itime) + ".plt";
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
		::WriteToPltIfVelocityOutOfBound(Output &output,
			SPHBodyVector bodies, Real velocity_bound)
		: WriteBodyStatesToPlt(output, bodies)
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
		::WriteRelaxBodyMeshToPlt(Output &output, RelaxBody *relax_body)
		: WriteBodyStates(output, relax_body), relax_body_(relax_body)
	{
		filefullpath_ = output_folder_ + "/" + relax_body->GetBodyName() + "_background_mesh.dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file.close();

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
		::WriteObservedFluidPressure(Output &output,
			ObserverBody* observer, WeaklyCompressibleFluidBody *interacting_body)
		: WriteBodyStates(output, observer), observer_(observer),
		ObserveFluidPressure(observer, interacting_body)
	{
		filefullpath_ = output_folder_ + "/" + observer->GetBodyName() + "_fluid_pressure.dat";
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
			out_file << "  " << pressures_[i] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	WriteObservedElasticDisplacement
		::WriteObservedElasticDisplacement(Output &output,
			ObserverBody* observer, ElasticBody *interacting_body)
		: WriteBodyStates(output, observer), observer_(observer),
		ObserveElasticDisplacement(observer, interacting_body)
	{
		dimension_ = Vecd(0).size();
		filefullpath_ = output_folder_ + "/" + observer->GetBodyName() + "_elastic_displacement.dat";
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
				out_file << "  " << displacements_[i][j] << " ";
		}
		out_file << "\n";
		out_file.close();
	}
	//===============================================================//
	WriteWaterMechanicalEnergy
		::WriteWaterMechanicalEnergy(Output &output, WeaklyCompressibleFluidBody* water_block, ExternalForce* external_force)
		: WriteBodyStates(output, water_block), TotalMechanicalEnergy(water_block, external_force)
	{
		filefullpath_ = output_folder_ + "/" + water_block->GetBodyName() +"_water_mechnical_energy.dat";
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
		::WriteTotalViscousForceOnSolid(Output &output, SolidBody *solid_body)
		: WriteBodyStates(output, solid_body), TotalViscousForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = output_folder_ + "/total_viscous_force_on_" + solid_body->GetBodyName() + ".dat";
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
		::WriteTotalForceOnSolid(Output &output, SolidBody *solid_body)
		: WriteBodyStates(output, solid_body), TotalForceOnSolid(solid_body)
	{
		Vecd zero(0);
		dimension_ = zero.size();

		filefullpath_ = output_folder_ + "/total_force_on_" + solid_body->GetBodyName() + ".dat";
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
	WriteWaterFront
		::WriteWaterFront(Output &output, WeaklyCompressibleFluidBody* water_block)
		:WriteBodyStates(output, water_block), DambreakWaterFront(water_block)
	{
		filefullpath_ = output_folder_ + "/" + water_block->GetBodyName() + "_water_front.dat";
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << "\"run_time\"" << "   ";
		out_file << water_block->GetBodyName() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	void WriteWaterFront::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);
		std::ofstream out_file(filefullpath_.c_str(), ios::app);
		out_file << time << "   ";
		out_file << parallel_exec() << "   ";
		out_file << "\n";
		out_file.close();
	};
	//===============================================================//
	WriteBodyStatesToXml
		::WriteBodyStatesToXml(Output &output, SPHBodyVector bodies)
		: WriteBodyStates(output, bodies)
	{

	}
	//===============================================================//
	void WriteBodyStatesToXml::WriteToFile(Real time)
	{
		int Itime = int(time*1.0e4);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = output_folder_ + "/SPHBody_" + body->GetBodyName() + "_" + std::to_string(Itime) + ".xml";
			std::ofstream out_file(filefullpath_.c_str(), ios::app);

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			body->WriteParticlesToXmlFile(filefullpath);
		}
	}
	//===============================================================//
	WriteRestartFileToXml::WriteRestartFileToXml(Output &output, SPHBodyVector bodies)
		: WriteBodyStates(output, bodies)
	{

	}
	//===============================================================//
	void WriteRestartFileToXml::WriteToFile(Real time)
	{
		int Itime = int(time);

		for (SPHBody* body : bodies_)
		{
			std::string filefullpath = restart_folder_ + "/SPHBody_" + body->GetBodyName() + "_rst_" + std::to_string(Itime) + ".xml";

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			body->WriteParticlesToXmlForRestart(filefullpath);
		}
	}
}