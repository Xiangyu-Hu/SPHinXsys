/**
 * @file 	fluid_particles.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;

namespace SPH
{
	//=================================================================================================//
	void FluidParticles::WriteParticlesToPltFile(ofstream& output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\",\"NumberDensity\" , \"density\", \"u\", \"v\", \"Vorticity\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< i << "  "
				<< base_particle_data_[i].sigma_0_ << " "
				<< fluid_particle_data_[i].rho_n_ << " "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< fluid_particle_data_[i].vorticity_[0] << " \n";
		}

	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		FluidParticles::WriteParticlesToPltFile(output_file);
	}
	//=================================================================================================//
}