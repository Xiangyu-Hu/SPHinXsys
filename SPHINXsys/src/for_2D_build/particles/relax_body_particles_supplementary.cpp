/**
 * @file 	relax_body_particles.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "relax_body_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;

namespace SPH {
	//===========================================================//
	void RelaxBodyParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"Sigma\",\"ID\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].sigma_0_ << "  "
				<< i << "\n ";
		}
	}
	//===========================================================//
}
