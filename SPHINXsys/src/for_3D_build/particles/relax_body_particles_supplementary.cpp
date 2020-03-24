#include "relax_body_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;

namespace SPH {
	//===========================================================//
	void RelaxBodyParticles::WriteParticlesToPltFile(ofstream& output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"f_x\",\"f_y\", \"f_z\", \"s_x\",\"s_y\", \"s_z\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< i << " "
				<< relax_body_data_[i].f_[0] << "  "
				<< relax_body_data_[i].f_[1] << "  "
				<< relax_body_data_[i].f_[2] << "  "
				<< relax_body_data_[i].s_[0] << "  "
				<< relax_body_data_[i].s_[1] << "  "
				<< relax_body_data_[i].s_[2] << "\n";

		}
	}
	//=================================================================================================//

}
