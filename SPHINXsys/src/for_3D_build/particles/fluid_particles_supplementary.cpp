#include "fluid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;

namespace SPH
{
	//=================================================================================================//
	void FluidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"density\", \"Numberdensity\", \"u\", \"v\",\"w\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< pos_n_[i][2] << "  "
				<< i << "  "
				<< rho_n_[i] << " "
				<< sigma_0_[i] << " "
				<< vel_n_[i][0] << " "
				<< vel_n_[i][1] << " "
				<< vel_n_[i][2] << "\n ";
		}
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		FluidParticles::writeParticlesToPltFile(output_file);
	}
	//=================================================================================================//
}