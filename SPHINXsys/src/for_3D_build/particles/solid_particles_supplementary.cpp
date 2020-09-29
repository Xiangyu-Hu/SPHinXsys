#include "solid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;
//=================================================================================================//
namespace SPH {
//=================================================================================================//
	void SolidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< pos_n_[i][2] << "  "
				<< i << "  "
				<< n_[i][0] << "  "
				<< n_[i][1] << "  "
				<< n_[i][2] << "\n ";
		}
	}
//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		Real J = rho_0_[particle_i] / rho_n_[particle_i];
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_[particle_i];
		Mat3d sigma = (F* stress*~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmazz = sigma(2, 2);
		Real sigmaxy = sigma(0, 1);
		Real sigmaxz = sigma(0, 2);
		Real sigmayz = sigma(1, 2);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz 
			- sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz
			+ 3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
	}
//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"u\", \"v\", \"w\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\", \"von Mieses\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file 
				<< pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< pos_n_[i][2] << "  "
				<< vel_n_[i][0] << "  "
				<< vel_n_[i][1] << "  "
				<< vel_n_[i][2] << "  "
				<< i << "  "
				<< n_[i][0] << "  "
				<< n_[i][1] << "  "
				<< n_[i][2] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
//=================================================================================================//
	void ActiveMuscleParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"Vx\", \"Vy\", \"Vz\", \"Ta\" ,\"von Mieses \" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file 
				<< pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< pos_n_[i][2] << "  "
				<< i << "  "
				<< vel_n_[i][0] << " "
				<< vel_n_[i][1] << " "
				<< vel_n_[i][2] << " "
				<< active_contraction_stress_[i] << "  "
				<< von_Mises_stress(i)              << "\n ";
		}
	}
//=================================================================================================//
}
