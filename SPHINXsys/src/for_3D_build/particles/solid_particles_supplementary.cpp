#include "solid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;
//=================================================================================================//
namespace SPH {
//=================================================================================================//
	void SolidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< i << "  "
				<< solid_body_data_[i].n_[0] << "  "
				<< solid_body_data_[i].n_[1] << "  "
				<< solid_body_data_[i].n_[2] << "\n ";
		}
	}
//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		ElasticSolidParticleData &elastic_data_i
			= elastic_body_data_[particle_i];
		Real J = elastic_data_i.rho_0_ / elastic_data_i.rho_n_;
		Mat3d F = elastic_data_i.F_;
		Mat3d stress = elastic_data_i.stress_;
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
	void ElasticSolidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"u\", \"v\", \"w\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\", \"von Mieses\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< base_particle_data_[i].vel_n_[0] << "  "
				<< base_particle_data_[i].vel_n_[1] << "  "
				<< base_particle_data_[i].vel_n_[2] << "  "
				<< i << "  "
				<< solid_body_data_[i].n_[0] << "  "
				<< solid_body_data_[i].n_[1] << "  "
				<< solid_body_data_[i].n_[2] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
//=================================================================================================//
	void ActiveMuscleParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"Vx\", \"Vy\", \"Vz\", \"Ta\" ,\"von Mieses \" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< i << "  "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< base_particle_data_[i].vel_n_[2] << " "
				<< active_muscle_data_[i].active_contraction_stress_ << "  "
				<< von_Mises_stress(i)              << "\n ";
		}
	}
//=================================================================================================//
}
