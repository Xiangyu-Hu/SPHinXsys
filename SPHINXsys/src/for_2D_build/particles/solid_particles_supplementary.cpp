/**
 * @file 	solid_particles_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_particles.h"
#include "base_body.h"

using namespace std;

namespace SPH {
//=================================================================================================//
	void SolidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"v_x\", \"v_y\", \"x_norm\", \"y_norm\" \n";
		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
						<< base_particle_data_[i].pos_n_[1] << "  "
						<< base_particle_data_[i].vel_n_[0] << " "
						<< base_particle_data_[i].vel_n_[1] << " "
						<< i << "  "
						<< solid_body_data_[i].n_[0] << "  "
						<< solid_body_data_[i].n_[1] << "\n ";
		}
	}
//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		SolidParticleData& solid_body_data_i = solid_body_data_[particle_i];
		ElasticSolidParticleData &elastic_data_i = elastic_body_data_[particle_i];
		Real J = solid_body_data_i.rho_0_ / solid_body_data_i.rho_n_;
		Mat2d F = elastic_data_i.F_;
		Mat2d stress = elastic_data_i.stress_;
		Mat2d sigma = (F* stress*~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy
			+ 3.0 * sigmaxy * sigmaxy);
	}
//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\",\"v_x\", \"v_y\", \"x_norm\", \"y_norm\", \"von Mises stress\" \n";
		size_t number_of_particles = body_->number_of_particles_;

		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< i << "  "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< solid_body_data_[i].n_[0] << "  "
				<< solid_body_data_[i].n_[1] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
//=================================================================================================//
	void ActiveMuscleParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"Vx\", \"Vy\", \"von Mieses\", \" Ta \" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< i << "  "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< von_Mises_stress(i) << " "
				<< active_muscle_data_[i].active_contraction_stress_ << "\n ";
		}
	}
//=================================================================================================//
}
