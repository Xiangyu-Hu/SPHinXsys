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
	void SolidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"v_x\", \"v_y\", \"x_norm\", \"y_norm\" \n";
		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << pos_n_[i][0] << "  "
						<< pos_n_[i][1] << "  "
						<< vel_n_[i][0] << " "
						<< vel_n_[i][1] << " "
						<< i << "  "
						<< n_[i][0] << "  "
						<< n_[i][1] << "\n ";
		}
	}
//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		Real J = rho_0_[particle_i] / rho_n_[particle_i];
		Mat2d F = F_[particle_i];
		Mat2d stress = stress_[particle_i];
		Mat2d sigma = (F* stress*~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy
			+ 3.0 * sigmaxy * sigmaxy);
	}
//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\",\"v_x\", \"v_y\", \"x_norm\", \"y_norm\", \"von Mises stress\" \n";
		size_t number_of_particles = body_->number_of_particles_;

		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< i << "  "
				<< vel_n_[i][0] << " "
				<< vel_n_[i][1] << " "
				<< n_[i][0] << "  "
				<< n_[i][1] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
//=================================================================================================//
	void ActiveMuscleParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"Vx\", \"Vy\", \"von Mieses\", \" Ta \" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< i << "  "
				<< vel_n_[i][0] << " "
				<< vel_n_[i][1] << " "
				<< von_Mises_stress(i) << " "
				<< active_contraction_stress_[i] << "\n ";
		}
	}
//=================================================================================================//
}
