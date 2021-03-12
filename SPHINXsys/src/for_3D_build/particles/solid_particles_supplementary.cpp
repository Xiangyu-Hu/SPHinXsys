#include "solid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;
//=================================================================================================//
namespace SPH {
	//=============================================================================================//
	void SolidParticles::ParticleTranslationAndRotation(Transformd& transform) 
	{
			std::cout << "\n Error: the function ParticleTranslationAndRotation in 3d is not defined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
	}
	//=================================================================================================//
	void SolidParticles::writeParticlesToPltFile(ofstream& output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\" \n";

		for (size_t i = 0; i != total_real_particles_; ++i)
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
		Real J = rho_0_ / rho_n_[particle_i];
		Mat3d F = F_[particle_i];
		Mat3d stress = stress_PK1_[particle_i];
		Mat3d sigma = (stress * ~F) / J;

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
	void ElasticSolidParticles::writeParticlesToPltFile(ofstream& output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"u\", \"v\", \"w\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\", \"von Mises\" \n";

		for (size_t i = 0; i != total_real_particles_; ++i)
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
	void ActiveMuscleParticles::writeParticlesToPltFile(ofstream& output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"Vx\", \"Vy\", \"Vz\", \"Ta\" ,\"von Mises \" \n";
		for (size_t i = 0; i != total_real_particles_; ++i)
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
				<< von_Mises_stress(i) << "\n ";
		}
	}
	//=================================================================================================//
	void ShellParticles::writeParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"u\", \"v\", \"w\", \"x_angle\", \"y_angle\", \"x_angle_v\", \"y_angle_v\", \"ID\", \"x_pseudo_norm\", \"y_pseudo_norm\", \"z_pseudo_norm\", \"von Mises\" \n";

		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			output_file
				<< pos_n_[i][0] << "  "
				<< pos_n_[i][1] << "  "
				<< pos_n_[i][2] << "  "
				<< vel_n_[i][0] << "  "
				<< vel_n_[i][1] << "  "
				<< vel_n_[i][2] << "  "
				<< rotation_[i][0] << "  "
				<< rotation_[i][1] << "  "
				<< angular_vel_[i][0] << "  "
				<< angular_vel_[i][1] << "  "
				<< i << "  "
				<< pseudo_n_[i][0] << "  "
				<< pseudo_n_[i][1] << "  "
				<< pseudo_n_[i][2] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
	//=================================================================================================//
}
