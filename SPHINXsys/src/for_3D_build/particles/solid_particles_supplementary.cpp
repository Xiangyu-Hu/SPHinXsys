#include "solid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;
//=================================================================================================//
namespace SPH {
//=================================================================================================//
	void SolidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		//define output iterators
		std::vector<Real>::const_iterator scalar_iterater;
		std::vector<Vecd>::const_iterator vector_iterater;
		std::vector<int>::const_iterator int_iterater;

		size_t number_of_particles = body_->number_of_particles_;
		output_file << "  <Piece Name =\"" <<  body_name_ << "\" NumberOfPoints=\"" << number_of_particles << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << base_particle_data_[i].pos_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << solid_body_data_[i].n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </PointData>\n";

		//write empty cells
		output_file << "   <Cells>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "   </Cells>\n";

		output_file << "  </Piece>\n";
	}
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
	void ElasticSolidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		//define output iterators
		std::vector<Real>::const_iterator scalar_iterater;
		std::vector<Vecd>::const_iterator vector_iterater;
		std::vector<int>::const_iterator int_iterater;

		size_t number_of_particles = body_->number_of_particles_;
		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << base_particle_data_[i].pos_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";

		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << elastic_body_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) <<  base_particle_data_[i].vel_n_[0] << " " <<  base_particle_data_[i].vel_n_[1] << " " <<  base_particle_data_[i].vel_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << solid_body_data_[i].n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </PointData>\n";


		//write empty cells
		output_file << "   <Cells>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "   </Cells>\n";

		output_file << "  </Piece>\n";
	}
//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"u\", \"v\", \"w\", \"ID\", \"x_norm\", \"y_norm\", \"z_norm\" \n";

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
				<< solid_body_data_[i].n_[2] << "\n ";
		}
	}
//=================================================================================================//
	void MuscleParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		//define output iterators
		std::vector<Real>::const_iterator scalar_iterater;
		std::vector<Vecd>::const_iterator vector_iterater;
		std::vector<int>::const_iterator int_iterater;
		size_t number_of_particles = body_->number_of_particles_;

		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << base_particle_data_[i].pos_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";

		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Voltage\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << muscle_body_data_[i].voltage_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"T_a\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << muscle_body_data_[i].T_a_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "   </PointData>\n";

		//write empty cells
		output_file << "   <Cells>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
		output_file << "    </DataArray>\n";
		output_file << "   </Cells>\n";

		output_file << "  </Piece>\n";
	}
//=================================================================================================//
	void MuscleParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"Vx\", \"Vy\", \"Vz\", \"Voltage\",\"Ta\" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< i << "  "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< base_particle_data_[i].vel_n_[2] << " "
				<< muscle_body_data_[i].voltage_n_ << "  "
				<< muscle_body_data_[i].T_a_ << "\n ";
		}
	}
//=================================================================================================//
}
