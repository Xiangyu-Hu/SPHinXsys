#include "weakly_compressible_fluid_particles.h"

#include <iterator>

using namespace std;

namespace SPH
{
	//===============================================================//
	void WeaklyCompressibleFluidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		//define output iterators
		std::vector<Real>::const_iterator scalar_iterater;
		std::vector<Vecd>::const_iterator vector_iterater;
		std::vector<int>::const_iterator int_iterater;

		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles_ << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << base_particle_data_[i].pos_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].particle_ID_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << fluid_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].vel_n_[0] << " " << base_particle_data_[i].vel_n_[1] << " " << base_particle_data_[i].vel_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Vorticity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << fluid_data_[i].vorticity_[0] << " " << fluid_data_[i].vorticity_[1] << " " << fluid_data_[i].vorticity_[2] << " ";
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
	//===============================================================//
	void WeaklyCompressibleFluidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"density\", \"u\", \"v\",\"w\" \n";
		for (size_t i = 0; i != number_of_particles_; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< fluid_data_[i].rho_n_ << " "
				<< base_particle_data_[i].vel_n_[0] << " " 
				<< base_particle_data_[i].vel_n_[1] << " " 
				<< base_particle_data_[i].vel_n_[2] << "\n ";
		}
	}
	//===============================================================//
	void WeaklyCompressibleFluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
	void WeaklyCompressibleFluidParticles::InitialParticleFromRestartXmlFile(std::string &filefullpath)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
	void Oldroyd_B_FluidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		//define output iterators
		std::vector<Real>::const_iterator scalar_iterater;
		std::vector<Vecd>::const_iterator vector_iterater;
		std::vector<int>::const_iterator int_iterater;

		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles_ << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << base_particle_data_[i].pos_n_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].particle_ID_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << fluid_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].vel_n_[0] << " " << base_particle_data_[i].vel_n_[1] << " " << base_particle_data_[i].vel_n_[2] << " ";
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
	//===============================================================//
	void Oldroyd_B_FluidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\",\"z\", \"ID\", \"density\", \"u\", \"v\",\"w\" \n";
		for (size_t i = 0; i != number_of_particles_; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].pos_n_[2] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< fluid_data_[i].rho_n_ << " "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< base_particle_data_[i].vel_n_[2] << "\n ";
		}
	}
	//===============================================================//
	void Oldroyd_B_FluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
	void Oldroyd_B_FluidParticles::InitialParticleFromRestartXmlFile(std::string &filefullpath)
	{
		cout << "\n This function is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
}