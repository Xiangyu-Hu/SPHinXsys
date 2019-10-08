#include "fluid_particles.h"
#include "base_body.h"

#include <iterator>

using namespace std;

namespace SPH
{
	//===============================================================//
	void FluidParticles::WriteParticlesToVtuFile(ofstream &output_file)
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
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].particle_ID_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].vel_n_[0] << " " << base_particle_data_[i].vel_n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Vorticity\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].vort_2d_ << " ";
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
	void FluidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"density\", \"u\", \"v\", \"Vorticity\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< fluid_particle_data_[i].rho_n_ << " "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " "
				<< fluid_particle_data_[i].vort_2d_ << " \n";
		}

	}
	//===============================================================//
	void FluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		XmlEngine* restart_xml = new XmlEngine(xml_name, ele_name);
		
		size_t number_of_particles = body_->number_of_particles_;
		for(size_t i = 0; i != number_of_particles; ++i)
  		{
  			restart_xml->CreatXmlElement("particle");
    		restart_xml->AddAttributeToElement("ID",base_particle_data_[i].particle_ID_);
    		restart_xml->AddAttributeToElement("Position",base_particle_data_[i].pos_n_);
    		restart_xml->AddAttributeToElement("Volume",base_particle_data_[i].Vol_);
    		restart_xml->AddAttributeToElement("Velocity",base_particle_data_[i].vel_n_);
    		restart_xml->AddAttributeToElement("Density",fluid_particle_data_[i].rho_n_);
    		restart_xml->AddElementToXmlDoc();
  		}
  		restart_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void FluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = pos_;
			Real rst_Vol_ = read_xml->GetRequiredAttributeRealValue(ele_ite_, "Volume");
			base_particle_data_[number_of_particles].Vol_ = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Velocity");
			base_particle_data_[number_of_particles].vel_n_ = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeRealValue(ele_ite_, "Density");
			fluid_particle_data_[number_of_particles].rho_n_ = rst_rho_n_;
			fluid_particle_data_[number_of_particles].vel_trans_(0);
			base_particle_data_[number_of_particles].dvel_dt_(0);
			number_of_particles++;
		}
	}
	//===============================================================//
	void ViscoelasticFluidParticles::WriteParticlesToVtuFile(ofstream &output_file)
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
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write data of particles
		output_file << "   <PointData  Vectors=\"vector\">\n";
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].particle_ID_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].vel_n_[0] << " " << base_particle_data_[i].vel_n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Vorticity\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].vort_2d_ << " ";
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
	void ViscoelasticFluidParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"density\", \"u\", \"v\" \n";

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< fluid_particle_data_[i].rho_n_ << " "
				<< base_particle_data_[i].vel_n_[0] << " "
				<< base_particle_data_[i].vel_n_[1] << " \n";
		}
	}
	//===============================================================//
	void ViscoelasticFluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");
		size_t number_of_particles = base_particle_data_.size();
		for(size_t i = 0; i != number_of_particles; ++i)
  		{
  			restart_xml->CreatXmlElement("particle");
    		restart_xml->AddAttributeToElement("ID",base_particle_data_[i].particle_ID_);
    		restart_xml->AddAttributeToElement("Position",base_particle_data_[i].pos_n_);
    		restart_xml->AddAttributeToElement("Volume",base_particle_data_[i].Vol_);
    		restart_xml->AddElementToXmlDoc();
  		}
  		restart_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void ViscoelasticFluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function ViscoelasticFluidParticles::ReadParticleFromXmlForRestart is not done. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
}