#include "solid_body_particles.h"

#include <iterator>

using namespace std;

namespace SPH {
	//===============================================================//
	void SolidBodyParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		output_file << "  <Piece Name =\"" <<  body_name_ << "\" NumberOfPoints=\"" << number_of_particles_ << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; 	i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
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

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << 0.0 << " ";
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
	void SolidBodyParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"x_norm\", \"y_norm\" \n";
		for (size_t i = 0; i != number_of_particles_; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
						<< base_particle_data_[i].pos_n_[1] << "  "
						<< base_particle_data_[i].particle_ID_ << "  "
						<< solid_body_data_[i].n_[0] << "  "
						<< solid_body_data_[i].n_[1] << "\n ";
		}
	}

	//===============================================================//
	void SolidBodyParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");
		
  		for(size_t i = 0; i != number_of_particles_; ++i)
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
	void SolidBodyParticles::InitialParticleFromRestartXmlFile(std::string &filefullpath)
	{
		/** Nothing should be done for non-moving BCs. */
	}

	Real ElasticBodyParticles::von_Mises_stress(size_t particle_i)
	{
		ElasticBodyParticleData &elastic_data_i
			= elastic_body_data_[particle_i];
		Real J = elastic_data_i.rho_0_ / elastic_data_i.rho_n_;
		Mat2d F = elastic_data_i.F_;
		Mat2d stress = elastic_data_i.stress_;
		Mat2d sigma = (F* stress*~F) / J;

		Real sigmaxx = sigma(0, 0);
		Real sigmayy = sigma(1, 1);
		Real sigmaxy = sigma(0, 1);

		return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy
			+ 3.0 * sigmaxy * sigmaxy);
	}
	//===============================================================//
	void ElasticBodyParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles_ << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
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
			output_file << elastic_body_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles_; ++i) {
			output_file << base_particle_data_[i].vel_n_[0] << " " << base_particle_data_[i].vel_n_[1] << " " << 0.0 << " ";
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
	void ElasticBodyParticles::WriteParticlesToPltFile(ofstream &output_file)
	{
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"x_norm\", \"y_norm\", \"von Mises stress\" \n";
		for (size_t i = 0; i != number_of_particles_; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< solid_body_data_[i].n_[0] << "  "
				<< solid_body_data_[i].n_[1] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
	//===============================================================//
	void ElasticBodyParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");
		
  		for(size_t i = 0; i != number_of_particles_; ++i)
  		{
  			restart_xml->CreatXmlElement("particle");
    		restart_xml->AddAttributeToElement("ID",base_particle_data_[i].particle_ID_);
    		restart_xml->AddAttributeToElement("Position",base_particle_data_[i].pos_n_);
    		restart_xml->AddAttributeToElement("Volume",base_particle_data_[i].Vol_);
    		restart_xml->AddAttributeToElement("Density", elastic_body_data_[i].rho_n_);
    		restart_xml->AddAttributeToElement("Velocity",base_particle_data_[i].vel_n_);
    		restart_xml->AddAttributeToElement("Displacement",elastic_body_data_[i].pos_temp_); 		
    		restart_xml->AddAttributeToElement("DefTensor", elastic_body_data_[i].F_);
    		restart_xml->AddElementToXmlDoc();
  		}
  		restart_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void ElasticBodyParticles::InitialParticleFromRestartXmlFile(std::string &filefullpath)
	{
		if (!fs::exists(filefullpath))
		{
			std::cout << "\n Error: the input file:"<< filefullpath << " is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}else{
			size_t number_of_particles = 0;
			XmlEngine* read_xml = new XmlEngine();
			read_xml->LoadXmlFile(filefullpath);
			SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
			for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
			{
				Vecd pos_ 			= read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Position");
				base_particle_data_[number_of_particles].pos_n_ = pos_;
				Real rst_Vol_ 		= read_xml->GetRequiredAttributeRealValue(ele_ite_, "Volume");
				base_particle_data_[number_of_particles].Vol_ = rst_Vol_;
				Vecd rst_vel_ 		= read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Velocity");
				base_particle_data_[number_of_particles].vel_n_ = rst_vel_;
				Real rst_rho_n_ 	= read_xml->GetRequiredAttributeRealValue(ele_ite_, "Density");
				elastic_body_data_[number_of_particles].rho_n_ = rst_rho_n_;	
				Vecd rst_dis_ 		= read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Displacement");
				elastic_body_data_[number_of_particles].pos_temp_ = rst_dis_;	
				Matd rst_F_ 	= read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
				elastic_body_data_[number_of_particles].F_ = rst_F_;		
				number_of_particles++;
			}
		}
	}
	//===============================================================//
}
