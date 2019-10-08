#include "solid_particles.h"
#include "base_body.h"

using namespace std;

namespace SPH {
//=================================================================================================//
	void SolidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		output_file << "  <Piece Name =\"" <<  body_name_ << "\" NumberOfPoints=\"" << number_of_particles << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; 	i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
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

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << 0.0 << " ";
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
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"x_norm\", \"y_norm\" \n";
		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
						<< base_particle_data_[i].pos_n_[1] << "  "
						<< base_particle_data_[i].particle_ID_ << "  "
						<< solid_body_data_[i].n_[0] << "  "
						<< solid_body_data_[i].n_[1] << "\n ";
		}
	}
//=================================================================================================//
	void SolidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");

		size_t number_of_particles = body_->number_of_particles_;
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
//=================================================================================================//
	void SolidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		/** Nothing should be done for non-moving BCs. */
	}
//=================================================================================================//
	Real ElasticSolidParticles::von_Mises_stress(size_t particle_i)
	{
		ElasticSolidParticleData &elastic_data_i
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
//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;
		output_file << "  <Piece Name =\"" << body_name_ << "\" NumberOfPoints=\"" << number_of_particles << "\" NumberOfCells=\"0\">\n";

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << base_particle_data_[i].pos_n_[0] << " " << base_particle_data_[i].pos_n_[1] << " " << 0.0 << " ";
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
			output_file << fixed << setprecision(9) << elastic_body_data_[i].rho_n_ << " ";
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
			output_file << fixed << setprecision(9) << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << 0.0 << " ";
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
		output_file << " VARIABLES = \" x \", \"y\", \"ID\", \"x_norm\", \"y_norm\", \"von Mises stress\" \n";
		size_t number_of_particles = body_->number_of_particles_;

		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< solid_body_data_[i].n_[0] << "  "
				<< solid_body_data_[i].n_[1] << "  "
				<< von_Mises_stress(i) << "\n ";
		}
	}
//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");
		
		size_t number_of_particles = body_->number_of_particles_;
		for(size_t i = 0; i != number_of_particles; ++i)
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
//=================================================================================================//
	void ElasticSolidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
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
			elastic_body_data_[number_of_particles].rho_n_ = rst_rho_n_;
			Vecd rst_dis_ = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Displacement");
			elastic_body_data_[number_of_particles].pos_temp_ = rst_dis_;
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			elastic_body_data_[number_of_particles].F_ = rst_F_;
			number_of_particles++;
		}
	}
//=================================================================================================//
	void MuscleParticles::WriteParticlesToVtuFile(ofstream &output_file)
	{
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

		output_file << "    <DataArray Name=\"Volume \" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << base_particle_data_[i].Vol_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Voltage \" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << muscle_body_data_[i].voltage_n_ << " ";
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
		/*
		output_file << " VARIABLES = \" x \", \"y\", \"ID\",\"Volume\", \" Voltage \" \n";
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			output_file << base_particle_data_[i].pos_n_[0] << "  "
				<< base_particle_data_[i].pos_n_[1] << "  "
				<< base_particle_data_[i].particle_ID_ << "  "
				<< base_particle_data_[i].Vol_ << "  "
				<< muscle_body_data_[i].voltage_n_ << "\n ";
		}
		*/
		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "Phi, " << "\n";
		output_file << "zone i=" << 50 << "  j=" << 50 << "  k=" << 1
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";
	
		for (size_t j = 0; j != 50; ++j)
		{
			for (size_t i = 0; i != 50; ++i)
			{
				int index_id = j * 50 + i;
				output_file << base_particle_data_[index_id].pos_n_[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != 50; ++j)
		{
			for (size_t i = 0; i != 50; ++i)
			{
				int index_id = j * 50 + i;
				output_file << base_particle_data_[index_id].pos_n_[1] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != 50; ++j)
		{
			for (size_t i = 0; i != 50; ++i)
			{
				int index_id = j * 50 + i;
				output_file << muscle_body_data_[index_id].voltage_n_ << " ";

			}
			output_file << " \n";
		}
	}
//=================================================================================================//
}
