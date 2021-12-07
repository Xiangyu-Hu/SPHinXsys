/**
 * @file base_particles.cpp
 * @brief Definition of functions declared in base_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "base_particles.h"

#include "base_body.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody* body, BaseMaterial* base_material) :
		base_material_(base_material), speed_max_(0.0), signal_speed_max_(0.0),
		total_real_particles_(0), real_particles_bound_(0), total_ghost_particles_(0),
		body_(body), body_name_(body->getBodyName()),
		restart_xml_engine_("xml_restart", "particles"),
		reload_xml_engine_("xml_particle_reload", "particles")
	{
		body->assignBaseParticles(this);
		rho0_ = base_material->ReferenceDensity();
		sigma0_ = body->particle_adaptation_->ReferenceNumberDensity();
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(pos_n_, "Position");
		registerAVariable<indexVector, Vecd>(vel_n_, "Velocity");
		registerAVariable<indexVector, Vecd>(dvel_dt_, "Acceleration");
		registerAVariable<indexVector, Vecd>(dvel_dt_prior_, "OtherAcceleration");
		registerAVariable<indexScalar, Real>(Vol_, "Volume");
		registerAVariable<indexScalar, Real>(rho_n_, "Density");
		registerAVariable<indexScalar, Real>(mass_, "Mass");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<indexVector, Vecd>("Velocity");
		addAVariableToWrite<indexVector, Vecd>("Acceleration");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Position");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Velocity");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Acceleration");
		addAVariableNameToList<indexScalar, Real>(variables_to_restart_, "Volume");

		ParticleGenerator* particle_generator = body_->particle_generator_;
		particle_generator->initialize(body_);
		particle_generator->createBaseParticles(this);
		real_particles_bound_ = total_real_particles_;

		body->particle_adaptation_->assignBaseParticles(this);
		base_material->assignBaseParticles(this);
	}
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody* body)
		: BaseParticles(body, new BaseMaterial()) {}
	//=================================================================================================//
	void BaseParticles::initializeABaseParticle(Vecd pnt, Real Vol_0)
	{
		total_real_particles_++;
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		pos_n_.push_back(pnt);
		vel_n_.push_back(Vecd(0));
		dvel_dt_.push_back(Vecd(0));
		dvel_dt_prior_.push_back(Vecd(0));

		Vol_.push_back(Vol_0);
		rho_n_.push_back(rho0_);
		mass_.push_back(rho0_ * Vol_0);
	}
	//=================================================================================================//
	void BaseParticles::addAParticleEntry()
	{
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		loopParticleData<addAParticleDataValue>(all_particle_data_);
	}
	//=================================================================================================//
	void BaseParticles::addBufferParticles(size_t buffer_size)
	{
		for (size_t i = 0; i != buffer_size; ++i)
		{
			addAParticleEntry();
		}
		real_particles_bound_ += buffer_size;
	}
	//=================================================================================================//
	void BaseParticles::copyFromAnotherParticle(size_t this_index, size_t another_index)
	{
		updateFromAnotherParticle(this_index, another_index);
	}
	//=================================================================================================//
	void BaseParticles::updateFromAnotherParticle(size_t this_index, size_t another_index)
	{
		loopParticleData<copyAParticleDataValue>(all_particle_data_, this_index, another_index);
	}
	//=================================================================================================//
	size_t BaseParticles::insertAGhostParticle(size_t index_i)
	{
		total_ghost_particles_ += 1;
		size_t expected_size = real_particles_bound_ + total_ghost_particles_;
		size_t expected_particle_index = expected_size - 1;
		if (expected_size <= pos_n_.size()) {
			copyFromAnotherParticle(expected_particle_index, index_i);
			/** For a ghost particle, its sorted id is that of corresponding real particle. */
			sorted_id_[expected_particle_index] = index_i;

		}
		else {
			addAParticleEntry();
			copyFromAnotherParticle(expected_particle_index, index_i);
			/** For a ghost particle, its sorted id is that of corresponding real particle. */
			sorted_id_[expected_particle_index] = index_i;
		}
		return expected_particle_index;
	}
	//=================================================================================================//
	void BaseParticles::switchToBufferParticle(size_t index_i)
	{
		size_t last_real_particle_index = total_real_particles_ - 1;
		updateFromAnotherParticle(index_i, last_real_particle_index);
		unsorted_id_[index_i] = unsorted_id_[last_real_particle_index];
		total_real_particles_ -= 1;
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToVtuFile(std::ostream& output_file)
	{
		size_t total_real_particles = total_real_particles_;

		//write current/final particle positions first
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		//write sorted particles ID
		output_file << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write unsorted particles ID
		output_file << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << unsorted_id_[i] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write matrices
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd>& variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				Mat3d matrix_value = upgradeToMatrix3D(variable[i]);
				for (int k = 0; k != 3; ++k) {
					Vec3d col_vector = matrix_value.col(k);
					output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
				}
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write vectors
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write scalars
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write integers
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int>& variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writeSurfaceParticlesToVtuFile(std::ostream& output_file, ShapeSurface& surface_particles)
	{
		size_t total_surface_particles = surface_particles.body_part_particles_.size();

		//write current/final particle positions first
		// precision: 3 - 0.1 mm accuracy
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vec3d particle_position = upgradeToVector3D(pos_n_[particle_i]);
			output_file << std::fixed << std::setprecision(1) << particle_position[0] << " "
						<< std::fixed << std::setprecision(1) << particle_position[1] << " "
						<< std::fixed << std::setprecision(1) << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		/* // REMOVED FOR PRODUCTION
		//write matrices
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd>& variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				Mat3d matrix_value = upgradeToMatrix3D(variable[particle_i]);
				for (int k = 0; k != 3; ++k) {
					Vec3d col_vector = matrix_value.col(k);
					output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
				}
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
		*/

		//write vectors
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				Vec3d vector_value = upgradeToVector3D(variable[particle_i]);
				output_file << std::fixed << std::setprecision(2) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
		/*
		//write scalars
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				output_file << std::fixed << std::setprecision(9) << variable[particle_i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write integers
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int>& variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				output_file << std::fixed << std::setprecision(9) << variable[particle_i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
		*/
	}
	//=================================================================================================//
	void BaseParticles::writePltFileHeader(std::ofstream& output_file)
	{
		output_file << " VARIABLES = \"x\",\"y\",\"z\",\"ID\"";
		
		for (size_t l = 0; l != variables_to_write_[indexInteger].size(); ++l) 
		{
			std::string variable_name = variables_to_write_[indexInteger][l].first;
			output_file << ",\"" << variable_name << "\"";
		};


		for (size_t l = 0; l != variables_to_write_[indexVector].size(); ++l) {
			std::string variable_name = variables_to_write_[indexVector][l].first;
			output_file << ",\"" << variable_name << "_x\"" << ",\"" << variable_name << "_y\"" << ",\"" << variable_name << "_z\"";
		};
		for (size_t l = 0; l != variables_to_write_[indexScalar].size(); ++l) {
			std::string variable_name = variables_to_write_[indexScalar][l].first;
			output_file << ",\"" << variable_name << "\"";
		};
	}
	//=================================================================================================//
	void BaseParticles::writePltFileParticleData(std::ofstream& output_file, size_t index_i)
	{
		//write particle positions and index first
		Vec3d particle_position = upgradeToVector3D(pos_n_[index_i]);
		output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
			<< index_i << " ";
			
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int>& variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << variable[index_i] << " ";
		};

		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			Vec3d vector_value = upgradeToVector3D(variable[index_i]);
			output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
		};

		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << variable[index_i] << " ";
		};
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToPltFile(std::ofstream& output_file)
	{
		writePltFileHeader(output_file);
		output_file << "\n";

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i) {
			writePltFileParticleData(output_file, i);
			output_file << "\n";
		};
	}
	//=================================================================================================//
	void BaseParticles::resizeXmlDocForParticles(XmlEngine& xml_engine)
	{
		size_t total_elements = xml_engine.SizeOfXmlDoc();

		if (total_elements <= total_real_particles_)
		{
			for (size_t i = total_elements; i != total_real_particles_; ++i)
				xml_engine.addElementToXmlDoc("particle");
		}
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToXmlForRestart(std::string& filefullpath)
	{
		resizeXmlDocForParticles(restart_xml_engine_);
		WriteAParticleVariableToXml write_variable_to_xml(restart_xml_engine_, total_real_particles_);
		loopParticleData<loopVariabaleNameList>(all_particle_data_, variables_to_restart_, write_variable_to_xml);
		restart_xml_engine_.writeToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readParticleFromXmlForRestart(std::string& filefullpath)
	{
		restart_xml_engine_.loadXmlFile(filefullpath);
		ReadAParticleVariableFromXml read_variable_from_xml(restart_xml_engine_, total_real_particles_);
		loopParticleData<loopVariabaleNameList>(all_particle_data_, variables_to_restart_, read_variable_from_xml);
	}
	//=================================================================================================//
	void BaseParticles::writeToXmlForReloadParticle(std::string& filefullpath)
	{
		resizeXmlDocForParticles(reload_xml_engine_);
		SimTK::Xml::element_iterator ele_ite = reload_xml_engine_.root_element_.element_begin();
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			reload_xml_engine_.setAttributeToElement(ele_ite, "Position", pos_n_[i]);
			reload_xml_engine_.setAttributeToElement(ele_ite, "Volume", Vol_[i]);
			ele_ite++;
		}
		reload_xml_engine_.writeToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readFromXmlForReloadParticle(std::string& filefullpath)
	{
		reload_xml_engine_.loadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite = reload_xml_engine_.root_element_.element_begin();
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			reload_xml_engine_.getRequiredAttributeValue<Vecd>(ele_ite, "Position", pos_n_[i]);
			reload_xml_engine_.getRequiredAttributeValue<Real>(ele_ite, "Volume", Vol_[i]);
			ele_ite++;
		}

		if (reload_xml_engine_.SizeOfXmlDoc() != total_real_particles_)
		{
			std::cout << "\n Error: reload particle number does not match!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
}
