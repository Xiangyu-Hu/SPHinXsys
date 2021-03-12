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
		total_real_particles_(0), 
		real_particles_bound_(0), total_ghost_particles_(0),
		body_(body), body_name_(body->getBodyName())
	{
		body->assignBaseParticles(this);
		base_material->assignBaseParticles(this);
		rho_0_ = base_material->ReferenceDensity();
		sigma_0_ = body->computeReferenceNumberDensity(Vecd(0));
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(pos_n_, "Position");
		registerAVariable<indexVector, Vecd>(vel_n_, "Velocity", true);
		registerAVariable<indexVector, Vecd>(dvel_dt_, "Acceleration");
		registerAVariable<indexVector, Vecd>(dvel_dt_others_, "OtherAcceleration");

		registerAVariable<indexScalar, Real>(Vol_, "Volume");
		registerAVariable<indexScalar, Real>(rho_n_, "Density", true);
		registerAVariable<indexScalar, Real>(mass_, "Mass");
		registerAVariable<indexScalar, Real>(h_ratio_, "SmoothingLengthRatio");

		ParticleGenerator* particle_generator = body_->particle_generator_;
		particle_generator->initialize(body_);
		particle_generator->createBaseParticles(this);

		real_particles_bound_ = total_real_particles_;
	}
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody* body)
		: BaseParticles(body, new BaseMaterial()) {}
	//=================================================================================================//
	void BaseParticles::initializeABaseParticle(Vecd pnt, Real Vol_0)
	{
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		pos_n_.push_back(pnt);
		vel_n_.push_back(Vecd(0));
		dvel_dt_.push_back(Vecd(0));
		dvel_dt_others_.push_back(Vecd(0));

		Vol_.push_back(Vol_0);
		rho_n_.push_back(rho_0_);
		mass_.push_back(rho_0_ * Vol_0);
		h_ratio_.push_back(1.0);
	}
	//=================================================================================================//
	void BaseParticles::addABufferParticle()
	{
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		//update registered data in particle dynamics
		for (size_t i = 0; i != std::get<indexScalar>(all_particle_data_).size(); ++i)
			std::get<indexScalar>(all_particle_data_)[i]->push_back(Real(0));
		for (size_t i = 0; i != std::get<indexVector>(all_particle_data_).size(); ++i)
			std::get<indexVector>(all_particle_data_)[i]->push_back(Vecd(0));
		for (size_t i = 0; i != std::get<indexMatrix>(all_particle_data_).size(); ++i)
			std::get<indexMatrix>(all_particle_data_)[i]->push_back(Matd(0));
		for (size_t i = 0; i != std::get<indexInteger>(all_particle_data_).size(); ++i)
			std::get<indexInteger>(all_particle_data_)[i]->push_back(int(0));
		for (size_t i = 0; i != std::get<indexBoolean>(all_particle_data_).size(); ++i)
			std::get<indexBoolean>(all_particle_data_)[i]->push_back(bool(0));
	}
	//=================================================================================================//
	void BaseParticles::copyFromAnotherParticle(size_t this_index, size_t another_index)
	{
		updateFromAnotherParticle(this_index, another_index);
	}
	//=================================================================================================//
	void BaseParticles::updateFromAnotherParticle(size_t this_index, size_t another_index)
	{
		//update registered data in particle dynamics
		for (size_t i = 0; i != std::get<0>(all_particle_data_).size(); ++i)
			(*std::get<0>(all_particle_data_)[i])[this_index] = (*std::get<0>(all_particle_data_)[i])[another_index];
		for (size_t i = 0; i != std::get<1>(all_particle_data_).size(); ++i)
			(*std::get<1>(all_particle_data_)[i])[this_index] = (*std::get<1>(all_particle_data_)[i])[another_index];
		for (size_t i = 0; i != std::get<2>(all_particle_data_).size(); ++i)
			(*std::get<2>(all_particle_data_)[i])[this_index] = (*std::get<2>(all_particle_data_)[i])[another_index];
		for (size_t i = 0; i != std::get<3>(all_particle_data_).size(); ++i)
			(*std::get<3>(all_particle_data_)[i])[this_index] = (*std::get<3>(all_particle_data_)[i])[another_index];
		for (size_t i = 0; i != std::get<4>(all_particle_data_).size(); ++i)
			(*std::get<4>(all_particle_data_)[i])[this_index] = (*std::get<4>(all_particle_data_)[i])[another_index];
	}
	//=================================================================================================//
	size_t BaseParticles ::insertAGhostParticle(size_t index_i)
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
			addABufferParticle();
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
	void BaseParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		size_t total_real_particles = total_real_particles_;

		//write particle positions first
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

		//write vectors
		for (size_t l = 0; l != variables_to_write_[indexVector].size(); ++l) {
			string variable_name = variables_to_write_[indexVector][l];
			size_t variable_index = name_index_maps_[indexVector][variable_name];
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[variable_index]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << fixed << setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
		
		//write scalars
		for (size_t l = 0; l != variables_to_write_[indexScalar].size(); ++l) {
			string variable_name = variables_to_write_[indexScalar][l];
			size_t variable_index_ = name_index_maps_[indexScalar][variable_name];
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[variable_index_]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << fixed << setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write booleans
		for (size_t l = 0; l != variables_to_write_[indexBoolean].size(); ++l) {
			string variable_name = variables_to_write_[indexBoolean][l];
			size_t variable_index_ = name_index_maps_[indexBoolean][variable_name];
			StdLargeVec<bool>& variable = *(std::get<indexBoolean>(all_particle_data_)[variable_index_]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int8\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToPltFile(ofstream& output_file)
	{
		size_t total_real_particles = total_real_particles_;
		output_file << " VARIABLES = \"x\",\"y\",\"z\",\"ID\"";
		for (size_t l = 0; l != variables_to_write_[indexVector].size(); ++l) {
			string variable_name = variables_to_write_[indexVector][l];
			output_file << ",\"" << variable_name << "_x\"" << ",\"" << variable_name << "_y\"" << ",\"" << variable_name << "_z\"";
		};
		for (size_t l = 0; l != variables_to_write_[indexScalar].size(); ++l) {
			string variable_name = variables_to_write_[indexScalar][l];
			output_file << ",\"" << variable_name << "\"";
		};
		output_file << "\n";

		//write particle positions first
		for (size_t i = 0; i != total_real_particles; ++i) {
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
				<< i << " ";


			for (size_t l = 0; l != variables_to_write_[indexVector].size(); ++l) {
				string variable_name = variables_to_write_[indexVector][l];
				size_t variable_index = name_index_maps_[indexVector][variable_name];
				StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[variable_index]);
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			};

			for (size_t l = 0; l != variables_to_write_[indexScalar].size(); ++l) {
				string variable_name = variables_to_write_[indexScalar][l];
				size_t variable_index = name_index_maps_[indexScalar][variable_name];
				StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[variable_index]);
				output_file << variable[i] << " ";
			};
			output_file << "\n";
		};
	}
	//=================================================================================================//
	void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		unique_ptr<XmlEngine> reload_xml(new XmlEngine(xml_name, ele_name));

		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			reload_xml->creatXmlElement("particle");
			reload_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			reload_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			reload_xml->AddElementToXmlDoc();
		}
		reload_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[total_real_particles] = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[total_real_particles] = volume;
			total_real_particles++;
		}

		if(total_real_particles != pos_n_.size())
		{
			std::cout << "\n Error: reload particle number does not match!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
}
