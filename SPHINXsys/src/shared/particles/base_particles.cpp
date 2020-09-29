/**
 * @file base_particles.cpp
 * @brief Definition of funcitons declared in base_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
 /**
  * @file 	base_particles.cpp
  * @author	Luhui Han, Chi ZHang and Xiangyu Hu
  * @version	0.1
  */

#include "base_particles.h"
#include "base_material.h"
#include "base_body.h"
#include "all_particle_generators.h"

namespace SPH
{
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody* body, BaseMaterial* base_material) : 
		base_material_(base_material), speed_max_(0.0), signal_speed_max_(0.0),
		real_particles_bound_(0), number_of_ghost_particles_(0),
		body_(body), body_name_(body->GetBodyName())
	{
		body->assignBaseParticle(this);
		base_material->assignBaseParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(pos_n_, registered_vectors_, vectors_map_, vectors_to_write_, "Position", false);
		registerAVariable(vel_n_, registered_vectors_, vectors_map_, vectors_to_write_, "Velocity", true);
		registerAVariable(dvel_dt_, registered_vectors_, vectors_map_, vectors_to_write_, "Acceleration", false);
		registerAVariable(dvel_dt_others_, registered_vectors_, vectors_map_, vectors_to_write_, "OtherAcceleration", false);

		registerAVariable(Vol_, registered_scalars_, scalars_map_, scalars_to_write_, "Volume", false);
		registerAVariable(Vol_0_, registered_scalars_, scalars_map_, scalars_to_write_, "InitialVolume", false);
		registerAVariable(sigma_0_, registered_scalars_, scalars_map_, scalars_to_write_, "InitialNumberDensity", false);
		registerAVariable(rho_n_, registered_scalars_, scalars_map_, scalars_to_write_, "Density", true);
		registerAVariable(rho_0_, registered_scalars_, scalars_map_, scalars_to_write_, "InitialDensity", false);
		registerAVariable(mass_, registered_scalars_, scalars_map_, scalars_to_write_, "Mass", false);
		registerAVariable(smoothing_length_, registered_scalars_, scalars_map_, scalars_to_write_, "SmoothingLength", false);

		ParticleGenerator* particle_generator = body_->particle_generator_;
		particle_generator->initialize(body_);
		particle_generator->CreateBaseParticles(this);
		delete particle_generator;

		real_particles_bound_ = body_->number_of_particles_;
	}
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody* body)
		: BaseParticles(body, new BaseMaterial()) 
	{

	}
	//=================================================================================================//
	void BaseParticles::initializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0)
	{
		particle_id_.push_back(pos_n_.size());
		is_sortable_.push_back(true);

		pos_n_.push_back(pnt);
		vel_n_.push_back(Vecd(0));
		dvel_dt_.push_back(Vecd(0));
		dvel_dt_others_.push_back(Vecd(0));

		Vol_0_.push_back(Vol_0);
		Vol_.push_back(Vol_0);
		sigma_0_.push_back(sigma_0);
		Real rho = base_material_->ReferenceDensity();
		rho_0_.push_back(rho);
		rho_n_.push_back(rho);
		mass_.push_back(rho * Vol_0);
		smoothing_length_.push_back(0);
	}
	//=================================================================================================//
	void BaseParticles::addABufferParticle()
	{
		particle_id_.push_back(pos_n_.size());
		is_sortable_.push_back(true);

		//update registered data in particle dynamics
		for (size_t i = 0; i != registered_matrices_.size(); ++i) 
			registered_matrices_[i]->push_back(Matd(0));
		for (size_t i = 0; i != registered_vectors_.size(); ++i) 
			registered_vectors_[i]->push_back(Vecd(0));
		for (size_t i = 0; i != registered_scalars_.size(); ++i) 
			registered_scalars_[i]->push_back(Real(0));
	}
	//=================================================================================================//
	void BaseParticles::copyFromAnotherParticle(size_t this_index, size_t another_index)
	{
		is_sortable_[this_index] = is_sortable_[another_index];

		updateFromAnotherParticle(this_index, another_index);
	}
	//=================================================================================================//
	void BaseParticles::updateFromAnotherParticle(size_t this_index, size_t another_index)
	{
		//update registered data in particle dynamics
		for (size_t i = 0; i != registered_matrices_.size(); ++i)
			(*registered_matrices_[i])[this_index] = (*registered_matrices_[i])[another_index];
		for (size_t i = 0; i != registered_vectors_.size(); ++i)
			(*registered_vectors_[i])[this_index] = (*registered_vectors_[i])[another_index];
		for (size_t i = 0; i != registered_scalars_.size(); ++i)
			(*registered_scalars_[i])[this_index] = (*registered_scalars_[i])[another_index];
	}
	//=================================================================================================//
	bool BaseParticles::isSwappingAllowed(size_t this_particle_index, size_t that_particle_index)
	{
		return  is_sortable_[this_particle_index] && is_sortable_[that_particle_index];
	}
	//=================================================================================================//
	size_t BaseParticles ::insertAGhostParticle(size_t index_i)
	{
		number_of_ghost_particles_ += 1;
		size_t expected_size = real_particles_bound_ + number_of_ghost_particles_;
		size_t expected_particle_index = expected_size - 1;
		if (expected_size <= pos_n_.size()) {
			copyFromAnotherParticle(expected_particle_index, index_i);
			particle_id_[expected_particle_index] = particle_id_[index_i];

		}
		else {
			addABufferParticle();
			copyFromAnotherParticle(expected_particle_index, index_i);
			particle_id_[expected_particle_index] = particle_id_[index_i];
		}
		return expected_particle_index;
	}
	//=================================================================================================//
	void  BaseParticles::mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction)
	{
		pos_n_[particle_index_i][axis_direction]
			= 2.0 * body_bound[axis_direction] - pos_n_[particle_index_i][axis_direction];
		vel_n_[particle_index_i][axis_direction] *= -1.0;
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;

		//write particle positions first
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		//write particles ID
		output_file << "    <DataArray Name=\"Particle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write vectors
		for (size_t l = 0; l != vectors_to_write_.size(); ++l) {
			string variable_name = vectors_to_write_[l];
			size_t variable_index_ = vectors_map_[variable_name];
			StdLargeVec<Vecd>& variable = *(registered_vectors_[variable_index_]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != number_of_particles; ++i) {
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
		
		//write scalars
		for (size_t l = 0; l != scalars_to_write_.size(); ++l) {
			string variable_name = scalars_to_write_[l];
			size_t variable_index_ = scalars_map_[variable_name];
			StdLargeVec<Real>& variable = *(registered_scalars_[variable_index_]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != number_of_particles; ++i) {
				output_file << fixed << setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		unique_ptr<XmlEngine> reload_xml(new XmlEngine(xml_name, ele_name));

		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			reload_xml->CreatXmlElement("particle");
			reload_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			reload_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			reload_xml->AddAttributeToElement<Real>("NumberDensity", sigma_0_[i]);
			reload_xml->AddElementToXmlDoc();
		}
		reload_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[number_of_particles] = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[number_of_particles] = volume;
			Real sigma = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "NumberDensity");
			sigma_0_[number_of_particles] = sigma;
			number_of_particles++;
		}

		if(number_of_particles != pos_n_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	BaseParticles* BaseParticles::pointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
}
