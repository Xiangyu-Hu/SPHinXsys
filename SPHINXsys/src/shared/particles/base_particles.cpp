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
	//===============================================================//
	size_t BaseParticleData::total_number_of_particles_ = 0;
	//===============================================================//
	BaseParticleData::BaseParticleData()
		: vel_n_(0), pos_n_(0), Vol_(0), Vol_0_(0), sigma_0_(0),
		dvel_dt_others_(0), dvel_dt_(0), is_sortable_(true),
		particle_id_(0)
	{
		total_number_of_particles_++;
	}
	//===============================================================//
	BaseParticleData::BaseParticleData(Vecd position, Real Vol_0, Real sigma_0) 
		: vel_n_(0), pos_n_(position), Vol_(Vol_0), Vol_0_(Vol_0), sigma_0_(sigma_0),
		dvel_dt_others_(0), dvel_dt_(0), is_sortable_(true), particle_id_(0)
	{
		total_number_of_particles_++;
	}
	//===============================================================//
	BaseParticles::BaseParticles(SPHBody* body, BaseMaterial* base_material)
		: body_(body), base_material_(base_material),
		body_name_(body->GetBodyName()), speed_max_(0.0)
	{
		body->base_particles_ = this;
		base_material->AssignParticles(this);

		switch (body->particle_generator_op_)
		{
		case ParticlesGeneratorOps::lattice: {
			particle_generator_ = new ParticleGeneratorLattice(*body, body->body_region_);
			break;
		}

		case ParticlesGeneratorOps::direct: {
			particle_generator_ = new ParticleGeneratorDirect(*body);
			break;
		}

		default: {
			std::cout << "\n FAILURE: the type of particle generator is undefined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
			break;
		}
		}

		particle_generator_->CreateBaseParticles();
		real_particles_bound_ = body_->number_of_particles_;
		number_of_ghost_particles_ = 0;
	}
	//===============================================================//
	void BaseParticles::InitializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0)
	{
		size_t particle_index = base_particle_data_.size();
		base_particle_data_.push_back(BaseParticleData(pnt, Vol_0, sigma_0));
		base_particle_data_[particle_index].particle_id_ = particle_index;
	}
	//===============================================================//
	void BaseParticles::AddABufferParticle()
	{
		size_t particle_index = base_particle_data_.size();
		base_particle_data_.push_back(BaseParticleData());
		base_particle_data_[particle_index].particle_id_ = particle_index;
	}
	//===============================================================//
	void BaseParticles::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		size_t particle_id = base_particle_data_[this_particle_index].particle_id_;
		base_particle_data_[this_particle_index] = base_particle_data_[another_particle_index];
		base_particle_data_[this_particle_index].particle_id_ = particle_id;
	}
	//===============================================================//
	void BaseParticles::UpdateFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		base_particle_data_[this_particle_index].pos_n_ = base_particle_data_[another_particle_index].pos_n_;
		base_particle_data_[this_particle_index].vel_n_ = base_particle_data_[another_particle_index].vel_n_;
	}
	//===============================================================//
	void BaseParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		std::swap(base_particle_data_[this_particle_index], base_particle_data_[that_particle_index]);
	}
	//===============================================================//
	bool BaseParticles::allowSwapping(size_t this_particle_index, size_t that_particle_index)
	{
		return  base_particle_data_[this_particle_index].is_sortable_
			    && base_particle_data_[that_particle_index].is_sortable_;
	}
	//=================================================================================================//
	size_t BaseParticles ::insertAGhostParticle(size_t index_particle_i)
	{
		number_of_ghost_particles_ += 1;
		size_t expected_size = real_particles_bound_ + number_of_ghost_particles_;
		size_t expected_particle_index = expected_size - 1;
		if (expected_size <= base_particle_data_.size()) {
			CopyFromAnotherParticle(expected_particle_index, index_particle_i);
			base_particle_data_[expected_particle_index].particle_id_
				= base_particle_data_[index_particle_i].particle_id_;

		}
		else {
			AddABufferParticle();
			CopyFromAnotherParticle(expected_particle_index, index_particle_i);
			base_particle_data_[expected_particle_index].particle_id_
				= base_particle_data_[index_particle_i].particle_id_;
		}
		return expected_particle_index;
	}
	//===============================================================//
	void BaseParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		size_t number_of_particles = body_->number_of_particles_;

		//write coordinates of particles
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			Vec3d particle_position = upgradeToVector3D(base_particle_data_[i].pos_n_);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
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

		output_file << "    <DataArray Name=\"Number Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << base_particle_data_[i].sigma_0_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			Vec3d particle_velocity = upgradeToVector3D(base_particle_data_[i].vel_n_);
			output_file << particle_velocity[0] << " " << particle_velocity[1] << " " << particle_velocity[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//===============================================================//
	void BaseParticles::WriteToXmlForReloadParticle(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		XmlEngine* reload_xml = new XmlEngine(xml_name, ele_name);

		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			reload_xml->CreatXmlElement("particle");
			reload_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			reload_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			reload_xml->AddAttributeToElement<Real>("NumberDensity", base_particle_data_[i].sigma_0_);
			reload_xml->AddElementToXmlDoc();
		}
		reload_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void BaseParticles::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			base_particle_data_[number_of_particles].Vol_ = volume;
			Real sigma = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "NumberDensity");
			base_particle_data_[number_of_particles].sigma_0_ = sigma;
			number_of_particles++;
		}

		if(number_of_particles != base_particle_data_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//===============================================================//
	BaseParticles* BaseParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
	void  BaseParticles
		::mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction)
	{
		BaseParticleData & base_particle_data_i = base_particle_data_[particle_index_i];
		base_particle_data_i.pos_n_[axis_direction]
			= 2.0 * body_bound[axis_direction] - base_particle_data_i.pos_n_[axis_direction];
		base_particle_data_i.vel_n_[axis_direction] *= -1.0;
	}
	//===============================================================//

}