/**
 * @file 	base_body.cpp
 * @brief 	This is the base classes of SPH bodies. The real body is for 
 *			that with cell linked list and the fivtitious one doesnot.     
 * 			Before the defination of the SPH bodies, the shapes with complex 
 *			geometries, i.e. those are produced by adavnced binary operation, 
 * 			such as intersection, should be produced first.
 * 			Then, all shapes used in body definition should be either contain 
 * 			or not contain each other. 
 *			Partial overlap between them are not premitted.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "base_body.h"
#include "sph_system.h"
#include "base_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"
#include "particle_generator_lattice.h"

namespace SPH
{
	//===========================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, string body_name,
		Particles &base_particles, int refinement_level, ParticlesGeneratorOps op)
		: sph_system_(sph_system), body_region_(body_name), 
		body_name_(body_name), base_particles_(base_particles), 
		refinement_level_(refinement_level), particle_generator_op_(op)
	{	
		sph_system_.AddBody(this);

		particle_spacing_ 	= RefinementLevelToParticleSpacing();
		kernel_ 			= sph_system_.GenerateAKernel(particle_spacing_);
		rst_step_ 			= sph_system_.restart_step_;
		mesh_cell_linked_list_
							= new MeshCellLinkedList(sph_system.lower_bound_,
									sph_system_.upper_bound_, kernel_->GetCutOffRadius());
		Vecd zero;
		by_cell_lists_particle_indexes_ 
							= new StdVec<IndexVector>[powern(3, zero.size())];
	}
	//===========================================================//
	Real SPHBody::RefinementLevelToParticleSpacing()
	{
		return sph_system_.particle_spacing_ref_	
			/ powern(2.0, refinement_level_);
	}
	//===========================================================//
	string SPHBody::GetBodyName()
	{
		return body_name_;
	}
	//===========================================================//
	void  SPHBody::AllocateMemoriesForConfiguration()
	{
		Allocate1dArray(indexes_contact_particles_, contact_map_.second.size());
		Allocate1dArray(reference_inner_configuration_, number_of_particles_);
		Allocate1dArray(current_inner_configuration_, number_of_particles_);
		Allocate2dArray(reference_contact_configuration_, {contact_map_.second.size(), number_of_particles_});
		Allocate2dArray(current_contact_configuration_, { contact_map_.second.size(), number_of_particles_ });
		//reserve memeory for concurrent vectors
		for (size_t i = 0; i != contact_map_.second.size(); ++i) {
			indexes_contact_particles_[i].reserve(number_of_particles_);
		}

	}
	//===========================================================//
	bool SPHBody::BodyContain(Vecd pnt)
	{
		return body_region_.contain(pnt);
	}
	//===========================================================//
	void SPHBody::ClosestPointOnBodySurface(Vecd input_pnt, Vecd& closest_pnt, Real& phi)
	{
		body_region_.closestpointonface(input_pnt, closest_pnt, phi);
	}
	//===========================================================//
	void SPHBody::BodyBounds(Vecd &lower_bound, Vecd &upper_bound)
	{
		body_region_.regionbound(lower_bound, upper_bound);
	}
	//===========================================================//
	void  SPHBody::SetContactMap(SPHBodyContactMap &contact_map)
	{
		contact_map_ = contact_map;
		//reserve memeory for concurrent vectors
		for (size_t i = 0; i != contact_map_.second.size(); ++i) {
			base_particles_contact_bodies_.push_back(&(contact_map_.second[i]->base_particles_));
		}
	}
	//===========================================================//
	void SPHBody::CreateParticelsInSpecificManner()
	{
		switch (particle_generator_op_)
		{
			case ParticlesGeneratorOps::lattice: {
				ParticleGeneratorLattice* lattice_particle_generator_ 
						= new ParticleGeneratorLattice(*this);
				lattice_particle_generator_->CreateParticles(*this);
				break;
			}

			case ParticlesGeneratorOps::direct: {
				ParticleGeneratorDirect* direct_particle_generator_ 
						= new ParticleGeneratorDirect(*this);
				direct_particle_generator_->CreateParticles(*this);
				break;
			}

			case ParticlesGeneratorOps::relax: {

				ReadRelaxedParticlsFromXmlFile* relax_particle_reader_
						= new ReadRelaxedParticlsFromXmlFile(*this);
				relax_particle_reader_->CreateParticles(*this);
				break;
			}

			case ParticlesGeneratorOps::restart: {
				ReadRestartParticlsFromXmlFile* restart_particle_reader_
						= new ReadRestartParticlsFromXmlFile(*this);
				restart_particle_reader_->CreateParticles(*this);
				break;
			}

			default: {
				std::cout << "\n FAILURE: the type of particle generator is undefined!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
				break;
			}
		}
	}
	//===========================================================//
	void SPHBody::GenerateAParticle(Vecd pnt, Real particle_volume)
	{
		base_particles_.InitializeAParticle(pnt, particle_volume);
	}
	//===========================================================//
	void SPHBody::WriteParticlesToVtuFile(ofstream &output_file)
	{
		base_particles_.WriteParticlesToVtuFile(output_file);
	}
	//===========================================================//
	void SPHBody::WriteParticlesToPltFile(ofstream &output_file)
	{
		base_particles_.WriteParticlesToPltFile(output_file);
	}
	//===========================================================//
	void SPHBody::WriteParticlesToXmlFile(std::string &filefullpath)
	{
		base_particles_.WriteParticlesToXmlFile(filefullpath);
	}
	//===============================================================//
	void SPHBody::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		base_particles_.WriteParticlesToXmlForRestart(filefullpath);
	}
	//===========================================================//
	SPHBody* SPHBody::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	void SPHBody::InitialConditionFromRestartFile()
	{
		std::string restart_folder_ = "./rstfile";
		if (!fs::exists(restart_folder_))
		{
			std::cout << "\n Error: the input file:"<< restart_folder_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		std::string filefullpath = restart_folder_ + "/SPHBody_" + GetBodyName() + "_rst_" + std::to_string(rst_step_) + ".xml";
		base_particles_.InitialParticleFromRestartXmlFile(filefullpath);
	}
	//===========================================================//
	void SPHBody::OffsetInitialParticlePosition(Vecd offset)
	{
		for (int i = 0; i < base_particles_.number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= base_particles_.base_particle_data_[i];

			base_particle_data_i.pos_n_ += offset;
		}
	}
	//===========================================================//
	RealBody::RealBody(SPHSystem &sph_system, string body_name,
		Particles &base_particles, int refinement_level, ParticlesGeneratorOps op)
		: SPHBody(sph_system, body_name, base_particles, refinement_level, op)
	{
		sph_system.AddRealBody(this);

		mesh_cell_linked_list_->AllocateMeshDataMatrix();
	}
	//===========================================================//
	void RealBody::AllocateMeoemryCellLinkedList()
	{
		mesh_cell_linked_list_->AllocateMeshDataMatrix();
	}
	//===========================================================//
	void RealBody::UpdateCellLinkedList()
	{
		mesh_cell_linked_list_->UpdateCellLists(*this);
	}
	//===========================================================//
	void RealBody::UpdateInnerConfiguration()
	{
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this);
	}
	//===========================================================//
	void RealBody::UpdateContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===========================================================//
	void RealBody
		::UpdateInteractionConfiguration(SPHBodyVector interacting_bodies)
	{
		mesh_cell_linked_list_
			->UpdateInteractionConfiguration(*this, interacting_bodies);
	}
	//===========================================================//
	void RealBody::SetAllParticleAtRest()
	{
		for (int i = 0; i < number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= base_particles_.base_particle_data_[i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
		}
	}
	//===========================================================//
	RealBody* RealBody::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	FictitiousBody::FictitiousBody(SPHSystem &system, string body_name,
		Particles &base_particles, 	int refinement_level, ParticlesGeneratorOps op)
		: SPHBody(system, body_name, base_particles, refinement_level, op)
	{
		system.AddFictitiousBody(this);
	}
	//===============================================================//
	void FictitiousBody::BuildInnerConfiguration()
	{
		UpdateInnerConfiguration();
	}
	//===========================================================//
	void FictitiousBody::UpdateCellLinkedList()
	{
		mesh_cell_linked_list_->UpdateParticleCellLocation(*this);
	}
	//===========================================================//
	void FictitiousBody::UpdateInnerConfiguration()
	{
		//do nothing
	}
	//===========================================================//
	void FictitiousBody::UpdateContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===========================================================//
	void FictitiousBody::UpdateInteractionConfiguration(SPHBodyVector interacting_bodies)
	{
		mesh_cell_linked_list_
			->UpdateInteractionConfiguration(*this, interacting_bodies);
	}
	//===========================================================//
	FictitiousBody* FictitiousBody::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	void FictitiousBody::GlobalBasicParameters(ofstream &out_file)
	{
		//noting done here
	}
}