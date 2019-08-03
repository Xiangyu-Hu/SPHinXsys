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
#include "in_output.h"
#include "base_particles.h"
#include "base_material.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"
#include "particle_generator_lattice.h"

namespace SPH
{
	//===========================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, string body_name, Material &base_material,
		int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op)
		: sph_system_(sph_system), body_region_(body_name), body_name_(body_name), base_material_(base_material),
		refinement_level_(refinement_level), smoothinglength_ratio_(smoothinglength_ratio), particle_generator_op_(op)
	{	
		sph_system_.AddBody(this);

		particle_spacing_ 	= RefinementLevelToParticleSpacing();
		kernel_ 			= sph_system_.GenerateAKernel(particle_spacing_*smoothinglength_ratio_);
		mesh_cell_linked_list_
							= new MeshCellLinkedList(sph_system.lower_bound_,
									sph_system_.upper_bound_, kernel_->GetCutOffRadius());
		number_of_by_cell_lists_ = powern(3, Vecd(0).size());
		/** I will use concurrent vector here later after tests. */
		by_cell_lists_particle_indexes_ = new StdVec<IndexVector>[number_of_by_cell_lists_];
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
		Allocate2dArray(reference_contact_configuration_, {contact_map_.second.size(), number_of_particles_ });
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
	}
	//===========================================================//
	void SPHBody::WriteParticlesToVtuFile(ofstream &output_file)
	{
		base_particles_->WriteParticlesToVtuFile(output_file);
	}
	//===========================================================//
	void SPHBody::WriteParticlesToPltFile(ofstream &output_file)
	{
		base_particles_->WriteParticlesToPltFile(output_file);
	}
	//===============================================================//
	void SPHBody::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		base_particles_->WriteParticlesToXmlForRestart(filefullpath);
	}
	//===========================================================//
	void SPHBody::ReadParticlesFromXmlForRestart(std::string &filefullpath) 
	{
		base_particles_->ReadParticleFromXmlForRestart(filefullpath);
	}
	//===============================================================//
	void SPHBody::WriteToXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->WriteToXmlForReloadParticle(filefullpath);
	}
	//===========================================================//
	void SPHBody::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->ReadFromXmlForReloadParticle(filefullpath);
	}
	//===========================================================//
	SPHBody* SPHBody::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	RealBody::RealBody(SPHSystem &sph_system, string body_name, Material &base_material,
		int refinement_level, Real smoothinglength_ratio, ParticlesGeneratorOps op)
		: SPHBody(sph_system, body_name, base_material, refinement_level, smoothinglength_ratio, op)
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
	RealBody* RealBody::PointToThisObject()
	{
		return this;
	}
	//===========================================================//
	FictitiousBody::FictitiousBody(SPHSystem &system, string body_name, int refinement_level, 
		Real smoothinglength_ratio, ParticlesGeneratorOps op)
		: SPHBody(system, body_name, *(new Material("FictitiousMaterial")), refinement_level, smoothinglength_ratio, op)
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
}