/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "solid_body.h"
#include "mesh_cell_linked_list.h"
#include "elastic_solid.h"
#include "solid_body_particles.h"
#include "sph_system.h"
#include "base_kernel.h"

namespace SPH {
	//===============================================================//
	SolidBody::SolidBody(SPHSystem &system, string body_name, 
		SolidBodyParticles &solid_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, solid_particles, refinement_level, op),
		solid_particles_(solid_particles)
	{
	
	}
	//===============================================================//
	void SolidBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceInnerConfiguration(*this);
	}
	//===============================================================//
	void SolidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceContactConfiguration(*this);
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
	void SolidBody::SetAllParticleAtRest()
	{
		for (int i = 0; i < solid_particles_.number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= solid_particles_.base_particle_data_[i];
			SolidBodyParticleData &solid_body_data_i
				= solid_particles_.solid_body_data_[i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
			base_particle_data_i.dvel_dt_ = zero;
			solid_body_data_i.vel_ave_ = zero;
			solid_body_data_i.dvel_dt_ave_ = zero;
		}
	}
	//===========================================================//
	void SolidBody::OffsetInitialParticlePosition(Vecd offset)
	{
		for (int i = 0; i < solid_particles_.number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= solid_particles_.base_particle_data_[i];
			SolidBodyParticleData &solid_particle_data_i
				= solid_particles_.solid_body_data_[i];

			base_particle_data_i.pos_n_ += offset;
			solid_particle_data_i.pos_0_ += offset;
		}
	}
	//===============================================================//
	void SolidBody::GlobalBasicParameters(ofstream &out_file)
	{
		out_file << "Body_Name :" << "   " << body_name_ << "\n";
		out_file << "Refinement_Level :" << "   " << refinement_level_ << "\n";
		out_file << "Particle_Spacing :" << "   " << particle_spacing_ << "\n";
		out_file << "Number_Of_Particles : " << "   "<< number_of_particles_ << "\n";
	}
	//===============================================================//
	SolidBodyPart::SolidBodyPart(SolidBody *solid_body, string soild_body_part_name)
		: LagrangianBodyPart(solid_body, soild_body_part_name),
		solid_body_(solid_body), soild_body_part_region_(soild_body_part_name)
	{

	}
	//===============================================================//
	void SolidBodyPart::TagBodyPartParticles()
	{
		for (size_t i = 0; i < solid_body_->number_of_particles_; ++i)
		{
			BaseParticleData &base_particle_data_i
				= solid_body_->base_particles_.base_particle_data_[i];

			if (soild_body_part_region_.contain(base_particle_data_i.pos_n_))
			{
				body_part_particles_.push_back(i);
			}
		}
	}
	//===============================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SolidBody *solid_body,
			string soild_body_part_name, Real solid_body_density)
		: SolidBodyPart(solid_body, soild_body_part_name), 
		solid_body_density_(solid_body_density)
	{

	}
	//===============================================================//
	ElasticBody::ElasticBody(SPHSystem &system, string body_name,
		ElasticSolid* material, ElasticBodyParticles &elastic_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, elastic_particles, refinement_level, op),
		material_(material), elastic_particles_(elastic_particles)
	{

	}
	//===============================================================//
	void ElasticBody::SetAllParticleAtRest()
	{
		for (int i = 0; i < elastic_particles_.number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= elastic_particles_.base_particle_data_[i];
			SolidBodyParticleData &solid_particle_data_i
				= elastic_particles_.solid_body_data_[i];
			ElasticBodyParticleData &elastic_particle_data_i
				= elastic_particles_.elastic_body_data_[i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
			elastic_particle_data_i.rho_0_ = material_->rho_0_;
			elastic_particle_data_i.rho_n_ = material_->rho_0_;
			elastic_particle_data_i.mass_
				= material_->rho_0_*base_particle_data_i.Vol_;
		}
	}
	//===============================================================//
	void ElasticBody::InitializeLocalMaterialProperties()
	{
		for (int i = 0; i < number_of_particles_; ++i) {
			ElasticBodyParticleData &elastic_data_i
				= elastic_particles_.elastic_body_data_[i];

			elastic_data_i.local_G_
				= material_->GetShearModulus(material_->E_0_, material_->nu_);
			elastic_data_i.local_c_
				= material_->GetSoundSpeed(material_->rho_0_, material_->E_0_, material_->nu_);
			elastic_data_i.local_lambda_
				= material_->GetLambda(material_->E_0_, material_->nu_);
			elastic_data_i.local_eta_
				= material_->eta_0_ + material_->GetArtificalViscosity(material_->rho_0_,
					elastic_data_i.local_c_, kernel_->GetSmoothingLength());
		}
	}
	//===============================================================//
	void ElasticBody::GlobalBasicParameters(ofstream &out_file)
	{
		SolidBody::GlobalBasicParameters(out_file);
		out_file << "Youngs Modulus : " << "   " << material_->E_0_ << "\n";
		out_file << "Poisson Ration : " << "   " << material_->nu_ << "\n";
		out_file << "Physical Viscosity : " << "   " << material_->eta_0_ << "\n";
		out_file << "Reference Density : " << "   " << material_->rho_0_ << "\n";

	}
	//===============================================================//
	Matd ElasticBody::GetElasticStress(Matd &deform_grad, size_t index_particle_i)
	{
		ElasticBodyParticleData &elastic_data_i
			= elastic_particles_.elastic_body_data_[index_particle_i];
		return material_->ConstitutiveRelation(deform_grad,
			elastic_data_i.local_G_, elastic_data_i.local_lambda_);
	}
	//===============================================================//
	MuscleBody::MuscleBody(SPHSystem &system, string body_name,
		Muscle* material, MuscleBodyParticles &muscle_particles, 
		int refinement_level, ParticlesGeneratorOps op)
		: ElasticBody(system, body_name, material, muscle_particles, refinement_level, op),
		muscle_particles_(muscle_particles_)
	{

	}
	//===============================================================//
	void MuscleBody::InitializeLocalMaterialProperties()
	{
		ElasticBody::InitializeLocalMaterialProperties();

		for (int i = 0; i < number_of_particles_; ++i) {
			MuscleBodyData &muscle_data_i
				= muscle_particles_.muscle_body_data_[i];
			
			for (size_t i = 0; i != 4; ++i)
			{
				muscle_data_i.local_a_[i] = material_->a_0_[i];
				muscle_data_i.local_b_[i] = material_->b_0_[i];
			}

		}
	}
	//===============================================================//
	void MuscleBody::GlobalBasicParameters(ofstream &out_file)
	{
		ElasticBody::GlobalBasicParameters(out_file);
		out_file << "Consitutive Parameters : " << "   " << material_->a_0_ << "\n";
		out_file << "Consitutive Parameters : " << "   " << material_->b_0_ << "\n";
	}
	//===============================================================//
	Matd MuscleBody::GetElasticStress(Matd &deform_grad, size_t index_particle_i)
	{
		ElasticBodyParticleData &elastic_data_i
			= elastic_particles_.elastic_body_data_[index_particle_i];
		MuscleBodyData &muscle_data_i
			= muscle_particles_.muscle_body_data_[index_particle_i];

		return material_->ConstitutiveRelation(deform_grad,
			elastic_data_i.local_lambda_, muscle_data_i.local_a_, muscle_data_i.local_b_,
			muscle_data_i.local_f0_, muscle_data_i.local_s0_, muscle_data_i.local_f0f0_,
			muscle_data_i.local_s0s0_, muscle_data_i.local_f0s0_);
	}
	//===============================================================//
}