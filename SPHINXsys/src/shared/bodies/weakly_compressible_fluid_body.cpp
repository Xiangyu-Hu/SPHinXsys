#include "weakly_compressible_fluid_body.h"
#include "mesh_cell_linked_list.h"
#include "weakly_compressible_fluid.h"
#include "weakly_compressible_fluid_particles.h"
#include "sph_system.h"
#include "base_kernel.h"

namespace SPH {
	//===============================================================//
	WeaklyCompressibleFluidBody
		::WeaklyCompressibleFluidBody(SPHSystem &system, string body_name,
		WeaklyCompressibleFluid* material, 
		WeaklyCompressibleFluidParticles &weakly_compressible_fluid_particles,
			int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, weakly_compressible_fluid_particles, 
			refinement_level, op), material_(material),
		weakly_compressible_fluid_particles_(weakly_compressible_fluid_particles)
	{

	}
	//===============================================================//
	void WeaklyCompressibleFluidBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this);
	}
	//===============================================================//
	void WeaklyCompressibleFluidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
	void WeaklyCompressibleFluidBody::SetAllParticleAtRest()
	{
		for (int i = 0; i < number_of_particles_; ++i) {
			BaseParticleData &base_particle_data_i
				= weakly_compressible_fluid_particles_.base_particle_data_[i];
			WeaklyCompressibleFluidParticleData &fluid_data_i
				= weakly_compressible_fluid_particles_.fluid_data_[i];

			fluid_data_i.p_ = 0.0; 
			base_particle_data_i.vel_n_(0);
			fluid_data_i.vel_trans_(0);
			base_particle_data_i.dvel_dt_(0);
			fluid_data_i.rho_0_
				= material_->ReinitializeRho(fluid_data_i.p_);
			fluid_data_i.rho_n_ = fluid_data_i.rho_0_;
			fluid_data_i.mass_
				= fluid_data_i.rho_0_*base_particle_data_i.Vol_;
		}
	}
	//===============================================================//
	void WeaklyCompressibleFluidBody::GlobalBasicParameters(ofstream &out_file)
	{
		out_file << "Sphbody_Body_Name :" << "   " << body_name_ << "\n";
		out_file << "Sph_System_Particle_Spacing_Ref :" << "   " << sph_system_.particle_spacing_ref_ << "\n";
		out_file << "Sph_System_Lower_Bound :" << "   " << sph_system_.lower_bound_ << "\n";
		out_file << "Sph_System_Upper_Bound :" << "   " << sph_system_.upper_bound_ << "\n";
		out_file << "Refinement_Level :" << "   " << refinement_level_ << "\n";
		out_file << "Speed_Max :" << "   " << speed_max_ << "\n";
		out_file << "Base_Particles_.Particle_Spacing :" << "   " << particle_spacing_ << "\n";
		out_file << "Reference Density :" << "   " << material_->rho_0_ << "\n";
		out_file << "Viscosity :" << "   " << material_->mu_ << "\n";
		out_file << "Sound Speed : " << "   " << material_->c_0_ << "\n";
		out_file << "Thermal Condution rate : " << "   " << material_->k_ << "\n";
		out_file << "Signal_Speed_max : " << "   " << signal_speed_max_ << "\n";
		out_file << "Number_Of_Particles : " << "   "<< number_of_particles_ << "\n";
	}
	//===============================================================//
	FluidBodyPart
		::FluidBodyPart(WeaklyCompressibleFluidBody *fluid_body, string fluid_body_part_name)
		: EulerianBodyPart(fluid_body, fluid_body_part_name),
		fluid_body_(fluid_body), fluid_body_part_region_(fluid_body_part_name)
	{

	}
	//===============================================================//
	Oldroyd_B_FluidBody::
		Oldroyd_B_FluidBody(SPHSystem &system, string body_name,
			Oldroyd_B_Fluid* oldroyd_b_material,
			Oldroyd_B_FluidParticles &oldroyd_b_fluid_particles,
			int refinement_level, ParticlesGeneratorOps op)
		: WeaklyCompressibleFluidBody(system, body_name, oldroyd_b_material, 
			oldroyd_b_fluid_particles, refinement_level, op),
		oldroyd_b_material_(oldroyd_b_material),
		oldroyd_b_fluid_particles_(oldroyd_b_fluid_particles)
	{

	}
	//===============================================================//
	void Oldroyd_B_FluidBody::GlobalBasicParameters(ofstream &out_file)
	{
		WeaklyCompressibleFluidBody::GlobalBasicParameters(out_file);
		out_file << "Relaxation Time: " << "   " << oldroyd_b_material_->lambda_ << "\n";
		out_file << "Olymeric Viscosity : " << "   " << oldroyd_b_material_->mu_p_ << "\n";
	}
	//===============================================================//
}