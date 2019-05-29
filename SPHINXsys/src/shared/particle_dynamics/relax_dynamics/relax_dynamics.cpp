#include "relax_dynamics.h"
#include "relax_body.h"
#include "relax_body_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"

#include "math.h"

using namespace SimTK;

namespace SPH
{
	namespace relax_dynamics
	{
		//===========================================================//
		GetTimeStepSize::GetTimeStepSize(RelaxBody* body)
			: ParticleDynamicsMinimum<RelaxBody, RelaxBodyParticles>(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			Real eta = 0.1 * 1.3 * body->particle_spacing_;
			initial_reference_ = 0.125 * smoothing_length_* smoothing_length_ * 1.0 / (eta + 1.0e-15);
		}
		//===========================================================//
		Real GetTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			return 0.25 * SMIN(sqrt(smoothing_length_ / (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				smoothing_length_ / (40 * base_particle_data_i.vel_n_.norm() + 1.0e-15));
		}
		//===========================================================//
		PhysicsRelaxationInner::PhysicsRelaxationInner(RelaxBody *body)
			: ParticleDynamicsInner1Level<RelaxBody, RelaxBodyParticles>(body)
		{
			smoothing_length_ = body_->kernel_->GetSmoothingLength();
			sound_speed_ = 1.0;

			//time setep size due to linear viscosity
			eta_ = 0.1 * 1.3 * body_->particle_spacing_;
			p_star_ = 2.0;
		}
		//===========================================================//
		void PhysicsRelaxationInner::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void PhysicsRelaxationInner::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData &relax_data_i
				= particles_->relax_body_data_[index_particle_i];

			//including gravity and force from fluid
			Vecd acceleration(0);

			StdVec<NeighboringParticle>  &neighors
				= (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;

				BaseParticleData &base_particle_data_j
					= particles_->base_particle_data_[index_particle_j];
				RelaxBodyParticleData &relax_data_j
					= particles_->relax_body_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * neighboring_particle.gradW_ij_
					*base_particle_data_j.Vol_ * base_particle_data_i.Vol_;
				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle.r_ij_.norm() + 0.01 * smoothing_length_ );
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_* neighboring_particle.dW_ij_;
				if(isnan(acceleration.norm()))
				{
					std::cout << "\n Error: the viscous acceleration is not valid" << std::endl;
					std::cout << "\n Particle ID :" << index_particle_i << "\n Particle ID :" << index_particle_j << std::endl;
					std::cout << "\n Particle distance :" << neighboring_particle.r_ij_.norm()<< std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
			}

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			//Vecd norm  = SGN(phi) *  elastic_body_.mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			Vecd dist_2_face(0);
			Vecd norm(0);
			dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			norm = -dist_2_face / (dist_2_face.norm() + 1.0e-15);
			
			if ((phi - 0.5 * body_->particle_spacing_) <= smoothing_length_)
			{
				Real sigma_ = (phi - 0.5 * body_->particle_spacing_) / smoothing_length_;
				Real gamma_ = 1.0 + (0.0625 - 0.0531 * sigma_) * (sigma_ - 2.0) * (sigma_ - 2.0) * (sigma_ - 2.0);
				Real grad_gamma_ = (0.2937 - 0.2124*sigma_) * (sigma_ - 2.0) * (sigma_ - 2.0);
				acceleration += p_star_ * base_particle_data_i.Vol_ / smoothing_length_ / gamma_
					* grad_gamma_ * norm;
				if ((phi - 0.5 * body_->particle_spacing_) <= back_mesh_spacing_)
				{
					acceleration -= dot(acceleration, norm) * norm;
				}
			}
			
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void PhysicsRelaxationInner::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData &relax_data_i
				= particles_->relax_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ = base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		PhysicsRelaxationComplex::PhysicsRelaxationComplex(RelaxBody *body, StdVec<RelaxBody*> interacting_bodies)
			: ParticleDynamicsComplex1Level<RelaxBody, RelaxBodyParticles, RelaxBody, RelaxBodyParticles>(body, interacting_bodies)
		{
			Real smoothing_length_ = body_->kernel_->GetSmoothingLength();
			eta_ = 0.25 * smoothing_length_;
			p_star_ = 1.0;
		}
		//===========================================================//
		void PhysicsRelaxationComplex::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//===========================================================//
		void PhysicsRelaxationComplex::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData &relax_data_i = particles_->relax_body_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];

			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				RelaxBodyParticleData &relax_data_j = particles_->relax_body_data_[index_particle_j];

				acceleration -= 2.0*p_star_*base_particle_data_j.Vol_ * base_particle_data_i.Vol_
					*neighboring_particle.gradW_ij_;
				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ neighboring_particle.r_ij_.norm();
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * neighboring_particle.dW_ij_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
		//===========================================================//
		void PhysicsRelaxationComplex
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData &relax_data_i = particles_->relax_body_data_[index_particle_i];

			Vecd acceleration(0);
			StdVec<NeighboringParticle>  &neighors = (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				RelaxBodyParticleData &relax_data_j
					= (*interacting_particles_[interacting_body_index]).relax_body_data_[index_particle_j];

				//pressure force
				acceleration -= 2.0 * p_star_ * base_particle_data_j.Vol_ * base_particle_data_i.Vol_
					* neighboring_particle.gradW_ij_;

				//viscous force
				Vecd vel_detivative = 2.0*base_particle_data_i.vel_n_ / neighboring_particle.r_ij_.norm();
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * neighboring_particle.dW_ij_;
			}
			base_particle_data_i.dvel_dt_ += acceleration;
		}
		//===========================================================//
		void PhysicsRelaxationComplex::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_ = base_particle_data_i.dvel_dt_* dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;

		}
		//===========================================================//
		void ConstriantSurfaceParticles ::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			Vecd dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			//base_particle_data_i.pos_n_ += dist_2_face;	
			Vecd norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
			if (phi >= 0.0) {
				base_particle_data_i.pos_n_ += (phi - 0.5 * body_->particle_spacing_) * norm;
			}
			else {
				base_particle_data_i.pos_n_ += (ABS(phi) + 0.5 * body_->particle_spacing_) * norm;
			}
		}
		//===========================================================//
	}
}
