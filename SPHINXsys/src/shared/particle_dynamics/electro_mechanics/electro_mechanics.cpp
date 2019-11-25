/**
 * @file 	electro_mechanics.cpp
 * @brief 	In is file, we define functions decleared in electro_mechanics.h
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
#include "electro_mechanics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "neighboring_particle.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "mesh_cell_linked_list.h"
#include "physiology_reaction.h"

using namespace SimTK;

namespace SPH
{
	namespace electro_mechanics
	{
        void OffsetInitialParticlePosition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_particle_data_i = particles_->solid_body_data_[index_particle_i];

			base_particle_data_i.pos_n_ += offset_;
			solid_particle_data_i.pos_0_ += offset_;
		}
//=================================================================================================//
		void ElectroMechanicsInitialCondition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_particle_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_particle_data_i = particles_->elastic_body_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			Vecd zero(0);
			base_particle_data_i.vel_n_ = zero;
			base_particle_data_i.dvel_dt_ = zero;

			solid_particle_data_i.vel_ave_ = zero;
			solid_particle_data_i.dvel_dt_ave_ = zero;

			elastic_particle_data_i.rho_0_ = material_->rho_0_;
			elastic_particle_data_i.rho_n_ = material_->rho_0_;
			elastic_particle_data_i.mass_ = material_->rho_0_*base_particle_data_i.Vol_;

            muscle_particle_data_i.voltage_n_ = 0.0;
		    muscle_particle_data_i.dvoltage_dt_ = 0.0;
		    muscle_particle_data_i.grad_voltage_ = 0.0;
		    muscle_particle_data_i.gate_var_ = 0.0;
            muscle_particle_data_i.T_a_ = 0.0;
			muscle_particle_data_i.active_stress_ = Matd(0.0);
		}
//=================================================================================================//
        computeActiveContractionStress::computeActiveContractionStress(SolidBody *body)
            : ElectroMechanicsSimple(body) 
        {
            reaction_model_ = dynamic_cast<ElectroPhysiology*>(body->base_reaction_->PointToThisObject());
        }
//=================================================================================================//
        void computeActiveContractionStress::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			/**
			 * We assume that f(y) = q - p y, then using the quasi-steady_state for ODE.
			 */
			Real q = reaction_model_->getProductionRateOfActiveContractionStress(muscle_particle_data_i.voltage_n_);
			Real p = reaction_model_->getLoseRateOfActiveContractionStress(muscle_particle_data_i.voltage_n_);
			Real T_a = muscle_particle_data_i.T_a_;
			muscle_particle_data_i.T_a_ = T_a * exp(-p * dt) + q * (1.0 - exp(-p * dt)) / p;
			muscle_particle_data_i.active_stress_ = 
				material_->ConstitutiveRelationOfActiveStress(elastic_data_i.F_, muscle_particle_data_i.T_a_, index_particle_i);
		}
//=================================================================================================//
        computeLinearActiveContractionStress::computeLinearActiveContractionStress(SolidBody *body)
            : ElectroMechanicsSimple(body) 
        {
            // We assume the active contraction stree is linear ot transmembrane potential
        }
//=================================================================================================//
        void computeLinearActiveContractionStress::Update(size_t index_particle_i, Real dt)
		{
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			muscle_particle_data_i.T_a_ = dt * muscle_particle_data_i.voltage_n_;
			muscle_particle_data_i.active_stress_ = 
				material_->ConstitutiveRelationOfActiveStress(elastic_data_i.F_, muscle_particle_data_i.T_a_, index_particle_i);
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationFirstStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];
			MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_*dt * 0.5;
			elastic_data_i.rho_n_ = elastic_data_i.rho_0_ / det(elastic_data_i.F_);
			elastic_data_i.stress_ = material_->ConstitutiveRelation(elastic_data_i.F_, index_particle_i)
				+ material_->NumericalDampingStress(elastic_data_i.F_, elastic_data_i.dF_dt_, numerical_viscosity_) 
				+ muscle_particle_data_i.active_stress_;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_*dt * 0.5;
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationFirstStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];
			MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			//including gravity and force from fluid
			// Vecd acceleration = base_particle_data_i.dvel_dt_others_ 
			// 	+ solid_data_i.force_from_fluid_/ elastic_data_i.mass_;
			Vecd acceleration(0);

			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				acceleration += (elastic_data_i.stress_ * solid_data_i.B_
					+ elastic_data_j.stress_* solid_data_j.B_)
					* neighboring_particle.dW_ij_ * neighboring_particle.e_ij_
					* base_particle_data_j.Vol_ / elastic_data_i.rho_0_;
			}
			base_particle_data_i.dvel_dt_ = acceleration;
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationFirstStep::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_ += base_particle_data_i.dvel_dt_* dt;
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationSecondStep::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i 	= particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationSecondStep::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];
			MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			Matd deformation_gradient_change_rate(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				SolidParticleData &solid_data_j = particles_->solid_body_data_[index_particle_j];
				ElasticSolidParticleData &elastic_data_j = particles_->elastic_body_data_[index_particle_j];

				//deformtion
				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				deformation_gradient_change_rate
					-= base_particle_data_j.Vol_
					*SimTK::outer((base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_), gradw_ij);
			}
			elastic_data_i.dF_dt_ = deformation_gradient_change_rate* solid_data_i.B_;
		}
//=================================================================================================//
		void ActivePassiveStressRelaxationSecondStep::Update(size_t index_particle_i, Real dt)
		{
			ElasticSolidParticleData &elastic_data_i = particles_->elastic_body_data_[index_particle_i];

			elastic_data_i.F_ += elastic_data_i.dF_dt_ * dt * 0.5;
		}
//=================================================================================================//
		Vecd SpringConstrainMuscleRegion::GetAcceleration(Vecd &disp, Real mass)
		{
			Vecd spring_force(0);
			for(int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i] / mass;
			}
			return spring_force;
		}
//=================================================================================================//
		void SpringConstrainMuscleRegion::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i 
				= particles_->elastic_body_data_[index_particle_i];

			Vecd disp_from_0 = base_particle_data_i.pos_n_ - solid_data_i.pos_0_;
			base_particle_data_i.vel_n_ 	+=  dt * GetAcceleration(disp_from_0, elastic_data_i.mass_);
			base_particle_data_i.pos_n_ 	+=  dt * dt * GetAcceleration(disp_from_0, elastic_data_i.mass_);
			//base_particle_data_i.dvel_dt_ 	+=  GetAcceleration(disp_from_0, elastic_data_i.mass_);
		}
//=================================================================================================//
    }
//=================================================================================================//
}
//=================================================================================================//