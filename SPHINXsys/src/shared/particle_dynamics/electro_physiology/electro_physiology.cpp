/**
 * @file 	electro_physiology.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */

#include "electro_physiology.h"
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
	namespace electro_physiology
	{
//=================================================================================================//
		getDiffusionTimeStepSize::getDiffusionTimeStepSize(SolidBody* body)
			: ElectroPhysiologyMinimum(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//time setep size due to unit diffusion parameter
			initial_reference_ = 0.25 *smoothing_length_;
		}
//=================================================================================================//
		Real getDiffusionTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

            Real diff_trace = material_->getDiffusionTensorTrace(index_particle_i);
			return SMIN(initial_reference_, 0.125 *smoothing_length_ * smoothing_length_ / diff_trace );
		}
//=================================================================================================//
		void CorrectConfiguration::Initialization(size_t index_particle_i, Real dt)
		{
            /* Nothing done here rightnow. */
		}
//=================================================================================================//
		void CorrectConfiguration::InnerInteraction(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];

			Matd local_configuration(0.0);
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				Vecd gradw_ij = neighboring_particle.dW_ij_ * neighboring_particle.e_ij_;
				Vecd r_ij = -neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ij, gradw_ij);
			}
			solid_data_i.temp_matrix_ = local_configuration;
		}
//=================================================================================================//
		void CorrectConfiguration::Update(size_t index_particle_i, Real dt)
		{
			SolidParticleData &solid_data_i = particles_->solid_body_data_[index_particle_i];
			/** note that the generalized inverse only works here*/
			solid_data_i.B_ = GeneralizedInverse(solid_data_i.temp_matrix_);
		}
//=================================================================================================//
		void DiffusionRelaxation::Initialization(size_t index_particle_i, Real dt)
		{
            /* Nothing done here rightnow. */
		}
//=================================================================================================//
		void DiffusionRelaxation::InnerInteraction(size_t index_particle_i, Real dt)
		{
            BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_particle_data_i = particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_particle_data_i = particles_->elastic_body_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			Real d_d_ = 0.0;
			StdVec<ReferenceNeighboringParticle>  &neighors = (*reference_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;

                BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
			    SolidParticleData &solid_particle_data_j = particles_->solid_body_data_[index_particle_j];
			    ElasticSolidParticleData &elastic_particle_data_j = particles_->elastic_body_data_[index_particle_j];
                MuscleParticleData &muscle_particle_data_j = particles_->muscle_body_data_[index_particle_j];

				// Matd cd_diff_ij = material_->getReferenceAverageDiffusionTensor(index_particle_i, n);
				// Vecd grad_ij = cd_diff_ij * (-neighboring_particle.e_ij_);
				Matd diff_cd_i = material_->getDiffussionTensor(index_particle_i);
				Matd diff_cd_j = material_->getDiffussionTensor(index_particle_j);
				Matd diff_cd_ij = getAverageValue(diff_cd_i, diff_cd_j);
				Vecd grad_ij = diff_cd_ij * (-neighboring_particle.e_ij_);

				Matd k_ij  = solid_particle_data_i.B_ * muscle_particle_data_i.voltage_n_ 
                           - solid_particle_data_j.B_ * muscle_particle_data_j.voltage_n_;
				Vecd v_ij = k_ij * (-neighboring_particle.e_ij_);
				Real grad_v_ij = dot(v_ij, -neighboring_particle.e_ij_);

				d_d_ += 2.0 * grad_v_ij * base_particle_data_j.Vol_ * neighboring_particle.dW_ij_ / 
					neighboring_particle.r_ij_ / dot(grad_ij,grad_ij);
			}
			muscle_particle_data_i.dvoltage_dt_ = d_d_;
		}
//=================================================================================================//
		void DiffusionRelaxation::Update(size_t index_particle_i, Real dt)
		{
			MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];
			/** note that the generalized inverse only works here*/
			muscle_particle_data_i.voltage_n_ += dt * muscle_particle_data_i.dvoltage_dt_;
		}
//=================================================================================================//
        TransmembranePotentialReaction::TransmembranePotentialReaction(SolidBody *body)
            : ElectroPhysiologySimple(body) 
        {
            reaction_model_ = dynamic_cast<ElectroPhysiology*>(body->base_reaction_->PointToThisObject());
        }
//=================================================================================================//
        void TransmembranePotentialReaction::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			/**
			 * We assume that f(y) = q - p y, then using the quasi-steady_state for ODE.
			 */
			Real q = reaction_model_->getProductionRateOfIonicCurrent(muscle_particle_data_i.voltage_n_, 
                                        muscle_particle_data_i.gate_var_);
			Real p = reaction_model_->getLoseRateOfIonicCurrent(muscle_particle_data_i.voltage_n_, 
                                        muscle_particle_data_i.gate_var_);
			Real voltage = muscle_particle_data_i.voltage_n_;
			muscle_particle_data_i.voltage_n_ = voltage * exp(-p * dt) + q * (1.0 - exp(-p * dt)) / p;
		}
//=================================================================================================//
		GateVariableReaction::GateVariableReaction(SolidBody *body)
				: ElectroPhysiologySimple(body) 
        {
            reaction_model_ = dynamic_cast<ElectroPhysiology*>(body->base_reaction_->PointToThisObject());
        }
//=================================================================================================//
        void GateVariableReaction::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

			Real q = reaction_model_->getProductionRateOfGateVarible(muscle_particle_data_i.voltage_n_, 
                                    muscle_particle_data_i.gate_var_);
			Real p = reaction_model_->getLoseRateOfGateVarible(muscle_particle_data_i.voltage_n_, 
                                    muscle_particle_data_i.gate_var_);
			Real gate = muscle_particle_data_i.gate_var_;
			muscle_particle_data_i.gate_var_ = gate * exp(-p * dt) + q * (1.0 - exp(-p * dt)) / p;
		}
//=================================================================================================//
    }
//=================================================================================================//
}