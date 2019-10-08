/**
 * @file 	diffusion_dynamics.cpp
 * @brief 	In is file, we define functions decleared in diffusion_dynamics.h
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
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
		void OffsetInitialParticlePosition::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_particle_data_i = particles_->solid_body_data_[index_particle_i];

			base_particle_data_i.pos_n_ += offset_;
			solid_particle_data_i.pos_0_ += offset_;
		}
//=================================================================================================//
		void ElectroPhysiologyInitialCondition::Update(size_t index_particle_i, Real dt)
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
		}
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

            Real diff_trace = material_->getDiffusionTensorTrace(base_particle_data_i.pos_n_);
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

				Matd diff_tensor_i =  material_->getDiffussionTensor(base_particle_data_i.pos_n_);
				Matd diff_tensor_j =  material_->getDiffussionTensor(base_particle_data_j.pos_n_);
				Matd diff_ij = getAverageValue(diff_tensor_i, diff_tensor_j);
				Matd cd_diff_ij = inverseCholeskyDecomposition(diff_ij);
				Vecd grad_ij = cd_diff_ij * (-neighboring_particle.e_ij_);

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
            reaction_model_ = dynamic_cast<ElectrophysiologyReaction*>(body->base_reaction_->PointToThisObject());
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
            reaction_model_ = dynamic_cast<ElectrophysiologyReaction*>(body->base_reaction_->PointToThisObject());
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
		ApplyStimulusCurrents::ApplyStimulusCurrents(SolidBody *body)
				: ElectroPhysiologySimple(body) 
        {
            /* Nothing done here rightnow */
        }
//=================================================================================================//
        void ApplyStimulusCurrents::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData &base_particle_data_i      = particles_->base_particle_data_[index_particle_i];
            MuscleParticleData &muscle_particle_data_i  = particles_->muscle_body_data_[index_particle_i];

            Real physical_time = GlobalStaticVariables::physical_time_;
            Real stimulus_currents(0);
		    if(physical_time <= 0.5)
		    {
			    if(base_particle_data_i.pos_n_[0] <= 0.02)
			    {
				stimulus_currents = 0.92;
			    }
		    }

		    if(60.0 <= physical_time &&  physical_time <= 61.25)
		    {
			    if(base_particle_data_i.pos_n_[0] <= 0.5 )
			    {
				    if(base_particle_data_i.pos_n_[1] <= 0.5)
				    {
					stimulus_currents = 0.95;
				    }
			    }
		    }
			muscle_particle_data_i.voltage_n_ += stimulus_currents * dt;
		}
//=================================================================================================//
    }
//=================================================================================================//
}