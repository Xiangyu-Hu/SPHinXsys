/**
 * @file 	eulerian_compressible_fluid_dynamics_complex.hpp
 * @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
 */
#include "eulerian_compressible_fluid_dynamics_complex.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace eulerian_compressible_fluid_dynamics
	{
		//=================================================================================================//
		template<class BaseRelaxationType>
		template<class BaseBodyRelationType>
		RelaxationWithWall<BaseRelaxationType>::
            RelaxationWithWall(BaseBodyRelationType &base_body_relation, 
            BaseBodyRelationContact &wall_contact_relation) :
				BaseRelaxationType(base_body_relation), CompressibleFluidWallData(wall_contact_relation)
		{
			if (base_body_relation.sph_body_ != wall_contact_relation.sph_body_)
			{
				std::cout << "\n Error: the two body_realtions do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
	
			for (size_t k = 0; k != CompressibleFluidWallData::contact_particles_.size(); ++k)
			{
				Real rho_0_k = CompressibleFluidWallData::contact_particles_[k]->rho0_;
				wall_inv_rho_0_.push_back(1.0 / rho_0_k);
				wall_mass_.push_back(&(CompressibleFluidWallData::contact_particles_[k]->mass_));
				wall_Vol_.push_back(&(CompressibleFluidWallData::contact_particles_[k]->Vol_));
				wall_vel_ave_.push_back(&(CompressibleFluidWallData::contact_particles_[k]->vel_ave_));
				wall_dvel_dt_ave_.push_back(&(CompressibleFluidWallData::contact_particles_[k]->dvel_dt_ave_));
				wall_n_.push_back(&(CompressibleFluidWallData::contact_particles_[k]->n_));
			}
		}
        //=================================================================================================//
        template<class BaseViscousAccelerationType>   	
		template<class BaseBodyRelationType>
		ViscousWithWall<BaseViscousAccelerationType>::
            ViscousWithWall(BaseBodyRelationType &base_body_relation, 
				BaseBodyRelationContact &wall_contact_relation) 
		: RelaxationWithWall<BaseViscousAccelerationType>(base_body_relation, wall_contact_relation) {}
		//=================================================================================================//
        template<class BaseViscousAccelerationType>
		void ViscousWithWall<BaseViscousAccelerationType>::Interaction(size_t index_i, Real dt)
		{
			BaseViscousAccelerationType::Interaction(index_i, dt);
			
			Real rho_i = this->rho_n_[index_i];
			Vecd& vel_i = this->vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			for (size_t k = 0; k < CompressibleFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->wall_Vol_[k]);
				StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*CompressibleFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0*(vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					acceleration += 2.0 * this->mu_ * vel_derivative 
								  * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
				}
			}
			
			this->dmom_dt_prior_[index_i] += rho_i * acceleration;
			this->dE_dt_prior_[index_i] += rho_i * dot(acceleration, vel_i);
		}
		//=================================================================================================//
        template<class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(ComplexBodyRelation &fluid_wall_relation) :
				BaseViscousAccelerationType(fluid_wall_relation.inner_relation_,
					fluid_wall_relation.contact_relation_) {}
        //=================================================================================================//
        template<class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(BaseBodyRelationInner &fluid_inner_relation, 
				BaseBodyRelationContact &wall_contact_relation) :
				BaseViscousAccelerationType(fluid_inner_relation,
					wall_contact_relation) {}
        //=================================================================================================//
        template<class BaseViscousAccelerationType>
		BaseViscousAccelerationWithWall<BaseViscousAccelerationType>::
			BaseViscousAccelerationWithWall(ComplexBodyRelation &fluid_complex_relation, 
				BaseBodyRelationContact &wall_contact_relation) :
				BaseViscousAccelerationType(fluid_complex_relation,
					wall_contact_relation) {}
       //=================================================================================================//
        template<class BasePressureRelaxationType>   	
		template<class BaseBodyRelationType>
		PressureRelaxation<BasePressureRelaxationType>::
            PressureRelaxation(BaseBodyRelationType &base_body_relation, 
				BaseBodyRelationContact &wall_contact_relation) :
				RelaxationWithWall<BasePressureRelaxationType>(base_body_relation, 
					wall_contact_relation) {}
       //=================================================================================================//
        template<class BasePressureRelaxationType>   	
        void PressureRelaxation<BasePressureRelaxationType>::Interaction(size_t index_i, Real dt)
		{
			BasePressureRelaxationType::Interaction(index_i, dt);

			CompressibleFluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i], this->E_[index_i]);
			Vecd momentum_change_rate(0.0);
			for (size_t k = 0; k < CompressibleFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->wall_Vol_[k]);
				StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
				Neighborhood& wall_neighborhood = (*CompressibleFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd& e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ij = wall_neighborhood.dW_ij_[n];
					Real r_ij = wall_neighborhood.r_ij_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;
					Real E_in_wall = state_i.E_;
					CompressibleFluidState state_j(rho_in_wall, vel_in_wall, p_in_wall, E_in_wall);
					CompressibleFluidState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Vecd vel_star = interface_state.vel_;
					Real p_star = interface_state.p_;
					Real rho_star = interface_state.rho_;

					momentum_change_rate -= 2.0 * Vol_k[index_j] *
						(SimTK::outer(rho_star * vel_star, vel_star) + p_star * Matd(1.0)) * e_ij * dW_ij;
				}
			}
			this->dmom_dt_[index_i] += momentum_change_rate;
		}
        //=================================================================================================//
        template<class BasePressureRelaxationType>
		BasePressureRelaxationWithWall<BasePressureRelaxationType>::
			BasePressureRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation) :
				BasePressureRelaxationType(fluid_wall_relation.inner_relation_,
					fluid_wall_relation.contact_relation_) {}
        //=================================================================================================//
        template<class BasePressureRelaxationType>
		BasePressureRelaxationWithWall<BasePressureRelaxationType>::
			BasePressureRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation, 
				BaseBodyRelationContact &wall_contact_relation) :
				BasePressureRelaxationType(fluid_inner_relation,
					wall_contact_relation) {}
        //=================================================================================================//
        template<class BasePressureRelaxationType>
		BasePressureRelaxationWithWall<BasePressureRelaxationType>::
			BasePressureRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation, 
				BaseBodyRelationContact &wall_contact_relation) :
				BasePressureRelaxationType(fluid_complex_relation,
					wall_contact_relation) {}
        //=================================================================================================//
        template<class BaseDensityAndenergyRelaxationType>
		template<class BaseBodyRelationType>
		DensityAndEnergyRelaxation<BaseDensityAndenergyRelaxationType>::
			DensityAndEnergyRelaxation(BaseBodyRelationType &base_body_relation,
				BaseBodyRelationContact &wall_contact_relation) :
				RelaxationWithWall<BaseDensityAndenergyRelaxationType>(base_body_relation, wall_contact_relation) {}
        //=================================================================================================//
        template<class BaseDensityAndenergyRelaxationType>
		void DensityAndEnergyRelaxation<BaseDensityAndenergyRelaxationType>::Interaction(size_t index_i, Real dt)
		{
			BaseDensityAndenergyRelaxationType::Interaction(index_i, dt);

			CompressibleFluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i], this->E_[index_i]);
			Real density_change_rate = 0.0;
			Real energy_change_rate = 0.0;
			for (size_t k = 0; k < CompressibleFluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->wall_Vol_[k]);
				StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
				Neighborhood& wall_neighborhood = (*CompressibleFluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd& e_ij = wall_neighborhood.e_ij_[n];
					Real r_ij = wall_neighborhood.r_ij_[n];
					Real dW_ij = wall_neighborhood.dW_ij_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;
					Real E_in_wall = state_i.E_;
					CompressibleFluidState state_j(rho_in_wall, vel_in_wall, p_in_wall, E_in_wall);
					CompressibleFluidState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Vecd vel_star = interface_state.vel_;
					Real p_star = interface_state.p_;
					Real rho_star = interface_state.rho_;
					Real E_star = interface_state.E_;

					density_change_rate -= 2.0 * Vol_k[index_j] * dot(rho_star * vel_star, e_ij) * dW_ij;
					energy_change_rate -= 2.0 * Vol_k[index_j] * dot(E_star * vel_star + p_star * vel_star, e_ij) * dW_ij;
				}
			}
			this->drho_dt_[index_i] += density_change_rate;
			this->dE_dt_[index_i] += energy_change_rate;
		}
        //=================================================================================================//
        template<class BaseDensityAndenergyRelaxationType>
		BaseDensityAndEnergyRelaxationWithWall<BaseDensityAndenergyRelaxationType>::
			BaseDensityAndEnergyRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation) :
			DensityAndEnergyRelaxation<BaseDensityAndenergyRelaxationType>(fluid_wall_relation.inner_relation_,
					fluid_wall_relation.contact_relation_) {}
        //=================================================================================================//
        template<class BaseDensityAndenergyRelaxationType>
		BaseDensityAndEnergyRelaxationWithWall<BaseDensityAndenergyRelaxationType>::
			BaseDensityAndEnergyRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation,
				BaseBodyRelationContact &wall_contact_relation) :
			DensityAndEnergyRelaxation<BaseDensityAndenergyRelaxationType>(fluid_inner_relation,
					wall_contact_relation) {}
        //=================================================================================================//
        template<class BaseDensityAndenergyRelaxationType>
		BaseDensityAndEnergyRelaxationWithWall<BaseDensityAndenergyRelaxationType>::
			BaseDensityAndEnergyRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation,
				BaseBodyRelationContact &wall_contact_relation) :
			DensityAndEnergyRelaxation<BaseDensityAndenergyRelaxationType>(fluid_complex_relation, wall_contact_relation) {}
		//=================================================================================================//		
    }
//=================================================================================================//
}
//=================================================================================================//