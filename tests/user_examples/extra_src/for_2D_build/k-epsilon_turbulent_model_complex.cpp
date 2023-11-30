#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
////=================================================================================================//
//		StandardWallFunctionCorrection::
//			StandardWallFunctionCorrection(ComplexRelation& complex_relation, Real offset_dist, const StdVec<int>& id_exclude, NearShapeSurface& near_surface)
//			: StandardWallFunctionCorrection(complex_relation.,
//				complex_relation.getContactRelation(), offset_dist, id_exclude, near_surface) {}
		//=================================================================================================//
		StandardWallFunctionCorrection::
			StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
				BaseContactRelation& contact_relation, Real offset_dist, const StdVec<int>& id_exclude, NearShapeSurface& near_surface)
			: LocalDynamics(inner_relation.getSPHBody()), FSIContactData(contact_relation), id_exclude_(id_exclude), 
			level_set_shape_(&near_surface.level_set_shape_),
			offset_dist_(offset_dist),vel_(particles_->vel_),pos_(particles_->pos_),dimension_(Vecd(0).size()),
			rho_(particles_->rho_), mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()) ,
			particle_spacing_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
			cutoff_radius_(inner_relation.getSPHBody().sph_adaptation_->getKernel()->CutOffRadius()),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			is_near_wall_P1_(*particles_->getVariableByName<int>("IsNearWallP1")),
			velocity_gradient_(*particles_->getVariableByName<Matd>("VelocityGradient")),
			k_production_(*particles_->getVariableByName<Real>("K_Production"))
		{
			
			particles_->registerVariable(dist_to_dmy_itfc_ls_, "DistanceToDummyInterfaceLS");
			particles_->registerSortableVariable<Real>("DistanceToDummyInterfaceLS");
			particles_->addVariableToWrite<Real>("DistanceToDummyInterfaceLS");

			particles_->registerVariable(dist_to_dmy_itfc_aver_, "DistanceToDummyInterfaceAver");
			particles_->registerSortableVariable<Real>("DistanceToDummyInterfaceAver");
			particles_->addVariableToWrite<Real>("DistanceToDummyInterfaceAver");

			particles_->registerVariable(y_p_, "Y_P");
			particles_->registerSortableVariable<Real>("Y_P");
			particles_->addVariableToWrite<Real>("Y_P");

			particles_->registerVariable(wall_Y_plus_, "WallYplus");
			particles_->registerSortableVariable<Real>("WallYplus");
			particles_->addVariableToWrite<Real>("WallYplus");

			particles_->registerVariable(wall_Y_star_, "WallYstar");
			particles_->registerSortableVariable<Real>("WallYstar");
			particles_->addVariableToWrite<Real>("WallYstar");

			particles_->registerVariable(is_near_wall_P2_, "IsNearWallP2");
			particles_->registerSortableVariable<int>("IsNearWallP2");
			particles_->addVariableToWrite<int>("IsNearWallP2");
			particles_->registerVariable(is_near_wall_P1_pre_, "IsNearWallP1Pre");
			particles_->registerSortableVariable<int>("IsNearWallP1Pre");
			particles_->addVariableToWrite<int>("IsNearWallP1Pre");
			particles_->registerVariable(is_migrate_, "IsMigrate");
			particles_->registerSortableVariable<int>("IsMigrate");
			particles_->addVariableToWrite<int>("IsMigrate");

			particles_->registerVariable(velo_friction_, "FrictionVelocity");
			particles_->registerSortableVariable<Vecd>("FrictionVelocity");
			particles_->addVariableToWrite<Vecd>("FrictionVelocity");

			particles_->registerVariable(velo_tan_, "TangentialVelocity");
			particles_->registerSortableVariable<Real>("TangentialVelocity");
			particles_->addVariableToWrite<Real>("TangentialVelocity");

			particles_->registerVariable(index_nearest, "NearestIndex");
			particles_->registerSortableVariable<int>("NearestIndex");
			particles_->addVariableToWrite<int>("NearestIndex");

			particles_->registerVariable(dist_to_dmy_interface_, "DistanceToDummyInterface");
			particles_->registerSortableVariable<Real>("DistanceToDummyInterface");
			particles_->addVariableToWrite<Real>("DistanceToDummyInterface");

			/*definition of near wall particles*/
			intial_distance_to_wall = 1.5 * particle_spacing_; //changed
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_n_.push_back(&(contact_particles_[k]->n_));
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
			/*definition the specific particles where normal vectors do not exist*/
			std::cout << "The id of excluded wall particle is: " << std::endl;
			for (int i = 0; i != id_exclude_.size(); ++i)
			{
				std::cout << id_exclude_[i] << std::endl;
			}
			//std::cout << "Are you sure? " << std::endl;
			//system("pause");
		};
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//