#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
//=================================================================================================//
		StandardWallFunctionCorrection::
			StandardWallFunctionCorrection(ComplexRelation& complex_relation)
			: StandardWallFunctionCorrection(complex_relation.getInnerRelation(),
				complex_relation.getContactRelation()) {}
		//=================================================================================================//
		StandardWallFunctionCorrection::
			StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
				BaseContactRelation& contact_relation)
			: LocalDynamics(inner_relation.getSPHBody()),FSIContactData(contact_relation), 
			vel_(particles_->vel_),pos_(particles_->pos_),dimension_(Vecd(0).size()), 
			rho_(particles_->rho_),mu_(particles_->fluid_.ReferenceViscosity()),
			particle_spacing_(inner_relation.getSPHBody().sph_adaptation_->ReferenceSpacing()),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity"))
		{
			particles_->registerVariable(wall_Y_plus_, "WallYplus");
			particles_->registerSortableVariable<Real>("WallYplus");
			particles_->addVariableToWrite<Real>("WallYplus");

			particles_->registerVariable(wall_Y_star_, "WallYstar");
			particles_->registerSortableVariable<Real>("WallYstar");
			particles_->addVariableToWrite<Real>("WallYstar");


			particles_->registerVariable(is_near_wall_P1_, "IsNearWallP1");
			particles_->registerSortableVariable<int>("IsNearWallP1");
			particles_->addVariableToWrite<int>("IsNearWallP1");
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
			particles_->registerSortableVariable<Real>("FrictionVelocity");
			particles_->addVariableToWrite<Real>("FrictionVelocity");

			particles_->registerVariable(velo_tan_, "TangentialVelocity");
			particles_->registerSortableVariable<Real>("TangentialVelocity");
			particles_->addVariableToWrite<Real>("TangentialVelocity");

			particles_->registerVariable(index_nearest, "NearestIndex");
			particles_->registerSortableVariable<int>("NearestIndex");
			particles_->addVariableToWrite<int>("NearestIndex");

			particles_->registerVariable(distance_to_wall, "DistanceToWall");
			particles_->registerSortableVariable<Real>("DistanceToWall");
			particles_->addVariableToWrite<Real>("DistanceToWall");

			//definition of near wall particles
			intial_distance_to_wall = 1.5 * particle_spacing_; //changed
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_n_.push_back(&(contact_particles_[k]->n_));
			}
		};

		//=================================================================================================//
		Real StandardWallFunctionCorrection::getFrictionVelo(Real a, Real b, Real e)
		{
			Real EPSILON = e;
			if (WallFunc(a) * WallFunc(b) >= 0)
			{
				std::cout << "You have not assumed right a and b\n";
				system("pause");
				return 0.0;
				exit(1);
			}
			Real c = a;
			while ((b - a) >= EPSILON)
			{
				// Find middle point
				c = (a + b) / 2;
				// Check if middle point is root
				if (WallFunc(c) == 0.0)
					break;
				// Decide the side to repeat the steps
				else if (WallFunc(c) * WallFunc(a) < 0)
					b = c;
				else
					a = c;
			}
			//cout << "The value of root is : " << c;
			return c;
		}
		//=================================================================================================//
		void StandardWallFunctionCorrection::checkFrictionVelo(Real velo_fric, Real e)
		{
			Real EPSILON = e;
			Real left_value = coefficientA / velo_fric;
			Real right_value = log(coefficientB * velo_fric);
			Real error = abs(left_value - right_value);
			if (error > EPSILON)
			{
				std::cout << "FrictionVelocity is not accurately calcualted" << std::endl;
				std::cout << "error=" << error << std::endl;
				std::cout << "velo_friction=" << velo_fric << std::endl;
				//system("pause");
			}
		}
		//=================================================================================================//

	}
	//=================================================================================================//
}
//=================================================================================================//