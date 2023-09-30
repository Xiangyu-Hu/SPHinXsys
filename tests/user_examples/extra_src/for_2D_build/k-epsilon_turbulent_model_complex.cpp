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
			StandardWallFunctionCorrection(ComplexRelation& complex_relation, Real y_p_theo)
			: StandardWallFunctionCorrection(complex_relation.getInnerRelation(),
				complex_relation.getContactRelation(), y_p_theo) {}
		//=================================================================================================//
		StandardWallFunctionCorrection::
			StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
				BaseContactRelation& contact_relation, Real y_p_theo)
			: LocalDynamics(inner_relation.getSPHBody()), FSIContactData(contact_relation),
			y_p_theo_(y_p_theo),vel_(particles_->vel_),pos_(particles_->pos_),dimension_(Vecd(0).size()),
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
			//** Calculate the offset distance for wall function wall *
			offset_dist_ = y_p_theo_ - 0.5 * particle_spacing_;

			particles_->registerVariable(y_p_, "Y_P");
			particles_->registerSortableVariable<Real>("Y_P");
			particles_->addVariableToWrite<Real>("Y_P");

			particles_->registerVariable(wall_Y_plus_, "WallYplus");
			particles_->registerSortableVariable<Real>("WallYplus");
			particles_->addVariableToWrite<Real>("WallYplus");

			particles_->registerVariable(wall_Y_star_, "WallYstar");
			particles_->registerSortableVariable<Real>("WallYstar");
			particles_->addVariableToWrite<Real>("WallYstar");


			//particles_->registerVariable(is_near_wall_P1_, "IsNearWallP1");
			//particles_->registerSortableVariable<int>("IsNearWallP1");
			//particles_->addVariableToWrite<int>("IsNearWallP1");
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

			particles_->registerVariable(distance_to_wall_, "DistanceToWall");
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
		Real StandardWallFunctionCorrection::getFrictionVelo(Real left_bound, Real right_bound, Real e, Real A, Real B)
		{
			Real EPSILON = e;
			Real left_value = WallFunc(left_bound,A,B);
			Real right_value = WallFunc(right_bound, A, B);

			if (left_value * right_value >= 0)
			{
				std::cout << "You have not assumed right a and b\n";
				system("pause");
				return 0.0;
				exit(1);
			}
			Real middle = left_bound;
			Real middle_value = 0.0;
			while ((right_bound - left_bound) >= EPSILON)
			{
				// Update middle point
				middle = (left_bound + right_bound) / 2;
				middle_value = WallFunc(middle, A, B);
				// Check if middle point is root
				if (abs(middle_value) < EPSILON)
					break;
				// Decide the side to repeat the steps
				else if (middle_value * left_value < 0)
					right_bound = middle;
				else
					left_bound = middle;
			}
			//cout << "The value of root is : " << c;
			return middle;
		}
		Real StandardWallFunctionCorrection::WallFunc(Real x, Real CA, Real CB)
		{
			return CA / x - log(CB * x);
		}
		//=================================================================================================//
		void StandardWallFunctionCorrection::checkFrictionVelo(Real velo_fric, Real e, Real A, Real B)
		{
			Real EPSILON = e;
			Real left_value = A / velo_fric;
			Real right_value = log(B * velo_fric);
			Real error = abs(left_value - right_value);
			if (error > EPSILON && GlobalStaticVariables::physical_time_>2.0)
			//if (error > EPSILON )
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