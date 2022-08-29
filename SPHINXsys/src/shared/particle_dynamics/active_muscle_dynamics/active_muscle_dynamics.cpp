/**
 * @file 	active_muscle_dynamics.cpp
 * @brief 	In is file, we define functions declared in active muscle dynamics.h
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "active_muscle_dynamics.h"

using namespace SimTK;

namespace SPH
{
	namespace active_muscle_dynamics
	{
		//=================================================================================================//
		MuscleActivation::
			MuscleActivation(SPHBody &sph_body) :
			LocalDynamics(sph_body), ElasticSolidDataSimple(sph_body),
			pos0_(particles_->pos0_), 
			active_contraction_stress_(*particles_->getVariableByName<Real>("ActiveContractionStress")) {};
		//=================================================================================================//
    }
}
