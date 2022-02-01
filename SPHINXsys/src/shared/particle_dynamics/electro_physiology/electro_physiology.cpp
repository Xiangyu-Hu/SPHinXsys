/**
 * @file 	electro_physiology.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 */

#include "electro_physiology.h"

using namespace SimTK;

namespace SPH
{
	namespace electro_physiology
	{
		//=================================================================================================//
		ElectroPhysiologyInitialCondition::
			ElectroPhysiologyInitialCondition(SolidBody &sph_body) :
			ParticleDynamicsSimple(sph_body),
			ElectroPhysiologyDataDelegateSimple(sph_body),
			pos_n_(particles_->pos_n_), species_n_(particles_->species_n_)
		{
		}
		//=================================================================================================//
	}
}