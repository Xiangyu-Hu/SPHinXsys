/**
 * @file 	particle_dynamics_configuration.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "particle_dynamics_configuration.h"
#include "sph_system.h"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	void ParticleSortingSplitting::ConfigurationInteraction(CellList* cell_list_here, Real dt)
	{
		cout << "\n The function "
			<< "ParticleSortingSplitting::ConfigurationInteractions"
			<< " is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
	void InnerConfigurationSplitting::ConfigurationInteraction(CellList* cell_list, Real dt)
	{
		cout << "\n The function "
			<< "ParticleSortingSplitting::ConfigurationInteraction"
			<< " is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
}
//=================================================================================================//