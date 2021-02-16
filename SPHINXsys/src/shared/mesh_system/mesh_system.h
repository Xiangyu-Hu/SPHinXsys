/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
#pragma once

#include "base_data_package.h"
#include "all_kernels.h"
#include "all_types_of_bodies.h"
#include "particle_generator.h"
#include "sph_system.h"
//#include "output.h"

namespace SPH
{

	class SPHSystem;

	/*
	Mesh system also managing objects in the system level
	It is dual system respected to SPH system
	Here a base class is defined for 2d and 3d 
	which will defined as dervied classes
	*/

	class MeshSystem
	{
		//creat a object for output
		virtual void create_output() = 0;
		//delete the object for output field data
		virtual void destroy_output() = 0;
		//domain bounds
		Vecd lower_bound_, upper_bound_;
		//reference mesh spacing
		Real mesh_spacing_;

	public:

		//SPHSystem constructor
		MeshSystem();
		virtual ~MeshSystem();

		//whether a domain added
		bool is_domain_added_;
		//Initialize TBB
		//task_scheduler_init *tbb_init_;

		//the bodies in the SPH system
		StdVec<SPHBody* > bodies_;
		//the real bodies
		StdVec<SPHBody* > real_bodies_;
		//the ghost bodies
		StdVec<SPHBody* > ghost_bodies_;

		//output object
		Output* output_;

		//public functions to build up the system
		//add the domain, which managing the geometric information 
		//add a SPH body in the system
		void addABody(SPHBody* body);
		//add a particle configuartion to the system 
		void AddNNPS(NNPS* nnps = nullptr, bool Eulerian = true);
		//write a output file at a given time
		void OutputFile(Real _time);
		//initialize the system
		void InitializeSystem();
	};
}