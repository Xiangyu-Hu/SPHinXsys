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
/**
 * @file 	net_work.cpp
 * @brief 	This is the example of generating a neural network on a sphere 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */
/** header file and namespace. */
#include "sphinxsys.h"
/** case file to setup the test case */
#include "sphere.h"

using namespace SPH;
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(domain_lower_bound, domain_upper_bound, dp_0);
	/** Output */
	In_Output in_output(system);
	/** Creat a sphere, corresponding material and particles. */
	MyPolygonBody *polygon_body = new MyPolygonBody(system, "Polygon", 0, 
		new ParticleGeneratorNetwork(Point(-1.0, 0.0, 0.0), Point(-0.964, 0.0, 0.266)));
	BodyMaterial* body_material = new BodyMaterial();
	ElasticSolidParticles 	body_particles(polygon_body, body_material);
	/** Write particle data. */
	WriteBodyStatesToVtu 		write_states(in_output, system.real_bodies_);
	write_states.WriteToFile(0.0);
	return 0;
}