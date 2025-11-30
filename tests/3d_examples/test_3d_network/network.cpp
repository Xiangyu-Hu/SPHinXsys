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
 * @file network.cpp
 * @brief This is the example of generating a neural network on a sphere
 * TODO: particles here better to linear particles.
 * However, the current implementation of linear particles is not consistent
 * with the particles used here.
 * @author Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;

Vec3d domain_lower_bound(-1.0, -1.0, -1.0);
Vec3d domain_upper_bound(1.0, 1.0, 1.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
/** Network starting point. */
Vecd starting_point(-1.0, 0.0, 0.0);
/** Network second point. */
Vecd second_point(-0.964, 0.0, 0.266);
/** Network iterative levels. */
int iteration_levels = 15;
/** Network defecting angle. */
Real grad_factor = 5.0;

int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    /** Creat a body, corresponding material and particles. */
    TreeBody tree_on_sphere(sph_system, makeShared<GeometricShapeBall>(Vec3d::Zero(), 1.0, "Sphere"));
    tree_on_sphere.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    tree_on_sphere.generateParticles<BaseParticles, Network>(starting_point, second_point, iteration_levels, grad_factor);
    /** Write particle data. */
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.writeToFile(0);

    return 0;
}