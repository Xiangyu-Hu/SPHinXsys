/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file particle_generator_mesh.h
 * @brief The mesh particle generator generates particles at the centroids
 * of mesh elements.
 * @author Zhentong Wang and Xiangyu Hu
 */

#ifndef PARTICLE_GENERATOR_MESH_H
#define PARTICLE_GENERATOR_MESH_H

#include "base_particle_generator.h"
#include "unstructured_mesh.h"

namespace SPH
{
template <> // Base class for generating particles from mesh centroids
class GeneratingMethod<UnstructuredMesh>
{
  public:
    explicit GeneratingMethod(ANSYSMesh &ansys_mesh);
    virtual ~GeneratingMethod(){};

  protected:
    StdLargeVec<Vecd> &elements_centroids_;
    StdLargeVec<Real> &elements_volumes_;
};

template <>
class ParticleGenerator<BaseParticles, UnstructuredMesh>
    : public ParticleGenerator<BaseParticles>, public GeneratingMethod<UnstructuredMesh>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, ANSYSMesh &ansys_mesh);
    virtual ~ParticleGenerator(){};
    virtual void prepareGeometricData() override;
};

} // namespace SPH
#endif // PARTICLE_GENERATOR_MESH_H
