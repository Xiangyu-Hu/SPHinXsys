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
 * @file 	secondary_structure.h
 * @brief 	A complex body is characterized with a secondary structure,
 * 			which can be imported externally or created according to specific rules.
 * 			The secondary structure will be used or even created by the corresponding
 * 			particle generator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SECONDARY_STRUCTURE_H
#define SECONDARY_STRUCTURE_H

#include "base_body.h"
#include "base_body_part.h"
#include "base_geometry.h"

namespace SPH
{
/**
 * @class SecondaryStructure
 * @brief Abstract class as interface for all secondary structures.
 * Currently, it provides interface on building inner configuration.
 * The interface can be extended.
 */
class SecondaryStructure
{
  public:
    explicit SecondaryStructure(){};
    virtual ~SecondaryStructure(){};

    virtual void buildParticleConfiguration(ParticleConfiguration &particle_configuration) = 0;
};
} // namespace SPH
#endif // SECONDARY_STRUCTURE_H