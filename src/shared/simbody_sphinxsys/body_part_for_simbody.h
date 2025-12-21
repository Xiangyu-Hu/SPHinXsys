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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    body_part_for_simbody.h
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Xiangyu Hu
 */

#ifndef BODY_PART_FOR_SIMBODY_H
#define BODY_PART_FOR_SIMBODY_H

#include "base_body_part.h"

#include "simtk_wrapper.h"
namespace SPH
{
/**
 * @class SolidBodyPartForSimbody
 * @brief A SolidBodyPart for coupling with Simbody.
 * The mass, origin, and unit inertial matrix are computed.
 * Note: In Simbody, all spatial vectors are three dimensional.
 */
class SolidBodyPartForSimbody : public BodyRegionByParticle
{
  protected:
    UniquePtrKeeper<SimTK::MassProperties> mass_properties_keeper_;

  public:
    Vecd initial_mass_center_;
    SimTK::MassProperties *body_part_mass_properties_;

    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape);
    SolidBodyPartForSimbody(SPHBody &body, SharedPtr<Shape> shape_ptr);
    virtual ~SolidBodyPartForSimbody() {};

  protected:
    Real rho0_;
    Real *Vol_;
    Vecd *pos_;

  private:
    void setMassProperties();
};
} // namespace SPH
#endif // BODY_PART_FOR_SIMBODY_H
