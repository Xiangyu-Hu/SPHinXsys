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
 * @file solid_dynamics_variable.h
 * @brief Here, we define the algorithm classes for computing derived solid dynamics variables.
 * @details These variable can be added into variable list for state output.
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef SOLID_DYNAMICS_VARIABLE_H
#define SOLID_DYNAMICS_VARIABLE_H

#include "base_general_dynamics.h"
#include "elastic_solid.h"

namespace SPH
{
/**
 * @class Displacement
 * @brief computing displacement from current and initial particle position
 */
class Displacement : public BaseDerivedVariable<Vecd>
{
  public:
    explicit Displacement(SPHBody &sph_body);
    virtual ~Displacement(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_, *pos0_;
};

/**
 * @class OffsetInitialPosition
 * @brief offset initial particle position
 */
class OffsetInitialPosition : public LocalDynamics
{
  public:
    explicit OffsetInitialPosition(SPHBody &sph_body, Vecd &offset);
    virtual ~OffsetInitialPosition(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd offset_;
    Vecd *pos_, *pos0_;
};

/**
 * @class TranslationAndRotation
 * @brief transformation on particle position and rotation
 */
class TranslationAndRotation : public LocalDynamics
{
  public:
    explicit TranslationAndRotation(SPHBody &sph_body, Transform &transform);
    virtual ~TranslationAndRotation(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Transform &transform_;
    Vecd *pos_, *pos0_;
};

class GreenLagrangeStrain : public BaseDerivedVariable<Matd>
{
  public:
    explicit GreenLagrangeStrain(SPHBody &sph_body);
    virtual ~GreenLagrangeStrain(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *F_;
};

/**
 * @class VonMisesStress
 * @brief computing von_Mises_stress
 */
class VonMisesStress : public BaseDerivedVariable<Real>
{
  public:
    explicit VonMisesStress(SPHBody &sph_body);
    virtual ~VonMisesStress(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real rho0_;
    Real *rho_;
    Matd *F_;
    ElasticSolid &elastic_solid_;
};

/**
 * @class VonMisesStrain
 * @brief computing von Mises strain
 */
class VonMisesStrain : public BaseDerivedVariable<Real>
{
  public:
    explicit VonMisesStrain(SPHBody &sph_body);
    virtual ~VonMisesStrain(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *F_;
};

/**
 * @class VonMisesStrain
 * @brief update von Mises strain
 */
class VonMisesStrainDynamic : public BaseDerivedVariable<Real>
{
  public:
    explicit VonMisesStrainDynamic(SPHBody &sph_body);
    virtual ~VonMisesStrainDynamic(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    ElasticSolid &elastic_solid_;
    Real poisson_ratio_;
    Matd *F_;
};

/**
 * @class MidSurfaceVonMisesStress
 * @brief computing mid-surface von Mises stress of shells
 */
class MidSurfaceVonMisesStress : public BaseDerivedVariable<Real>
{
  public:
    explicit MidSurfaceVonMisesStress(SPHBody &sph_body);
    virtual ~MidSurfaceVonMisesStress(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *mid_surface_cauchy_stress_;
};
} // namespace SPH
#endif // SOLID_DYNAMICS_VARIABLE_H