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
 * @file 	complex_solid.h
 * @brief 	These are classes for define complex solid materials.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef COMPLEX_SOLID_H
#define COMPLEX_SOLID_H

#include "elastic_solid.h"
#include <general_dynamics.h>

namespace SPH
{
/**
 * @class ActiveMuscle
 * @brief Here, the active response is considered.
 */
template <class MuscleType>
class ActiveMuscle : public MuscleType
{
  protected:
    StdLargeVec<Real> active_contraction_stress_; /**<  active contraction stress */

  public:
    template <typename... ConstructorArgs>
    explicit ActiveMuscle(ConstructorArgs &&...args);
    virtual ~ActiveMuscle(){};

    /** initialize the local properties, fiber and sheet direction. */
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    /** compute the stress through Constitutive relation. */
    virtual Matd StressPK2(Matd &deformation, size_t index_i) override;
    virtual ActiveMuscle<MuscleType> *ThisObjectPtr() override { return this; };
};

/**
 * @class CompositeSolid
 */
class CompositeSolid : public ElasticSolid
{
    StdLargeVec<int> material_id_;

  protected:
    UniquePtrsKeeper<ElasticSolid> composite_ptrs_keeper_;
    StdVec<ElasticSolid *> composite_materials_;

  public:
    explicit CompositeSolid(Real rho0);
    virtual ~CompositeSolid(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Matd StressPK2(Matd &deformation, size_t index_i) override;
    virtual Matd StressPK1(Matd &deformation, size_t index_i) override;
    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t index_i) override;
    virtual Real VolumetricKirchhoff(Real J) override { return 0.0; };
    virtual std::string getRelevantStressMeasureName() override { return "PK2"; };

    Real CompositeDensity(size_t index_i);

    template <class ElasticSolidType, typename... Args>
    void add(Args &&...args)
    {
        ElasticSolid *added_material =
            composite_ptrs_keeper_.createPtr<ElasticSolidType>(std::forward<Args>(args)...);
        composite_materials_.push_back(added_material);
        c0_ = SMAX(c0_, added_material->ReferenceSoundSpeed());
        setContactStiffness(c0_);
    };
};

/**
 * @class MaterialIdInitialization
 */
class MaterialIdInitialization
    : public LocalDynamics,
      public GeneralDataDelegateSimple
{
  public:
    explicit MaterialIdInitialization(SPHBody &sph_body);

  protected:
    StdLargeVec<int> &material_id_;
    StdLargeVec<Vecd> &pos0_;
};
} // namespace SPH
#endif // COMPLEX_SOLID_H
