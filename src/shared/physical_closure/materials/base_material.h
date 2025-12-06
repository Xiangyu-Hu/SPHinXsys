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
 * @file 	base_material.h
 * @brief 	This is the base classes of all materials.
 *		    A function in a derived material class returns a value with the inputs
 *          from the particle data.
 *			Basically, it is a interface from which
 *			one can access derived material by dynamic cast.
 *          Note that the derived material may have position dependent or
 *          local properties.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_MATERIAL_H
#define BASE_MATERIAL_H

#include "base_data_type_package.h"
#include "sphinxsys_containers.h"

namespace SPH
{
class BaseParticles;
class SPHSystem;

/** @class  BaseMaterial
 *  @brief Base of all materials
 *  @details Note that the case dependent parameters of the material properties
 *  will be defined in applications.
 */
class BaseMaterial
{
  public:
    explicit BaseMaterial(Real rho0)
        : material_type_name_("BaseMaterial"), rho0_(rho0) {};
    BaseMaterial() : BaseMaterial(1.0) {};
    virtual ~BaseMaterial() {};
    std::string MaterialType() { return material_type_name_; }
    Real ReferenceDensity() { return rho0_; };
    /**
     * This will be called after particles generation
     * and is important because particles are not defined yet when material is constructed.
     * For a composite material, i.e. there is a material pointer with another material,
     * one need assign the base particle to that material too.
     */
    void setLocalParameters(SPHSystem &sph_system, BaseParticles *base_particles);
    virtual void registerLocalParameters(BaseParticles *base_particles) {};
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) {};
    virtual void initializeLocalParameters(BaseParticles *base_particles) {};

  protected:
    std::string material_type_name_;
    Real rho0_; /**< reference density. */
};

/** @class  Fluid
 *  @brief  Base class of all fluids
 */
class Fluid : public BaseMaterial
{
  protected:
    Real c0_; /**< reference sound speed and pressure */

  public:
    explicit Fluid(Real rho0, Real c0);
    virtual ~Fluid() {};
    Real ReferenceSoundSpeed() { return c0_; };
    virtual Real getPressure(Real rho) = 0;
    virtual Real getPressure(Real rho, Real rho_e) { return getPressure(rho); };
    virtual Real DensityFromPressure(Real p) = 0;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;

    class EosKernel
    {
      public:
        EosKernel(Fluid &encloser);

      protected:
        Real c0_, rho0_;
    };
};

/** @class  Solid
 *  @brief Base class of all solid materials
 */
class Solid : public BaseMaterial
{
  public:
    Solid(Real rho0, Real contact_stiffness, Real contact_friction = 0.0)
        : BaseMaterial(rho0), contact_stiffness_(contact_stiffness),
          contact_friction_(contact_friction)
    {
        material_type_name_ = "Solid";
    };
    explicit Solid(Real rho0) : Solid(rho0, 1.0) {};
    Solid() : Solid(1.0) {};
    virtual ~Solid() {};

    Real ContactFriction() { return contact_friction_; };
    Real ContactStiffness() { return contact_stiffness_; };
    /** Get average velocity when interacting with fluid. */
    virtual Vecd *AverageVelocity(BaseParticles *base_particles);
    /** Get average acceleration when interacting with fluid. */
    virtual Vecd *AverageAcceleration(BaseParticles *base_particles);

    /** Get average velocity when interacting with fluid. */
    virtual DiscreteVariable<Vecd> *AverageVelocityVariable(BaseParticles *base_particles);
    /** Get average acceleration when interacting with fluid. */
    virtual DiscreteVariable<Vecd> *AverageAccelerationVariable(BaseParticles *base_particles);

  protected:
    Real contact_stiffness_; /**< contact-force stiffness related to bulk modulus*/
    Real contact_friction_;  /**< friction property mimic fluid viscosity*/
    void setContactStiffness(Real c0) { contact_stiffness_ = rho0_ * c0 * c0; };
};
} // namespace SPH
#endif // BASE_MATERIAL_H