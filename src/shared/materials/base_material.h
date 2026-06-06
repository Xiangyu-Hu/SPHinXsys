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

#include "data_type.h"
#include "sphinxsys_variable.h"

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
    BaseMaterial() : material_type_name_("BaseMaterial"){};
    virtual ~BaseMaterial();
    std::string Name();
    void setName(const std::string &name) { material_name_ = name; };
    std::string MaterialType() { return material_type_name_; }
    /**
     * This will be called after particles generation
     * and is important because particles are not defined yet when material is constructed.
     * For a composite material, i.e. there is a material pointer with another material,
     * one need assign the base particle to that material too.
     */
    void setLocalParameters(SPHSystem &sph_system, BaseParticles *base_particles);
    virtual void registerLocalParametersToReload(BaseParticles *base_particles) {};
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) {};
    virtual void initializeLocalParameters(BaseParticles *base_particles) {};

  protected:
    std::string material_type_name_;
    std::string material_name_;
    UniquePtrsKeeper<Quantity> unique_entity_ptrs_;
};

class MatterMaterial : public BaseMaterial
{
    DiscreteVariable<Real> *dv_rho_;
    DiscreteVariable<Real> *dv_mass_;

  public:
    explicit MatterMaterial();
    virtual ~MatterMaterial(){};
    DiscreteVariable<Real> *dvDensity() const { return dv_rho_; };
    DiscreteVariable<Real> *dvMass() const { return dv_mass_; };
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Real ReferenceDensity() const = 0;
};

/** @class  Fluid
 *  @brief  Base class of all fluids
 */
class Fluid : public MatterMaterial
{
  public:
    explicit Fluid();
    virtual ~Fluid(){};
    virtual Real getPressure(Real rho) = 0;
    virtual Real getPressure(Real rho, Real rho_e) { return getPressure(rho); };
    virtual Real DensityFromPressure(Real p) = 0;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;
    virtual Real ReferenceSoundSpeed() const = 0;
};

class SolidContact
{
  private:
    Real rho0_copy_; /**< reference density. */
  public:
    explicit SolidContact(Real rho0, Real contact_stiffness, Real contact_friction = 0.0);
    virtual ~SolidContact(){};
    Real ContactReferenceDensity() { return rho0_copy_; };
    Real ContactFriction() { return contact_friction_; };
    Real ContactStiffness() { return contact_stiffness_; };

  protected:
    Real contact_stiffness_; /**< contact-force stiffness related to bulk modulus*/
    Real contact_friction_;  /**< friction property mimic fluid viscosity*/
    void setContactStiffness(Real rho0, Real c0) { contact_stiffness_ = rho0 * c0 * c0; };
};

/** @class  Solid
 *  @brief Base class of all solid materials
 */
class Solid : public MatterMaterial, public SolidContact
{
  public:
    Solid(Real rho0, Real contact_stiffness, Real contact_friction = 0.0);
    explicit Solid(Real rho0) : Solid(rho0, 1.0){};
    Solid() : Solid(1.0){};
    virtual ~Solid(){};
    virtual Real ReferenceDensity() const override { return rho0_; };
    /** Get average velocity when interacting with fluid. */
    virtual Vecd *AverageVelocity(BaseParticles *base_particles);
    /** Get average acceleration when interacting with fluid. */
    virtual Vecd *AverageAcceleration(BaseParticles *base_particles);

    /** Get average velocity when interacting with fluid. */
    virtual DiscreteVariable<Vecd> *AverageVelocityVariable(BaseParticles *base_particles);
    /** Get average acceleration when interacting with fluid. */
    virtual DiscreteVariable<Vecd> *AverageAccelerationVariable(BaseParticles *base_particles);

  protected:
    Real rho0_; /**< reference density. */
};
} // namespace SPH
#endif // BASE_MATERIAL_H