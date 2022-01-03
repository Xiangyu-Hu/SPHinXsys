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

#include "base_data_package.h"
#include "xml_engine.h"

#include <string>

namespace SPH
{
	class BaseParticles;
	class FluidParticles;
	class SolidParticles;

	/** @class  BaseMaterial
	 *  @brief Base of all materials
	 *  @details Note that the case dependent parameters of the material properties
	 *  will be defined in applications.
	 */
	class BaseMaterial
	{
	protected:
		std::string material_type_;
		std::string parameters_name_;
		Real rho0_; /**< reference density. */
		BaseParticles *base_particles_;
		XmlEngine reload_material_xml_engine_;
		ParticleVariableList reload_local_parameters_;

	public:
		explicit BaseMaterial(Real rho0)
			: material_type_("BaseMaterial"),
			  parameters_name_("LocalMaterialParameters"),
			  rho0_(rho0), base_particles_(nullptr),
			  reload_material_xml_engine_("xml_material", "local_material_paramaters"){};
		BaseMaterial() : BaseMaterial(1.0){};
		virtual ~BaseMaterial(){};

		/** This will be called in BaseParticle constructor
		 * and is important because particles are not defined in SPHBody constructor.
		 * For a composite material, i.e. there is a material pointer with another material,
		 * one need assign the base particle to that material too. */
		virtual void assignBaseParticles(BaseParticles *base_particles);
		std::string MaterialType() { return material_type_; }
		std::string LocalParametersName() { return parameters_name_; }
		Real ReferenceDensity() { return rho0_; };

		virtual void writeToXmlForReloadLocalParameters(const std::string &filefullpath);
		virtual void readFromXmlForLocalParameters(const std::string &filefullpath);

		virtual BaseMaterial *ThisObjectPtr() { return this; };
	};

	/** @class  Fluid
	 *  @brief Base class of all fluids
	*/
	class Fluid : public BaseMaterial
	{
	protected:
		Real mu_; /**< reference sound speed, viscosity. */
		FluidParticles *fluid_particles_;

	public:
		explicit Fluid(Real rho0, Real mu)
			: BaseMaterial(rho0), mu_(mu), fluid_particles_(nullptr)
		{
			material_type_ = "Fluid";
		};
		virtual ~Fluid(){};

		void assignFluidParticles(FluidParticles *fluid_particles)
		{
			fluid_particles_ = fluid_particles;
		};

		Real ReferenceViscosity() { return mu_; };
		virtual Real getPressure(Real rho) = 0;
		virtual Real getPressure(Real rho, Real rho_e) { return getPressure(rho); };
		virtual Real DensityFromPressure(Real p) = 0;
		virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;
		virtual Fluid *ThisObjectPtr() override { return this; };
	};

	/** @class  Solid
	 *  @brief Base class of all solid materials
	*/
	class Solid : public BaseMaterial
	{
	public:
		Solid(Real rho0, Real contact_stiffness, Real contact_friction = 0.0)
		: BaseMaterial(rho0), contact_stiffness_(contact_stiffness),
			  contact_friction_(contact_friction), solid_particles_(nullptr)
		{
			material_type_ = "Solid";
		};
		explicit Solid(Real rho0) : Solid(rho0, 1.0) {};
		Solid() : Solid(1.0) {}; 
		virtual ~Solid(){};

		void assignSolidParticles(SolidParticles *solid_particles)
		{
			solid_particles_ = solid_particles;
		};

		Real ContactFriction() { return contact_friction_; };
		Real ContactStiffness() { return contact_stiffness_; };
		virtual Solid *ThisObjectPtr() override { return this; };

	protected:
		Real contact_stiffness_; /**< contact-force stiffness related to bulk modulus*/
		Real contact_friction_;	 /**< friction property mimic fluid viscosity*/
		SolidParticles *solid_particles_;

		void setContactStiffness(Real c0) { contact_stiffness_ = c0 * c0; };
	};
}
#endif //BASE_MATERIAL_H