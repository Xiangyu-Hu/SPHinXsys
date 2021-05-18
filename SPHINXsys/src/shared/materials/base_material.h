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
#include "base_particles.h"

#include <string>

namespace SPH {

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
		std::string material_name_;
		std::string parameters_name_;
		Real rho_0_; /**< reference density. */
		BaseParticles* base_particles_;
		XmlEngine reload_material_xml_engine_;
		//----------------------------------------------------------------------
		//		Generalized data for parameters
		//----------------------------------------------------------------------
		ParticleData	parameter_data_;
		ParticleDataMap parameter_maps_;

		virtual void assignDerivedMaterialParameters() {};
	public:
		BaseMaterial() : material_name_("BaseMaterial"), parameters_name_("LocalParameters"),
		rho_0_(1.0), base_particles_(NULL),
		reload_material_xml_engine_("xml_material", "material_paramaters") {};
		virtual ~BaseMaterial() {};

		/** This will be called in BaseParticle constructor
		 * and is important because particles are not defined in SPHBody constructor.  
		 * For a composite material, i.e. there is a material pointer with another material,
		 * one need assign the base particle to that material too. */
		virtual void assignBaseParticles(BaseParticles* base_particles) 
		{ 
			base_particles_ = base_particles; 
		};
		std::string MaterialName() { return material_name_;}
		std::string LocalParametersName() { return parameters_name_;}
		Real ReferenceDensity() { return rho_0_; };

		virtual void writeToXmlForReloadLocalParameters(std::string &filefullpath)
		{
			std::cout << "\n Material properties writing. " << std::endl;
			size_t total_real_particles = base_particles_->total_real_particles_;
			reload_material_xml_engine_.resizeXmlDocForParticles(total_real_particles);
			WriteAParticleVariableToXml write_variable_to_xml(reload_material_xml_engine_);
			loopParticleData<loopParticleDataMap>(parameter_data_, parameter_maps_, write_variable_to_xml);
			reload_material_xml_engine_.writeToXmlFile(filefullpath);
			std::cout << "\n Material properties writing finished. " << std::endl;
		};
		virtual void readFromXmlForLocalParameters(std::string &filefullpath)
		{
			reload_material_xml_engine_.loadXmlFile(filefullpath);
			ReadAParticleVariableFromXml read_variable_from_xml(reload_material_xml_engine_);
			loopParticleData<loopParticleDataMap>(parameter_data_, parameter_maps_, read_variable_from_xml);

			size_t total_real_particles = base_particles_->total_real_particles_;
			if(total_real_particles != reload_material_xml_engine_.SizeOfXmlDoc())
			{
				std::cout << "\n Error: reload material properties does not match!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}	
			else 
			{
				std::cout << "\n Material properties reading finished." << std::endl;
			}		
		};
		
		virtual void writeLocalParametersToVtuFile(std::ofstream& output_file) 
		{
			size_t total_real_particles = base_particles_->total_real_particles_;
			//write vectors
			for (auto const& name_index : parameter_maps_[indexVector])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<Vecd>& variable = *(std::get<indexVector>(parameter_data_)[name_index.second]);
				output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
				output_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i) {
					Vec3d vector_value = upgradeToVector3D(variable[i]);
					output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
				}
				output_file << std::endl;
				output_file << "    </DataArray>\n";
			}
		
			//write scalars
			for (auto const& name_index : parameter_maps_[indexScalar])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<Real>& variable = *(std::get<indexScalar>(parameter_data_)[name_index.second]);
				output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
				output_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i) {
					output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
				}
				output_file << std::endl;
				output_file << "    </DataArray>\n";
			}			
		};
		
		virtual BaseMaterial* ThisObjectPtr() { return this; };
	};

	/** @class  Fluid
	 *  @brief Base class of all fluids
	*/
	class Fluid : public BaseMaterial
	{
	protected:
		Real c_0_, mu_; /**< reference sound speed, viscosity. */
		FluidParticles* fluid_particles_;

		virtual void assignDerivedMaterialParameters() override 
		{
			BaseMaterial::assignDerivedMaterialParameters();
		};
	public:
		Fluid() : BaseMaterial(), c_0_(1.0), mu_(0.0),
			fluid_particles_(NULL) {
			material_name_ = "Fluid"; 
		};
		virtual ~Fluid() {};

		void assignFluidParticles(FluidParticles* fluid_particles) 
		{
			fluid_particles_ = fluid_particles;
		};

		Real ReferenceSoundSpeed() { return c_0_; };
		Real ReferenceViscosity() { return mu_; };
		virtual Real getPressure(Real rho) = 0;
		virtual Real getPressure(Real rho, Real rho_e) { return getPressure(rho); };
		virtual Real DensityFromPressure(Real p) = 0;
		virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;
		virtual Fluid* ThisObjectPtr() override { return this; };
	};

	/** @class  Solid
	 *  @brief Base class of all solid materials
	*/
	class Solid : public BaseMaterial
	{
	public:
		Solid() : BaseMaterial(), contact_stiffness_(1.0),
			contact_friction_(0.0), solid_particles_(NULL) 
		{
			material_name_ = "Solid";
		};
		virtual ~Solid() {};

		void assignSolidParticles(SolidParticles* solid_particles) 
		{
			solid_particles_ = solid_particles;
		};

		Real ContactFriction() { return contact_friction_; };
		Real ContactStiffness() { return contact_stiffness_; };
		virtual Solid* ThisObjectPtr() override { return this; };
	protected:
		Real contact_stiffness_; /**< contact-force stiffness related to bulk modulus*/
		Real contact_friction_; /**< friction property mimic fluid viscosity*/
		SolidParticles* solid_particles_;

		virtual void assignDerivedMaterialParameters() override 
		{
			BaseMaterial::assignDerivedMaterialParameters();
		};
	};
}
#endif //BASE_MATERIAL_H