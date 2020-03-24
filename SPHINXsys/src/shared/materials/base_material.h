/**
 * @file 	base_material.h
 * @brief 	This is the base classes of all materials. 
 *		    A function in a derived material class returns a value with the inputs
 *          from the particle data.
 *			Basically, it is a interface from which
 *			one can access devirved material by dynamic cast.
 *          Note that the derived material may have position dependent or 
 *          local properties.
* @author	Chi Zhang and Xiangyu Hu
* @version	0.1
* @version  0.2.1
*           Chi Zhang
*			add the electrophysiology to muscle body.
* @version  0.2.2
*           Chi Zhang
*           Add the electro-mechnaics and local properties of muscle material.
*/
#pragma once
#include <string>

#include "base_data_package.h"
#include "base_particles.h"

using namespace std;

namespace SPH {
	/** @class  Material
	 *  @brief Base of all materials
	 *  @details Note that the same material object can be shared by several
	 *  SPH bodies, and the case dependent material properties will defined in 
	 *  apliications.
	*/
	class BaseMaterial
	{
	protected:
		string material_name_;
		/** base partilce information for defining local material properties*/
		BaseParticles* base_particles_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() {};
	public:
		/** Default constructor */
		BaseMaterial() : material_name_("BaseMaterial"), base_particles_(NULL) {};
		/** constructor with material name. */
		BaseMaterial(string material_name) : BaseMaterial() {
			material_name_ = material_name;
		};
		virtual ~BaseMaterial() {};

		/** Assign base particles to this material
		  * for defining local material properties. 
		  */
		void AssignParticles(BaseParticles* base_particles) { base_particles_ = base_particles; };
		/** The interface for dynamical cast. */
		virtual BaseMaterial* PointToThisObject() { return this; };
		/**
		 * @brief Return the material name.
		 */
		string getMaterialName() { return material_name_;}
		/** 
		 * @brief Write the material property to xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void writeToXmlForReloadMaterialProperty(std::string &filefullpath) {};
		/** 
		 * @brief Read the material property from xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void readFromXmlForMaterialProperty(std::string &filefullpath) {};
	};


	/** @class  Fluid
	 *  @brief Base calss  of all fluids
	*/
	class Fluid : public BaseMaterial
	{
	protected:
		/** reference density, sound speed, viscosity. */
		Real rho_0_, c_0_, mu_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() {
			BaseMaterial::assignDerivedMaterialParameters();
		};
	public:
		/** constructor with material name. */
		Fluid(string fluid_name) : BaseMaterial(fluid_name), 
			rho_0_(1.0), c_0_(1.0), mu_(0.0) {};
		virtual ~Fluid() {};

		/** the interface for dynamical cast*/
		virtual Fluid* PointToThisObject() override { return this; };

		Real GetReferenceSoundSpeed() { return c_0_; };
		Real GetReferenceDensity() { return rho_0_; };
		Real getReferenceViscosity() { return mu_; };
		virtual Real GetPressure(Real rho) = 0;
		virtual Real GetPressure(Real rho, Real rho_e) { return GetPressure(rho); };
		virtual Real ReinitializeRho(Real p) = 0;
		virtual Real GetSoundSpeed(Real p = 0.0, Real rho = 1.0) = 0;
		virtual Real RiemannSolverForPressure(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;
		virtual Real RiemannSolverForVelocity(Real rhol, Real Rhor, Real pl, Real pr, Real ul, Real ur) = 0;
	};

	/** @class  Solid
	 *  @brief Base calss  of all solids
	*/
	class Solid : public BaseMaterial
	{
	protected:
		/** reference density */
		Real rho_0_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialProperties() {
			BaseMaterial::assignDerivedMaterialParameters();
		};
	public:
		/** constructor with material name. */
		Solid(string solid_name) 
			: BaseMaterial(solid_name), rho_0_(1.0) {};
		virtual ~Solid() {};

		/** Access to reference density. */
		Real getReferenceDensity() { return rho_0_; };
		/** the interface for dynamical cast*/
		virtual Solid* PointToThisObject() override { return this; };
		/** 
		 * @brief Write the material property to xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void writeToXmlForReloadMaterialProperty(std::string &filefullpath) override {};
		/** 
		 * @brief Read the material property from xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void readFromXmlForMaterialProperty(std::string &filefullpath) override {};
	};
	
	/**
	  * @class MaterialWithLocalProperties
	  * @brief The base class for a material with local properties
	  */
	template<class MaterialType>
	class MaterialWithLocalProperties : public MaterialType
	{
	protected:
		/** Assigning local material properties. */
		virtual void AssignLocalProperties() = 0;
	public:
		/** constructor with material name. */
		MaterialWithLocalProperties<MaterialType>(string material_name)
			: MaterialType(material_name) {};
		virtual ~MaterialWithLocalProperties() {};

		/** the interface for dynamical cast*/
		virtual MaterialWithLocalProperties* PointToThisObject() override { return this; };
	};

	/**
	  * @class CompositeMaterial
	  * @brief The base class for a material compaite with several
	  * similar materials
	  */
	template<class MaterialType>
	class CompositeMaterial : public MaterialType
	{
	protected:
		/** local material mapping. */
		StdLargeVec<MaterialType*> matrial_mapping_;
		/** Vector of materal objects. */
		StdVec<MaterialType*> reference_materials_;

		/** Assigning local material properties. */
		virtual void AssignComposite() = 0;
	public:
		CompositeMaterial<MaterialType>(string material_name,
			StdVec<MaterialType*> reference_materials) 
			: MaterialType(material_name),
			reference_materials_(reference_materials) {};
		virtual ~CompositeMaterial() {};

		/** the interface for dynamical cast*/
		virtual CompositeMaterial* PointToThisObject() override { return this; };
	}; 
}
