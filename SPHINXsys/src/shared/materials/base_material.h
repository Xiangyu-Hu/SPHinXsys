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
* @version	0.1
*/
#pragma once
#include <string>

#include "base_data_package.h"
#include "base_particles.h"

using namespace std;

namespace SPH {
	/** preclaimed classes */
	class FluidParticles;
	class SolidParticles;

	/** @class  BaseMaterial
	 *  @brief Base of all materials
	 *  @details Note that the case dependent material properties will defined in 
	 *  applications.
	*/
	class BaseMaterial
	{
	protected:
		string material_name_;
		/** inverse of dimension */
		Real inv_dimension_;
		/** base particle information for defining local material properties*/
		BaseParticles* base_particles_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() {};
	public:
		/** Default constructor */
		BaseMaterial() : material_name_("BaseMaterial"), 
			inv_dimension_(1.0 / (Real)Vecd(0).size()), base_particles_(NULL) {};
		virtual ~BaseMaterial() {};

		/** Assign base particles to this material for defining local material properties. */
		void assignParticles(BaseParticles* base_particles) { base_particles_ = base_particles; };
		/** The interface for dynamical cast. */
		virtual BaseMaterial* PointToThisObject() { return this; };
		/** access the material name */
		string getMaterialName() { return material_name_;}
		/** initialize the local property. */
		virtual void initializeLocalProperties(BaseParticles* base_particles) {};
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
		/**
		 * @brief Write local material properties particle data in VTU format for Paraview.
		 * @param[inout] output_file Ofstream of particle data.
		 */
		virtual void WriteMaterialPropertyToVtuFile(ofstream& output_file) {};
	};


	/** @class  Fluid
	 *  @brief Base class  of all fluids
	*/
	class Fluid : public BaseMaterial
	{
	protected:
		/** reference density, sound speed, viscosity. */
		Real rho_0_, c_0_, mu_;
		/** particles for this material */
		FluidParticles* fluid_particles_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override {
			BaseMaterial::assignDerivedMaterialParameters();
		};
	public:
		/** constructor with material name. */
		Fluid() : BaseMaterial(), rho_0_(1.0), c_0_(1.0), mu_(0.0),
			fluid_particles_(NULL) {
			material_name_ = "Fluid"; 
		};
		virtual ~Fluid() {};

		/** assign particles to this material */
		void assignFluidParticles(FluidParticles* fluid_particles) {
			fluid_particles_ = fluid_particles;
		};
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
	 *  @brief Base class  of all solids
	*/
	class Solid : public BaseMaterial
	{
	public:
		/** constructor with material name. */
		Solid() : BaseMaterial(), rho_0_(1.0), collision_stiffness_(1.0),
			collision_friction_(0.0), solid_particles_(NULL) {
			material_name_ = "Solid";
		};
		virtual ~Solid() {};

		/** assign particles to this material */
		void assignSolidParticles(SolidParticles* solid_particles) {
			solid_particles_ = solid_particles;
		};
		/** Access to reference density. */
		Real getReferenceDensity() { return rho_0_; };
		Real getFriction() { return collision_friction_; };
		Real getStiffness() { return collision_stiffness_; };
		/** the interface for dynamical cast*/
		virtual Solid* PointToThisObject() override { return this; };

	protected:
		/** reference density */
		Real rho_0_;
		/** artifical bulk modulus*/
		Real collision_stiffness_;
		/** friction property mimic fluid viscosity*/
		Real collision_friction_;

		/** particles for this material */
		SolidParticles* solid_particles_;

		/** assign derived material properties*/
		virtual void assignDerivedMaterialParameters() override {
			BaseMaterial::assignDerivedMaterialParameters();
		};
	};
	
}
