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
* @file 	elastic_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These classes are based on isotropic linear elastic solid.
* 			Several more complex materials, including neo-hookean, FENE noe-hookean
*			and anisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
*/
#pragma once

#include "base_material.h"
#include <fstream>

using namespace std;

namespace SPH {

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class ElasticSolidParticles;
	class ActiveMuscleParticles;

	/**
	* @class ElasticSolid
	* @brief Abstract class for a generalized elastic solid
	*/
	class ElasticSolid : public Solid
	{
	protected:
		Real eta_0_; 		/*> physical viscosity */
		Real c_0_; 			/*> speed of sound */
		Real lambda_0_; 	/*> First Lame parameter */
		ElasticSolidParticles* elastic_particles_;

		/** set contact stiffness for contact dynamics */
		void setContactStiffness() { contact_stiffness_ = c_0_* c_0_; };
	public:
		/** Constructor */
		ElasticSolid() : 
		Solid(), eta_0_(0.0), c_0_(1.0), lambda_0_(0.375), elastic_particles_(NULL) {};
		virtual ~ElasticSolid() {};

		/** assign particles to this material */
		void assignElasticSolidParticles(ElasticSolidParticles* elastic_particles);
		Real ReferenceSoundSpeed() { return c_0_; };
		Real getPhysicalViscosity() { return eta_0_; };
		/** the interface for dynamical cast*/
		virtual ElasticSolid* pointToThisObject() override { return this; };
		/**
		 * @brief Get the speed of sound.
		 * @param[in] particle_index_i Particle index
		 */
		 /** the speed of sound. */
		virtual Real SetSoundSpeed() = 0;
		/** Get viscous time step size. */
		virtual Real getViscousTimeStepSize(Real smoothing_length);
		/** Get numerical viscosity. */
		virtual Real getNumericalViscosity(Real smoothing_length);
		/**
		 * @brief compute the stress through Constitutive relation.
		 * In the present total Lagrangian formulation, second Piola-Kirchhoff stress is used.  
		 * @param[in] deform_grad deformation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) = 0;
		/**
		 * @brief Compute physical and numerical damping stress.
		 * In the present total Lagrangian formulation, the viscous stress is also in form of second Piola-Kirchhoff stress.  
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] numerical_viscosity numerical viscosity.
		 * @param[in] particle_index_i 	Particle index.
		 * */
		virtual Matd NumericalDampingStress(Matd& deform_grad, 
		Matd& deform_grad_rate, Real numerical_viscosity, size_t particle_index_i);
	};

	/**
	* @class LinearElasticSolid
	* @brief Isotropic linear elastic solid
	*/
	class LinearElasticSolid : public ElasticSolid
	{
	protected:
		Real E_0_;	/*> Youngs modulus */
		Real nu_;	/*> poisson ratio */
		/**Second Lame parameter, shear modulus. Derived material property. */
		Real G_0_; 

		/** Set the shear modulus. */
		virtual Real SetShearModulus();
		/** Set lambda. */
		virtual Real SetLambda();
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override;
	public:
		/** Constructor */
		LinearElasticSolid() : ElasticSolid(),
			E_0_(1.0), nu_(1.0 / 3.0), G_0_(0.75) {
			material_name_ = "LinearElasticSolid";
		};
		virtual ~LinearElasticSolid() {};
		/** the interface for dynamical cast*/
		virtual LinearElasticSolid* pointToThisObject() override { return this; };

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** Compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
	};

	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid
	*/
	class NeoHookeanSolid : public LinearElasticSolid
	{
	public:
		/** Constructor. */
		NeoHookeanSolid() : LinearElasticSolid() {
			material_name_ = "NeoHookeanSolid";
		};
		virtual ~NeoHookeanSolid() {};
		/** the interface for dynamical cast*/
		virtual NeoHookeanSolid* pointToThisObject() override { return this; };
		/**
		 * @brief compute the stress through Constitutive relation.
		 * @param[in] deform_grad deformation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
	};
	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid with finite extension
	*/
	class FeneNeoHookeanSolid : public LinearElasticSolid
	{
	protected:
		/** reference extension */
		Real j1_m_;

	public:
		/** Constructor */
		FeneNeoHookeanSolid() : LinearElasticSolid(), j1_m_(1.0) {
			material_name_ = "FeneNeoHookeanSolid";
		};
		virtual ~FeneNeoHookeanSolid() {};
		/** the interface for dynamical cast*/
		virtual FeneNeoHookeanSolid* pointToThisObject() override { return this; };
		/**
		 * @brief compute the stress through Constitutive relation.
		 * @param[in] deform_grad deformation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
	};

	/**
	* @class Muscle
	* @brief Globally orthotropic muscle. 
	*/
	class Muscle : public ElasticSolid
	{
	protected:
		Vecd f0_, s0_; 				/**< Reference fiber and sheet directions. */
		Matd f0f0_, s0s0_, f0s0_;	/**< Direct products of fiber and sheet directions. */
		/** constitutive parameters */
		Real a_0_[4], b_0_[4];
		/** reference stress to achieve weakly compressible condition */
		Real bulk_modulus_;

		/** Set lambda. */
		virtual Real SetLambda();
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override;
	public:
		/** Constructor */
		Muscle() : ElasticSolid(),
			f0_(0), s0_(0), f0f0_(0), s0s0_(0), f0s0_(0),
			a_0_{1.0, 0.0, 0.0, 0.0 }, b_0_{1.0, 0.0, 0.0, 0.0 }, bulk_modulus_(30.0) {
			material_name_ = "Muscle";
		};
		virtual ~Muscle() {};

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** obtain fiber matrix */
		virtual Matd getMuscleFiber(size_t particle_index_i) { return f0f0_; };
		/** the interface for dynamical cast*/
		virtual Muscle* pointToThisObject() override { return this; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
	};

	/**
	* @class LocallyOrthotropicMuscle
	* @brief muscle model is a anisotropic material in which
	* there are local fiber direction and cross-fiber sheet direction.
	* the model here is from
	* Holzapfel and Ogden, 2009, Phil. Trans. R. Soc. 367:3445-3475
	* we consider a neo-hookean model for the background isotropic contribution.
	*/
	class LocallyOrthotropicMuscle : public Muscle
	{
	protected:
		StdVec<Matd> local_f0f0_, local_s0s0_, local_f0s0_;			/**< Sheet direction. */
		virtual void assignDerivedMaterialParameters() override {
			Muscle::assignDerivedMaterialParameters();
		};
	public:
		/** local fiber direction. */
		StdVec<Vecd> local_f0_;
		/** local sheet direction. */
		StdVec<Vecd> local_s0_;

		/** Constructor */
		LocallyOrthotropicMuscle() : Muscle() {
			material_name_ = "LocallyOrthotropicMuscle";
		};
		virtual ~LocallyOrthotropicMuscle() {};

		/** obtain fiber matrix */
		virtual Matd getMuscleFiber(size_t particle_index_i) override { return local_f0f0_[particle_index_i]; };
		/** the interface for dynamical cast*/
		virtual LocallyOrthotropicMuscle* pointToThisObject() override { return this; };
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] Base particles of elastic solid. 
		 */
		virtual void initializeLocalProperties(BaseParticles* base_particles) override;
		/** 
		 * @brief Compute the stress through Constitutive relation.
		 * @param[in] deform_grad Deformation tensor.
		 * @param[in] particle_index_i Particle index.
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
		/** Write the material property to xml file. */
		virtual void writeToXmlForReloadMaterialProperty(std::string &filefullpath) override;
		/** Read the material property from xml file. */
		virtual void readFromXmlForMaterialProperty(std::string &filefullpath) override;
		/** Write local material properties particle data in VTU format for Paraview */
		virtual void WriteMaterialPropertyToVtuFile(ofstream& output_file) override;
	};

	/**
	* @class ActiveMuscle
	* @brief Here, the active reponse is considered.
	*/
	class ActiveMuscle : public ElasticSolid
	{
	protected:
		ActiveMuscleParticles* active_muscle_particles_;

		Muscle& muscle_;
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override {
			ElasticSolid::assignDerivedMaterialParameters();
		};
	public:
		/** Constructor. */
		ActiveMuscle(Muscle* muscle) : 
		ElasticSolid(*muscle),active_muscle_particles_(NULL), muscle_(*muscle) {
			material_name_ = "ActiveMuscle";
		};
		virtual ~ActiveMuscle() {};

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** assign particles to this material */
		void assignActiveMuscleParticles(ActiveMuscleParticles* active_muscle_particles);

		/** the interface for dynamical cast*/
		virtual ActiveMuscle* pointToThisObject() override { return this; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
		/** Write the material property to xml file */
		virtual void writeToXmlForReloadMaterialProperty(std::string& filefullpath) override;
		/** Read the material property from xml file. */
		virtual void readFromXmlForMaterialProperty(std::string& filefullpath) override;
		/** Write local material properties particle data in VTU format for Paraview */
		virtual void WriteMaterialPropertyToVtuFile(ofstream& output_file) override;
	};
}
