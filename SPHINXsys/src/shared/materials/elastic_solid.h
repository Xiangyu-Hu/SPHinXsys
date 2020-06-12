/**
* @file 	elastic_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These classes are based on isotropic linear elastic solid.
* 			Several more complex materials, including mneo-hookean, FENE noe-hookean
*			and anisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
* @version	0.1
* @version  0.2.1
*           Chi Zhang
*			add the electrophysiology to muscle body.
* @version  0.2.2
*           Chi Zhang
*           Add the electro-mechnaics and local properties of muscle material.
* @version  0.3
*           Xiangyu Hu
*           The relations between the materials are revised.
*/
#pragma once

#include "base_material.h"
#include <fstream>

using namespace std;

namespace SPH {

	/** preclaimed classes. */
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

	public:
		/** Constructor */
		ElasticSolid() : Solid(), elastic_particles_(NULL),
			c_0_(1.0), eta_0_(0.0), lambda_0_(0.375) {};
		virtual ~ElasticSolid() {};

		/** assign particles to this material */
		void assignElasticSolidParticles(ElasticSolidParticles* elastic_particles);
		/** Access to reference sound speed. */
		Real getReferenceSoundSpeed() { return c_0_; };
		/** the interface for dynamical cast*/
		virtual ElasticSolid* PointToThisObject() override { return this; };
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
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) = 0;
		/**
		 * @brief Compute phsical and numerical damping stress.
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] numerical_viscoisty numerical viscosity.
		 * @param[in] particle_index_i 	Particle index.
		 * */
		virtual Matd DampingStress(Matd& deform_grad, Matd& deform_grad_rate, Real numerical_viscoisty, size_t particle_index_i);
	};

	/**
	* @class LinearElasticSolid
	* @brief Isoptropic linear elastic solid
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
		virtual LinearElasticSolid* PointToThisObject() override { return this; };

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** Compute the stree through Constitutive relation. */
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
		virtual NeoHookeanSolid* PointToThisObject() override { return this; };
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
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
		virtual FeneNeoHookeanSolid* PointToThisObject() override { return this; };
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
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
		/** consitutive parameters */
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
			a_0_{1.0, 0.0, 0.0, 0.0 }, b_0_{1.0, 0.0, 0.0, 0.0 }, bulk_modulus_(30.0),
			f0_(0), s0_(0), f0f0_(0), f0s0_(0), s0s0_(0) {
			material_name_ = "Muscle";
		};
		virtual ~Muscle() {};

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** obtain fiber matrix */
		virtual Matd getMuscleFiber(size_t particle_index_i) { return f0f0_; };
		/** the interface for dynamical cast*/
		virtual Muscle* PointToThisObject() override { return this; };
		/** compute the stree through Constitutive relation. */
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
		/** lcoal fiber direction. */
		StdVec<Vecd> local_f0_;
		/** lcoal sheet direction. */
		StdVec<Vecd> local_s0_;

		/** Constructor */
		LocallyOrthotropicMuscle() : Muscle() {
			material_name_ = "LocallyOrthotropicMuscle";
		};
		virtual ~LocallyOrthotropicMuscle() {};

		/** obtain fiber matrix */
		virtual Matd getMuscleFiber(size_t particle_index_i) override { return local_f0f0_[particle_index_i]; };
		/** the interface for dynamical cast*/
		virtual LocallyOrthotropicMuscle* PointToThisObject() override { return this; };
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] Base particles of elastic solid. 
		 */
		virtual void initializeLocalProperties(BaseParticles* base_particles) override;
		/** 
		 * @brief Compute the stree through Constitutive relation.
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
		ActiveMuscle(Muscle* muscle) : ElasticSolid(*muscle), muscle_(*muscle),
			active_muscle_particles_(NULL) {
			material_name_ = "ActiveMuscle";
		};
		virtual ~ActiveMuscle() {};

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** assign particles to this material */
		void assignActiveMuscleParticles(ActiveMuscleParticles* active_muscle_particles);

		/** the interface for dynamical cast*/
		virtual ActiveMuscle* PointToThisObject() override { return this; };
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
