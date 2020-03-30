/**
* @file 	elastic_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These class based on isotropic linear elastic solid.
* 			Several more complex materials, including mneo-hookean, FENE noe-hookean
*			 and aisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
* @version	0.1
* @version  0.2.1
*           Chi Zhang
*			add the electrophysiology to muscle body.
* @version  0.2.2
*           Chi Zhang
*           Add the electro-mechnaics and local properties of muscle material.
*/
#pragma once

#include "base_material.h"
#include "diffusion_reaction.h"
#include "xml_engine.h"
#include <fstream>

using namespace std;

namespace SPH {
	/**
	* @class ElasticSolid
	* @brief general elastic solid
	*/
	class ElasticSolid : public Solid
	{
	protected:
		Real eta_0_; 		/*> physical viscosity */
		Real c_0_; 			/*> speed of sound */
		Real lambda_0_; 	/*> First Lame parameter */

		/** the speed of sound. */
		virtual Real SetSoundSpeed() = 0;
	public:
		/** Constructor */
		ElasticSolid(string elastic_solid_name) : Solid(elastic_solid_name),
			c_0_(1.0), eta_0_(0.0), lambda_0_(0.375) {};
		virtual ~ElasticSolid() {};
		/** Access to reference sound speed. */
		Real getReferenceSoundSpeed() { return c_0_; };
		/** the interface for dynamical cast*/
		virtual ElasticSolid* PointToThisObject() override { return this; };
		/**
		 * @brief Get the speed of sound.
		 * @param[in] particle_index_i Particle index
		 */
		virtual Real GetSoundSpeed(size_t particle_index_i);
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
		 * @brief Compute damping stress.
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] particle_index_i 	Particle index.
		 * */
		virtual Matd DampingStress(Matd& deform_grad, Matd& deform_grad_rate, size_t particle_index_i);
		/**
		 * @brief Numericla damping stress.
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] numerical_viscoisty 	Nuemrical visocsity. */
		virtual Matd NumericalDampingStress(Matd& deform_grad, Matd& deform_grad_rate, Real numerical_viscoisty);
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
	* @class LinearElasticSolid
	* @brief Isoptropic linear elastic solid
	*/
	class LinearElasticSolid : public ElasticSolid
	{
	protected:
		/** Youngs modulus, poisson ration. */
		Real E_0_, nu_;
	
		/**Second Lame parameter, shear modulus. Derived material property. */
		Real G_0_; 

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** Set the shear modulus. */
		virtual Real SetShearModulus();
		/** Set lambda. */
		virtual Real SetLambda();
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override;
	public:
		/** Constructor */
		LinearElasticSolid(string linear_elastic_solid_name)
			: ElasticSolid(linear_elastic_solid_name),
			E_0_(1.0), nu_(1.0 / 3.0), G_0_(0.75) {};
		virtual ~LinearElasticSolid() {};
		/** the interface for dynamical cast*/
		virtual LinearElasticSolid* PointToThisObject() override { return this; };
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
		NeoHookeanSolid(string neohookean_solid_name)
			: LinearElasticSolid(neohookean_solid_name) {};
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
		FeneNeoHookeanSolid(string feneneohookean_solid_name)
			: LinearElasticSolid(feneneohookean_solid_name), j1_m_(1.0) {};
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
		/** Poisson ratio. */
		Real nu_;

		/** the speed of sound. */
		virtual Real SetSoundSpeed() override;
		/** Set lambda. */
		virtual Real SetLambda();
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override;
	public:
		/** Constructor */
		Muscle(string muscle_name) : ElasticSolid(muscle_name),
			a_0_{1.0, 0.0, 0.0, 0.0 }, b_0_{1.0, 0.0, 0.0, 0.0 }, nu_(0.49),
			f0_(0), s0_(0), f0f0_(0), f0s0_(0), s0s0_(0) {};
		virtual ~Muscle() {};

		/** the interface for dynamical cast*/
		virtual Muscle* PointToThisObject() override { return this; };
		/** compute the stree through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
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
		StdVec<Vecd> local_s0_;										/**< Fiber direction. */
		StdVec<Matd> local_f0f0_, local_s0s0_, local_f0s0_;			/**< Sheet direction. */
		virtual void assignDerivedMaterialParameters() override
		{
			Muscle::assignDerivedMaterialParameters();
		};
	public:
		StdVec<Vecd> local_f0_;
		/** Constructor */
		LocallyOrthotropicMuscle(string muscle_name) : Muscle(muscle_name) {};
		virtual ~LocallyOrthotropicMuscle() {};

		/** the interface for dynamical cast*/
		virtual LocallyOrthotropicMuscle* PointToThisObject() override { return this; };
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] elasticsolid_particles Particles of elastic solid. 
		 */
		virtual void SetupLocalProperties(SPHBody *body) { };
		/** 
		 * @brief Compute the stree through Constitutive relation.
		 * @param[in] deform_grad Deformation tensor.
		 * @param[in] particle_index_i Particle index.
		 */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
		/** 
		 * @brief Write the material property to xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void writeToXmlForReloadMaterialProperty(std::string &filefullpath) override;
		/** 
		 * @brief Read the material property from xml file.
		 * @param[in] filefullpath Full path the the output file.
		 */
		virtual void readFromXmlForMaterialProperty(std::string &filefullpath) override;
	};

	/**
	* @class ActiveMuscle
	* @brief Here, the active reponse is considered.
	*/
	class ActiveMuscle : public Muscle
	{
	protected:
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override 
		{
		
			Muscle::assignDerivedMaterialParameters();
		};
	public:
		/** Constructor. */
		ActiveMuscle(string active_muscle_name)
			: Muscle(active_muscle_name) {};	
		virtual ~ActiveMuscle() {};

		/** the interface for dynamical cast*/
		virtual ActiveMuscle* PointToThisObject() override { return this; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
	};

	/**
	* @class ActiveMuscle
	* @brief Here, the active reponse is considered.
	*/
	class ActiveLocallyOrthotropicMuscle : public LocallyOrthotropicMuscle
	{
	protected:
		/** assign derived material properties. */
		virtual void assignDerivedMaterialParameters() override 
		{
		
			LocallyOrthotropicMuscle::assignDerivedMaterialParameters();
		};
	public:
		/** Constructor. */
		ActiveLocallyOrthotropicMuscle(string active_muscle_name)
			: LocallyOrthotropicMuscle(active_muscle_name) {};	
		virtual ~ActiveLocallyOrthotropicMuscle() {};

		/** the interface for dynamical cast*/
		virtual ActiveLocallyOrthotropicMuscle* PointToThisObject() override { return this; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deform_grad, size_t particle_index_i) override;
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] elasticsolid_particles Particles of elastic solid. 
		 */
		virtual void SetupLocalProperties(SPHBody *body) override {};
	};
}
