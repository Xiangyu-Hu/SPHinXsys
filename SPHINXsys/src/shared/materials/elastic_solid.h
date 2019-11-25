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
#include "solid_particles.h"

using namespace std;

namespace SPH {

	/**
	* @class ElasticSolid
	* @brief Isoptropic linear elastic solid
	*/
	class ElasticSolid : public Solid
	{
	public:
		ElasticSolid(string elastic_solid_name, SPHBody *body,
			Real rho_0 = 1.0, Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0);
		virtual ~ElasticSolid() {};

		/** the interface for dynamical cast*/
		virtual ElasticSolid* PointToThisObject() override { return this; };
		/**
		 * @brief Set the speed of sound.
		 * @param[in] rho_0 initial density.
		 * @param[in] Young modulus
		 * @param[in] nu Possion ratio
		 */
		virtual Real SetSoundSpeed(Real rho_0, Real E, Real nu);
		/**
		 * @brief Set the shear modulus
		 * @param[in] Young modulus
		 * @param[in] nu Possion ratio.
		 */
		virtual Real SetShearModulus(Real E, Real nu);
		/**
		 * @brief Set lambda
		 * @param[in] Young modulus
		 * @param[in] nu Possion ratio. 
		 */
		virtual Real SetLambda(Real E, Real nu);
		/**
		 * @brief Get the speed of sound.
		 * @param[in] particle_index_i Particle index
		 */
		virtual Real GetSoundSpeed(size_t particle_index_i);
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] elasticsolid_particles Particles of elastic solid. 
		 */
		virtual void SetupLocalProperties(SPHBody *body) {};
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i);
		/** 
		 * @brief Compute damping stress.
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] particle_index_i 	Particle index. 
		 * */
		virtual Matd DampingStress(Matd &deform_grad, Matd &deform_grad_rate, size_t particle_index_i);
		/** 
		 * @brief Compute artificial viscosity for numerial stability 
		 * @param[in] rho 	Density. 
		 * @param[in] sound_speed Speed of sound.
		 * @param[in] smoothing_length 	Smoothing lenght. 
		 */
		virtual Real GetArtificalViscosity(Real rho, Real sound_speed, Real smoothing_length);
		/** 
		 * @brief Numericla damping stress.
		 * @param[in] deform_grad 	Gradient of deformation.
		 * @param[in] deform_grad_rate 	Rate of gradient of deformation.
		 * @param[in] numerical_viscoisty 	Nuemrical visocsity.  
		 * */
		virtual Matd NumericalDampingStress(Matd &deform_grad, Matd &deform_grad_rate, Real numerical_viscoisty);
		/** Youngs modulus, poisson ration and physical viscosity */
		Real E_0_, nu_, eta_0_, c_0_;
		/** Lame parameters */
		Real lambda_0_, G_0_;
	};

	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid
	*/
	class NeoHookeanSolid : public ElasticSolid
	{
	public:
		NeoHookeanSolid(string elastic_solid_name, SPHBody *body, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0);
		virtual ~NeoHookeanSolid() {};
		/** the interface for dynamical cast*/
		virtual NeoHookeanSolid* PointToThisObject() override { return this; };
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i) override;
	};
	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid with finite extension
	*/
	class FeneNeoHookeanSolid : public ElasticSolid
	{
	protected:
		/** reference extension */
		Real j1_m_;

	public:
		FeneNeoHookeanSolid(string elastic_solid_name, SPHBody *body, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0, Real j1_m = 1.0);
		virtual ~FeneNeoHookeanSolid() {};
		/** the interface for dynamical cast*/
		virtual FeneNeoHookeanSolid* PointToThisObject() override { return this; };
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i) override;
	};

	/**
	* @class Muscle
	* @brief muscle model is a anisotropic material in which
	* there are local fiber direction and cross-fiber sheet direction.
	* the model here is from 
	* Holzapfel and Ogden, 2009, Phil. Trans. R. Soc. 367:3445-3475
	* we consider a neo-hookean model for the background isotropic contribution.
	*/
	//using ReferenceNeighborDiffusionList = StdVec<ReferenceNeighboringParticleDiffusion>;

	class Muscle : public ElasticSolid
	{
	protected:
		StdVec<Vecd> f0_, s0_;				/**< Fiber direction. */
		StdVec<Matd> f0f0_, s0s0_, f0s0_, diff_cd_0;	/**< Sheet direction. */
	public:
		/** consitutive parameters */
		Real a_0_[4], b_0_[4];
		Vecd d_0_;

		/** Reference inner diffusion tensor configurations for totoal Lagrangian formulation. */
		//StdVec<ReferenceNeighborDiffusionList> reference_inner_diffusionn_tensor_;

		Muscle(string elastic_solid_name, SPHBody *body, Real a_0[4], Real b_0[4], Vecd d_0,
			Real rho_0 = 1.0, Real poisson = 0.49, Real eta_0 = 0.0);
		virtual ~Muscle() {};

		/** the interface for dynamical cast*/
		virtual Muscle* PointToThisObject() override { return this; };
		/**
		 * @brief Setup the local properties, fiber and sheet direction.
		 * @param[in] elasticsolid_particles Particles of elastic solid. 
		 */
		virtual void SetupLocalProperties(SPHBody *body) override {};
		/**
		 * @brief compute the stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i) override;
		/**
		 * @brief compute the active stree through Constitutive relation.
		 * @param[in] deform_grad defromation gradient
		 * @param[in] T_a active contraction stress
		 * @param[in] particle_index_i Particle index
		 */
		virtual Matd ConstitutiveRelationOfActiveStress(Matd &deform_grad, Real T_a, size_t particle_index_i);
		/**
		 * @brief get the diffusion tensor of muscle material.
		 * @param[in]	pnt Position of particles.
		 */
		Matd getDiffussionTensor(size_t particle_index_i);
		/**
		 * @brief get the trace of diffusion matrix.
		 * @param[in]	pnt Position of particles.
		 */
		Real getDiffusionTensorTrace(size_t particle_index_i);
		// /**
		//  * @brief get the trace of diffusion matrix.
		//  * @param[in]	index_i Particle index i
		//  * @param[in] 	index_j Paeticle index j
		//  */
		// Matd getReferenceAverageDiffusionTensor(size_t index_i, size_t index_j);
	};
}
