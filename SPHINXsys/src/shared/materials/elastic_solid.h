/**
* @file 	elastic_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These class based on isotropic linear elastic solid.
* 			Several more complex materials, including mneo-hookean, FENE noe-hookean
*			 and aisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
* @version	0.1
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
	protected:
		virtual void SetupLocalProperties(ElasticSolidParticles &elasticsolid_particles) {};
	public:
		ElasticSolid(string elastic_solid_name, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0);
		virtual ~ElasticSolid() {};

		/** the interface for dynamical cast*/
		virtual ElasticSolid* PointToThisObject() override { return this; };

		virtual Real SetSoundSpeed(Real rho_0, Real E, Real nu);
		virtual Real SetShearModulus(Real E, Real nu);
		virtual Real SetLambda(Real E, Real nu);

		/** obtain saved sound speed */
		virtual Real GetSoundSpeed(size_t particle_index_i);
		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i);
		/** compute damping stress */
		virtual Matd DampingStress(Matd &deform_grad, Matd &deform_grad_rate, size_t particle_index_i);
		/** compute artificial viscosity for numerial stability */
		virtual Real GetArtificalViscosity(Real rho, Real sound_speed, Real smoothing_length);
		/** compute damping stress */
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
		NeoHookeanSolid(string elastic_solid_name, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0);
		virtual ~NeoHookeanSolid() {};

		/** the interface for dynamical cast*/
		virtual NeoHookeanSolid* PointToThisObject() override { return this; };

		/** compute elastic stress */
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
		FeneNeoHookeanSolid(string elastic_solid_name, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0, Real j1_m = 1.0);
		virtual ~FeneNeoHookeanSolid() {};

		/** the interface for dynamical cast*/
		virtual FeneNeoHookeanSolid* PointToThisObject() override { return this; };

		/** compute elastic stress */
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
	class Muscle : public ElasticSolid
	{
	protected:
		StdVec<Vecd> f0_, s0_;
		StdVec<Matd> f0f0_, s0s0_, f0s0_;

		virtual void SetupLocalProperties(ElasticSolidParticles &elasticsolid_particles) override;
	public:
		/** consitutive parameters */
		Real a_0_[4], b_0_[4];

		Muscle(string elastic_solid_name, Real a_0[4], Real b_0[4],
			Real rho_0 = 1.0, Real poisson = 0.49, Real eta_0 = 0.0);
		virtual ~Muscle() {};
		
		/** the interface for dynamical cast*/
		virtual Muscle* PointToThisObject() override { return this; };

		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, size_t particle_index_i) override;
	};
}
