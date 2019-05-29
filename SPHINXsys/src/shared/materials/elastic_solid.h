/**
* @file 	base_solid.h
* @brief 	These are classes for define properties of elastic solid materials.
*			These class based on isotropic linear elastic solid.
* 			Several more complex materials, including mneo-hookean, FENE noe-hookean
*			 and aisotropic muscle, are derived from the basic elastic solid class.
* @author	Xiangyu Hu and Chi Zhang
* @version	0.1
*/
#pragma once

#include "base_solid.h"

using namespace std;

namespace SPH {

	/**
	* @class ElasticSolid
	* @brief Isoptropic linear elastic solid
	*/
	class ElasticSolid : public Solid
	{
	public:

		ElasticSolid(string elastic_solid_name, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0);
		virtual ~ElasticSolid() {};

		virtual Real GetSoundSpeed(Real rho_0, Real E, Real nu);
		virtual Real GetShearModulus(Real E, Real nu);
		virtual Real GetLambda(Real E, Real nu);

		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, Real local_G, Real local_lambda);
		/** compute damping stress */
		virtual Matd DampingStress(Matd &deform_grad, Matd &deform_grad_rate, Real local_eta);
		/** compute artificial viscosity for numerial stability */
		virtual Real GetArtificalViscosity(Real rho, Real sound_speed, Real smoothing_length);

		//Youngs modulus, poisson ration and physical viscosity
		Real E_0_, nu_, eta_0_;
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

		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, Real local_G, Real local_lambda) override;
	};

	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid with finite extension
	*/
	class FeneNeoHookeanSolid : public ElasticSolid
	{
	protected:
		Real j1_m_;

	public:
		FeneNeoHookeanSolid(string elastic_solid_name, Real rho_0 = 1.0,
			Real E_0 = 1.0, Real poisson = 0.3, Real eta_0 = 0.0, Real j1_m = 1.0);
		virtual ~FeneNeoHookeanSolid() {};

		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, Real local_G, Real local_lambda) override;
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
	public:
		Muscle(string elastic_solid_name, Real a_0[4], Real b_0[4],
			Real rho_0 = 1.0, Real poisson = 0.49, Real eta_0 = 0.0);
		virtual ~Muscle() {};
		
		/** consitutive parameters */
		Real a_0_[4], b_0_[4];

		/** compute elastic stress */
		virtual Matd ConstitutiveRelation(Matd &deform_grad, 
			Real local_lambda, Real a[4], Real b[4], Vecd local_f0, Vecd local_s0,
			Matd local_f0f0, Matd local_s0s0, Matd local_f0s0);
	};
}
