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

#ifndef ELASTIC_SOLID_H
#define ELASTIC_SOLID_H


#include "base_material.h"
#include <fstream>

namespace SPH {

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class ElasticSolidParticles;

	/**
	* @class ElasticSolid
	* @brief Abstract class for a generalized elastic solid
	*/
	class ElasticSolid : public Solid
	{
	protected:
		Real c0_; 			/*< sound wave speed */
		Real ct0_; 			/*< tensile wave speed */
		Real cs0_;			/*< shear wave speed */
		Real E0_; 			/*< Youngs or tensile modules  */
		Real G0_; 			/*< shearmodules  */
		Real K0_; 			/*< bulkmodules  */
		Real nu_; 			/*< Poisson ratio  */
		ElasticSolidParticles* elastic_particles_;

		virtual void setReferenceSoundSpeed() = 0;
		virtual void setTensileWaveSpeed() = 0;
		virtual void setShearWaveSpeed() = 0;
		virtual void setYoungsModulus() = 0;
		virtual void setShearModulus() = 0;
		virtual void setBulkModulus() = 0;
		virtual void setPoissonRatio() = 0;
		void setContactStiffness() { contact_stiffness_ = c0_* c0_; };

		virtual void assignDerivedMaterialParameters() override;
	public:
		ElasticSolid() : Solid(), c0_(1.0), ct0_(1.0), cs0_(0.0717), 
			E0_(1.0), G0_(0.5), K0_(1.0), nu_(0.0), elastic_particles_(nullptr) {};
		virtual ~ElasticSolid() {};

		virtual void assignElasticSolidParticles(ElasticSolidParticles* elastic_particles);
		Real ReferenceSoundSpeed() { return c0_; };
		Real TensileWaveSpeed() { return ct0_; };
		Real ShearWaveSpeed() { return cs0_; };
		Real YoungsModulus() { return E0_; };
		Real ShearModulus() { return G0_; };
		Real BulkModulus() { return K0_; };
		Real PoissonRatio() { return nu_; };

		/** compute the stress through defoemation, which can be green-lagrangian tensor, left or right cauchy tensor. */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) = 0;
		//TODO: the NumericalViscosity and NumericalDampingStress to be delete after shell model done.
		virtual Real NumericalViscosity(Real smoothing_length);
		/** Compute numerical damping stress using right cauchy tensor. */
		virtual Matd NumericalDampingRightCauchy(Matd& deformation, Matd& deformation_rate, Real smoothing_length, size_t particle_index_i);
		/** Compute numerical damping stress using left cauchy tensor. */
		virtual Matd NumericalDampingLeftCauchy(Matd& deformation, Matd& deformation_rate, Real smoothing_length, size_t particle_index_i);
		/** numerical demaping is computed between particles i and j */
		virtual Real NumericalDamping(Real dE_dt_ij, Real smoothing_length);

		/** Deviatoric Kirchhoff stress related with the deviatoric part of left cauchy-green deformation tensor.
		 *  Note that, dependent of the normalizeation of the later, the returned stress can be normalized or non-normalized. */
		virtual Matd DeviatoricKirchhoff(const Matd& deviatoric_be);
		/** Volumetric Kirchhoff stress determinate */
		virtual Real VolumetricKirchhoff(Real J) = 0;

		virtual ElasticSolid* ThisObjectPtr() override {return this;};
	};

	/**
	* @class LinearElasticSolid
	* @brief Isotropic linear elastic solid.
	* Note that only basic parameters are used to set ElasticSolid parmaters 
	*/
	class LinearElasticSolid : public ElasticSolid
	{
	public:
		LinearElasticSolid() : ElasticSolid(), youngs_modulus_(1.0), poisson_ratio_(0), lambda0_(1.0)
		{
			material_name_ = "LinearElasticSolid";
		};
		LinearElasticSolid(Real rho_0, Real Youngs_modulus, Real poisson) : ElasticSolid()
		{
			material_name_ = "LinearElasticSolid";
			rho0_ = rho_0;
			youngs_modulus_ = Youngs_modulus;
			poisson_ratio_ = poisson;

			assignDerivedMaterialParameters();
		};
		virtual ~LinearElasticSolid() {};

		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
	protected:
		Real youngs_modulus_; 		/*< Youngs modules as basic inpiut parameter */
		Real poisson_ratio_; 		/*< Poisson ratio as basic inpiut parameter */
		Real lambda0_; 				/*< first Lame parameter */

		virtual void setReferenceSoundSpeed() override;
		virtual void setTensileWaveSpeed() override;
		virtual void setShearWaveSpeed() override;
		virtual void setYoungsModulus() override { E0_ = youngs_modulus_; };
		virtual void setShearModulus() override;
		virtual void setBulkModulus() override;
		virtual void setPoissonRatio() override { nu_ = poisson_ratio_; };
		virtual void assignDerivedMaterialParameters() override;
	private:
		Real getBulkModulus();
		Real getShearModulus();
		Real getLambda();
	};

	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid
	*/
	class NeoHookeanSolid : public LinearElasticSolid
	{
	public:
		NeoHookeanSolid() : LinearElasticSolid() 
		{
			material_name_ = "NeoHookeanSolid";
		};
		NeoHookeanSolid(Real rho_0, Real Youngs_modulus, Real poisson)
			: LinearElasticSolid(rho_0, Youngs_modulus, poisson)
		{
			material_name_ = "NeoHookeanSolid";
		};
		virtual ~NeoHookeanSolid() {};
	
		/** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
	};

	/**
	* @class FeneNeoHookeanSolid
	* @brief Neo-Hookean solid with finite extension
	*/
	class FeneNeoHookeanSolid : public LinearElasticSolid
	{
	protected:
		Real j1_m_; /**< reference extension as basic paramter */
	public:
		FeneNeoHookeanSolid() : LinearElasticSolid(), j1_m_(1.0) {
			material_name_ = "FeneNeoHookeanSolid";
		};
		virtual ~FeneNeoHookeanSolid() {};
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
	};

	/**
	* @class Muscle
	* @brief Globally orthotropic muscle. 
	*/
	class Muscle : public ElasticSolid
	{
	public:
		Muscle() : ElasticSolid(),
			f0_(0), s0_(0), f0f0_(0), s0s0_(0), f0s0_(0),
			a0_{ 1.0, 0.0, 0.0, 0.0 }, b0_{ 1.0, 0.0, 0.0, 0.0 }, bulk_modulus_(30.0)
		{
			material_name_ = "Muscle";
		};
		virtual ~Muscle() {};

		virtual Matd MuscleFiberDirection(size_t particle_index_i) { return f0f0_; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress determinate */
		virtual Real VolumetricKirchhoff(Real J) override;

		virtual Muscle* ThisObjectPtr() override { return this; };
	protected:
		Vecd f0_, s0_; 				/**< Reference fiber and sheet directions as basic parameter. */
		Matd f0f0_, s0s0_, f0s0_;	/**< Tensor products of fiber and sheet directions as basic parameter.. */
		Real a0_[4], b0_[4]; 		/**< constitutive parameters  as basic parameter.*/
		Real bulk_modulus_;			/**< to achieve weakly compressible condition  as basic parameter.*/
		Real lambda0_; 				/*< first Lame parameter */

		virtual void setReferenceSoundSpeed() override;
		virtual void setTensileWaveSpeed() override;
		virtual void setShearWaveSpeed() override;
		virtual void setYoungsModulus() override;
		virtual void setShearModulus() override;
		virtual void setBulkModulus() override { K0_ = bulk_modulus_; };
		virtual void setPoissonRatio() override;
		virtual void assignDerivedMaterialParameters() override;
	private:
		Real getPoissonRatio();
		Real getShearModulus();
		Real getYoungsModulus();
		Real getLambda();
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
		StdLargeVec<Matd> local_f0f0_, local_s0s0_, local_f0s0_;			/**< Sheet direction. */
		virtual void assignDerivedMaterialParameters() override 
		{
			Muscle::assignDerivedMaterialParameters();
		};
		/** initialize the local properties, fiber and sheet direction. */
		void initializeFiberAndSheet();
	public:
		StdLargeVec<Vecd> local_f0_; /**< local fiber direction. */
		StdLargeVec<Vecd> local_s0_; /**< local sheet direction. */

		LocallyOrthotropicMuscle() : Muscle()
		{
			material_name_ = "LocallyOrthotropicMuscle";
			parameters_name_ = "LocalFiberAndSheet";
		};
		virtual ~LocallyOrthotropicMuscle() {};

		virtual void assignElasticSolidParticles(ElasticSolidParticles* elastic_particles) override;
		virtual Matd MuscleFiberDirection(size_t particle_index_i) override { return local_f0f0_[particle_index_i]; };
		/** Compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		virtual void readFromXmlForLocalParameters(std::string &filefullpath) override;
	};
}
#endif //ELASTIC_SOLID_H
