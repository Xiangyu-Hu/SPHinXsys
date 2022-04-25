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

namespace SPH
{

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
		Real E0_;  /*< Youngs or tensile modules  */
		Real G0_;  /*< shearmodules  */
		Real K0_;  /*< bulkmodules  */
		Real nu_;  /*< Poisson ratio  */
		Real c0_;  /*< sound wave speed */
		Real ct0_; /*< tensile wave speed */
		Real cs0_; /*< shear wave speed */
		ElasticSolidParticles *elastic_particles_;

		void setSoundSpeeds();

	public:
		explicit ElasticSolid(Real rho0)
			: Solid(rho0), c0_(0.0), ct0_(0.0), cs0_(0.0),
			  E0_(0.0), G0_(0.0), K0_(0.0), nu_(0.0), elastic_particles_(nullptr)
		{
			material_type_ = "ElasticSolid";
		};
		virtual ~ElasticSolid(){};

		virtual void assignElasticSolidParticles(ElasticSolidParticles *elastic_particles);
		Real ReferenceSoundSpeed() { return c0_; };
		Real TensileWaveSpeed() { return ct0_; };
		Real ShearWaveSpeed() { return cs0_; };
		Real YoungsModulus() { return E0_; };
		Real ShearModulus() { return G0_; };
		Real BulkModulus() { return K0_; };
		Real PoissonRatio() { return nu_; };

		/** Compute the stress through deformation, which can be Green-Lagrangian tensor, left or right Cauchy tensor. */
		virtual Matd ConstitutiveRelation(Matd &deformation, size_t particle_index_i) = 0;
		/** Compute the Cauchy stress through Eulerian Almansi strain tensor. */
		virtual Matd EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i) = 0;
		/** Compute numerical damping stress using right Cauchy tensor. */
		virtual Matd NumericalDampingRightCauchy(Matd &deformation, Matd &deformation_rate, Real smoothing_length, size_t particle_index_i);
		/** Compute numerical damping stress using left Cauchy tensor. */
		virtual Matd NumericalDampingLeftCauchy(Matd &deformation, Matd &deformation_rate, Real smoothing_length, size_t particle_index_i);
		/** Numerical demaping is computed between particles i and j */
		virtual Real PairNumericalDamping(Real dE_dt_ij, Real smoothing_length);

		/** Deviatoric Kirchhoff stress related with the deviatoric part of left Cauchy-Green deformation tensor.
		 *  Note that, dependent of the normalizeation of the later, the returned stress can be normalized or non-normalized. */
		virtual Matd DeviatoricKirchhoff(const Matd &deviatoric_be);
		/** Volumetric Kirchhoff stress from determinate */
		virtual Real VolumetricKirchhoff(Real J) = 0;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() = 0;

		virtual ElasticSolid *ThisObjectPtr() override { return this; };
	};

	/**
	* @class LinearElasticSolid
	* @brief Isotropic linear elastic solid.
	* Note that only basic parameters are used to set ElasticSolid parmaters 
	*/
	class LinearElasticSolid : public ElasticSolid
	{
	public:
		explicit LinearElasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
		virtual ~LinearElasticSolid(){};

		virtual void assignElasticMaterialParameters(Real youngs_modulus, Real poisson_ratio);

		virtual Matd ConstitutiveRelation(Matd &deformation, size_t particle_index_i) override;
		virtual Matd EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress from determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() override { return "PK2"; };

		/** get methods */
		Real getYoungsModulus() { return E0_; };
		Real getPoissonRatio() { return nu_; };
		Real getDensity() { return rho0_; };

	protected:
		Real lambda0_; /*< first Lame parameter */
		Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
		Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
		Real getLambda(Real youngs_modulus, Real poisson_ratio);
	};

	/**
	* @class NeoHookeanSolid
	* @brief Neo-Hookean solid, Compressible formulation!
	*/
	class NeoHookeanSolid : public LinearElasticSolid
	{
	public:
		explicit NeoHookeanSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
			: LinearElasticSolid(rho0, youngs_modulus, poisson_ratio)
		{
			material_type_ = "NeoHookeanSolid";
		};
		virtual ~NeoHookeanSolid(){};

		/** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
		virtual Matd ConstitutiveRelation(Matd &deformation, size_t particle_index_i) override;
		virtual Matd EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress from determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
	};

	/**
	* @class NeoHookeanSolidIncompressible
	* @brief Neo-Hookean solid, Incomressible formulation!
	* Currently only works with KirchhoffStressRelaxationFirstHalf, not with StressRelaxationFirstHalf
	*/
	class NeoHookeanSolidIncompressible : public LinearElasticSolid
	{
	public:
		NeoHookeanSolidIncompressible(Real rho_0, Real Youngs_modulus, Real poisson)
			: LinearElasticSolid(rho_0, Youngs_modulus, poisson)
		{
			material_type_ = "NeoHookeanSolidIncompressible";
		};
		virtual ~NeoHookeanSolidIncompressible() {};
	
		/** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
		virtual Matd ConstitutiveRelation(Matd &deformation, size_t particle_index_i) override;
		virtual Matd EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress from determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
	};

	/**
	* @class OrthotropicSolid
	* @brief Ortothropic solid - generic definition with 3 orthogonal directions + 9 independent parameters, ONLY for 3D applications
	* @param "a" --> 3 principal direction vectors
	* @param "E" --> 3 principal Young's moduli
	* @param "G" --> 3 principal shear moduli
	* @param "poisson" --> 3 principal Poisson's ratios
	*/
	class OrthotropicSolid : public LinearElasticSolid
	{
	public:
		OrthotropicSolid(Real rho_0, std::array<Vecd, 3> a, std::array<Real, 3> E, std::array<Real, 3> G, std::array<Real, 3> poisson);

		/** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress determinate */
		virtual Real VolumetricKirchhoff(Real J) override;

	protected:
		// input data
		std::array<Vecd, 3> a_;
		std::array<Real, 3> E_;
		std::array<Real, 3> G_;
		std::array<Real, 3> poisson_;
		// calculated data
		Real Mu_[3];
		Matd Lambda_;
		Matd A_[3];

		virtual void CalculateAllMu();
		virtual void CalculateAllLambda();
		virtual void CalculateA0();
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
		explicit FeneNeoHookeanSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
			: LinearElasticSolid(rho0, youngs_modulus, poisson_ratio), j1_m_(1.0)
		{
			material_type_ = "FeneNeoHookeanSolid";
		};
		virtual ~FeneNeoHookeanSolid() {};
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
	};

	/**
	* @class Muscle
	* @brief Globally orthotropic muscle. 
	*/
	class Muscle : public NeoHookeanSolid
	{
	public:
		explicit Muscle(Real rho0, Real bulk_modulus,
						const Vecd &f0, const Vecd &s0, const Real (&a0)[4], const Real (&b0)[4])
			: NeoHookeanSolid(rho0, this->getYoungsModulus(bulk_modulus, a0, b0), this->getPoissonRatio(bulk_modulus, a0, b0)),
			  f0_(f0), s0_(s0), f0f0_(SimTK::outer(f0_, f0_)), s0s0_(SimTK::outer(s0_, s0_)),
			  f0s0_(SimTK::outer(f0_, s0_))
		{
			material_type_ = "Muscle";
			std::copy(a0, a0 + 4, a0_);
			std::copy(b0, b0 + 4, b0_);
		};
		virtual ~Muscle(){};

		virtual Matd MuscleFiberDirection(size_t particle_index_i) { return f0f0_; };
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd &deformation, size_t particle_index_i) override;
		/** Volumetric Kirchhoff stress form determinate */
		virtual Real VolumetricKirchhoff(Real J) override;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };

		virtual Muscle *ThisObjectPtr() override { return this; };

	protected:
		Vecd f0_, s0_;			  /**< Reference fiber and sheet directions as basic parameter. */
		Matd f0f0_, s0s0_, f0s0_; /**< Tensor products of fiber and sheet directions as basic parameter.. */
		Real a0_[4], b0_[4];	  /**< constitutive parameters  as basic parameter.*/

	private:
		Real getPoissonRatio(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4]);
		Real getShearModulus(const Real (&a0)[4], const Real (&b0)[4]);
		Real getYoungsModulus(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4]);
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
		StdLargeVec<Matd> local_f0f0_, local_s0s0_, local_f0s0_; /**< Sheet direction. */

		/** initialize the local properties, fiber and sheet direction. */
		void initializeFiberAndSheet();

	public:
		StdLargeVec<Vecd> local_f0_; /**< local fiber direction. */
		StdLargeVec<Vecd> local_s0_; /**< local sheet direction. */

		explicit LocallyOrthotropicMuscle(Real rho0, Real bulk_modulus,
						const Vecd &f0, const Vecd &s0, const Real (&a0)[4], const Real (&b0)[4])
			: Muscle(rho0, bulk_modulus, f0, s0, a0, b0)
		{
			material_type_ = "LocallyOrthotropicMuscle";
		};
		virtual ~LocallyOrthotropicMuscle(){};

		virtual void assignElasticSolidParticles(ElasticSolidParticles *elastic_particles) override;
		virtual Matd MuscleFiberDirection(size_t particle_index_i) override { return local_f0f0_[particle_index_i]; };
		/** Compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd& deformation, size_t particle_index_i) override;
		virtual void readFromXmlForLocalParameters(const std::string &filefullpath) override;
		/** Define the calculation of the stress matrix for postprocessing */
		virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
	};
}
#endif //ELASTIC_SOLID_H
