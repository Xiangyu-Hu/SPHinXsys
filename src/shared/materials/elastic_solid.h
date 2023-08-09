/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	elastic_solid.h
 * @brief 	These are classes for define properties of elastic solid materials.
 *			These classes are based on isotropic linear elastic solid.
 * 			Several more complex materials, including neo-Hookean, FENE neo-Hookean
 *			and anisotropic muscle, are derived from the basic elastic solid class.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ELASTIC_SOLID_H
#define ELASTIC_SOLID_H

#include "base_material.h"
#include <fstream>

namespace SPH
{
/**
 * @class ElasticSolid
 * @brief Abstract class for elastic solid material.
 */
class ElasticSolid : public Solid
{
  protected:
    Real c0_;  /*< sound wave speed */
    Real ct0_; /*< tensile wave speed */
    Real cs0_; /*< shear wave speed */
    Real E0_;  /*< Youngs or tensile modules  */
    Real G0_;  /*< shear modules  */
    Real K0_;  /*< bulk modules  */
    Real nu_;  /*< Poisson ratio  */

    void setSoundSpeeds();

  public:
    explicit ElasticSolid(Real rho0) : Solid(rho0), c0_(0.0), ct0_(0.0), cs0_(0.0),
                                       E0_(0.0), G0_(0.0), K0_(0.0), nu_(0.0)
    {
        material_type_name_ = "ElasticSolid";
    };
    virtual ~ElasticSolid(){};

    Real ReferenceSoundSpeed() { return c0_; };
    Real TensileWaveSpeed() { return ct0_; };
    Real ShearWaveSpeed() { return cs0_; };
    Real YoungsModulus() { return E0_; };
    Real ShearModulus() { return G0_; };
    Real BulkModulus() { return K0_; };
    Real PoissonRatio() { return nu_; };

    /** 1st Piola-Kirchhoff stress through deformation. */
    virtual Matd StressPK1(Matd &deformation, size_t particle_index_i) = 0;
    /** 2nd Piola-Kirchhoff stress through deformation. */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) = 0;
    /** Cauchy stress through Eulerian Almansi strain tensor. */
    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t particle_index_i) = 0;
    /** Numerical damping stress using right Cauchy tensor. */
    template <typename ScalingType>
    Matd NumericalDampingRightCauchy(const Matd &deformation, const Matd &deformation_rate, const ScalingType &scaling, size_t particle_index_i)
    {
        Matd strain_rate = 0.5 * (deformation_rate.transpose() * deformation + deformation.transpose() * deformation_rate);
        Matd normal_rate = getDiagonal(strain_rate);
        return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * scaling;
    }
    /** Numerical damping stress using left Cauchy tensor. */
    template <typename ScalingType>
    Matd NumericalDampingLeftCauchy(const Matd &deformation, const Matd &deformation_rate, const ScalingType &scaling, size_t particle_index_i)
    {
        Matd strain_rate = 0.5 * (deformation_rate * deformation.transpose() + deformation * deformation_rate.transpose());
        Matd normal_rate = getDiagonal(strain_rate);
        return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * scaling;
    }
    /** Numerical damping is computed between particles i and j */
    virtual Real PairNumericalDamping(Real dE_dt_ij, Real smoothing_length);

    /** Deviatoric Kirchhoff stress related with the deviatoric part of left Cauchy-Green deformation tensor.
     *  Note that, dependent of the normalization of the later, the returned stress can be normalized or non-normalized. */
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
 * 		 Note that only basic parameters are used to set ElasticSolid parameters
 */
class LinearElasticSolid : public ElasticSolid
{
  public:
    explicit LinearElasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
    virtual ~LinearElasticSolid(){};

    virtual Matd StressPK1(Matd &deformation, size_t particle_index_i) override;
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
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
 * @class SaintVenantKirchhoffSolid
 * @brief Every thing same as linear elastic but assume
 * 		 finite deformation (geometry nonlinearity).
 * 		 It is the simplest hyper-elastic material model.
 */
class SaintVenantKirchhoffSolid : public LinearElasticSolid
{
  public:
    explicit SaintVenantKirchhoffSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio)
    {
        material_type_name_ = "SaintVenantKirchhoffSolid";
    };
    virtual ~SaintVenantKirchhoffSolid(){};

    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
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
        material_type_name_ = "NeoHookeanSolid";
    };
    virtual ~NeoHookeanSolid(){};

    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
    /** Volumetric Kirchhoff stress from determinate */
    virtual Real VolumetricKirchhoff(Real J) override;
    /** Define the calculation of the stress matrix for postprocessing */
    virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
};

/**
 * @class NeoHookeanSolidIncompressible
 * @brief Neo-Hookean solid, Incompressible formulation!
 * Currently only works with DecomposedIntegration1stHalf, not with Integration1stHalf
 */
class NeoHookeanSolidIncompressible : public LinearElasticSolid
{
  public:
    NeoHookeanSolidIncompressible(Real rho_0, Real Youngs_modulus, Real poisson)
        : LinearElasticSolid(rho_0, Youngs_modulus, poisson)
    {
        material_type_name_ = "NeoHookeanSolidIncompressible";
    };
    virtual ~NeoHookeanSolidIncompressible(){};

    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    virtual Matd StressCauchy(Matd &almansi_strain, Matd &F, size_t particle_index_i) override;
    /** Volumetric Kirchhoff stress from determinate */
    virtual Real VolumetricKirchhoff(Real J) override;
};

/**
 * @class OrthotropicSolid
 * @brief Generic definition with 3 orthogonal directions + 9 independent parameters,
 * ONLY for 3D applications
 * @param "a" --> 3 principal direction vectors
 * @param "E" --> 3 principal Young's moduli
 * @param "G" --> 3 principal shear moduli
 * @param "poisson" --> 3 principal Poisson's ratios
 */
class OrthotropicSolid : public LinearElasticSolid
{
  public:
    OrthotropicSolid(Real rho_0, std::array<Vecd, Dimensions> a, std::array<Real, Dimensions> E,
                     std::array<Real, Dimensions> G, std::array<Real, Dimensions> poisson);

    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    /** Volumetric Kirchhoff stress determinate */
    virtual Real VolumetricKirchhoff(Real J) override;

  protected:
    // input data
    std::array<Vecd, Dimensions> a_;
    std::array<Real, Dimensions> E_;
    std::array<Real, Dimensions> G_;
    std::array<Real, Dimensions> poisson_;
    // calculated data
    Real Mu_[Dimensions];
    Matd Lambda_;
    Matd A_[Dimensions];

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
    Real j1_m_; /**< reference extension as basic parameter */
  public:
    explicit FeneNeoHookeanSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio), j1_m_(1.0)
    {
        material_type_name_ = "FeneNeoHookeanSolid";
    };
    virtual ~FeneNeoHookeanSolid(){};
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
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
        : NeoHookeanSolid(rho0, this->getYoungsModulus(bulk_modulus, a0, b0),
                          this->getPoissonRatio(bulk_modulus, a0, b0)),
          f0_(f0), s0_(s0), f0f0_(f0_ * f0_.transpose()), s0s0_(s0_ * s0_.transpose()),
          f0s0_(f0_ * s0_.transpose() + s0_ * f0_.transpose())
    {
        material_type_name_ = "Muscle";
        std::copy(a0, a0 + 4, a0_);
        std::copy(b0, b0 + 4, b0_);
    };
    virtual ~Muscle(){};

    virtual Matd MuscleFiberDirection(size_t particle_index_i) { return f0f0_; };
    /** compute the stress through Constitutive relation. */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    /** Volumetric Kirchhoff stress form determinate */
    virtual Real VolumetricKirchhoff(Real J) override;
    /** Define the calculation of the stress matrix for postprocessing */
    virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };

    virtual Muscle *ThisObjectPtr() override { return this; };

  protected:
    Vecd f0_, s0_;            /**< Reference fiber and sheet directions as basic parameter. */
    Matd f0f0_, s0s0_, f0s0_; /**< Tensor products of fiber and sheet directions as basic parameter.. */
    Real a0_[4], b0_[4];      /**< constitutive parameters  as basic parameter.*/

  private:
    Real getPoissonRatio(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4]);
    Real getShearModulus(const Real (&a0)[4], const Real (&b0)[4]);
    Real getYoungsModulus(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4]);
};

/**
 * @class LocallyOrthotropicMuscle
 * @brief muscle model is a anisotropic material in which
 * 		 there are local fiber direction and cross-fiber sheet direction.
 * 		 the model here is from
 * 		 Holzapfel and Ogden, 2009, Phil. Trans. R. Soc. 367:3445-3475
 * 		 we consider a neo-Hookean model for the background isotropic contribution.
 */
class LocallyOrthotropicMuscle : public Muscle
{
  protected:
    StdLargeVec<Matd> local_f0f0_, local_s0s0_, local_f0s0_; /**< Sheet direction. */

  public:
    StdLargeVec<Vecd> local_f0_; /**< local fiber direction. */
    StdLargeVec<Vecd> local_s0_; /**< local sheet direction. */

    explicit LocallyOrthotropicMuscle(Real rho0, Real bulk_modulus,
                                      const Vecd &f0, const Vecd &s0, const Real (&a0)[4], const Real (&b0)[4])
        : Muscle(rho0, bulk_modulus, f0, s0, a0, b0)
    {
        material_type_name_ = "LocallyOrthotropicMuscle";
    };
    virtual ~LocallyOrthotropicMuscle(){};

    virtual void registerReloadLocalParameters(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Matd MuscleFiberDirection(size_t particle_index_i) override { return local_f0f0_[particle_index_i]; };
    /** Compute the stress through Constitutive relation. */
    virtual Matd StressPK2(Matd &deformation, size_t particle_index_i) override;
    /** Define the calculation of the stress matrix for postprocessing */
    virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
};
} // namespace SPH
#endif // ELASTIC_SOLID_H
