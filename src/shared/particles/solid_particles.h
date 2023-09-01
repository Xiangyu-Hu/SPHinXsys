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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	solid_particles.h
 * @brief 	This is the derived class of base particles.
 * @author	Chi Zhang, Dong Wu and Xiangyu Hu
 */

#ifndef SOLID_PARTICLES_H
#define SOLID_PARTICLES_H

#include "base_particles.h"
#include "base_particles.hpp"
#include "elastic_solid.h"

#include "particle_generator_lattice.h"
namespace SPH
{
/**
 *	pre-claimed classes
 */
class Solid;
class ElasticSolid;

/**
 * @class SolidParticles
 * @brief A group of particles with solid body particle data.
 */
class SolidParticles : public BaseParticles
{
  public:
    SolidParticles(SPHBody &sph_body, Solid *solid);
    virtual ~SolidParticles(){};

    StdLargeVec<Vecd> pos0_; /**< initial position */
    StdLargeVec<Vecd> n_;    /**< normal direction */
    StdLargeVec<Vecd> n0_;   /**< initial normal direction */
    StdLargeVec<Matd> B_;    /**< configuration correction for linear reproducing */
    Solid &solid_;

    /** Get wall average velocity when interacting with fluid. */
    virtual StdLargeVec<Vecd> *AverageVelocity() { return &vel_; };
    /** Get wall average acceleration when interacting with fluid. */
    virtual StdLargeVec<Vecd> *AverageAcceleration() { return &acc_; };
    /** Initialized variables for solid particles. */
    virtual void initializeOtherVariables() override;
    /** Return this pointer. */
    virtual SolidParticles *ThisObjectPtr() override { return this; };
};

/**
 * @class ElasticSolidParticles
 * @brief A group of particles with elastic body particle data.
 */
class ElasticSolidParticles : public SolidParticles
{
  public:
    ElasticSolidParticles(SPHBody &sph_body, ElasticSolid *elastic_solid);
    virtual ~ElasticSolidParticles(){};

    StdLargeVec<Matd> F_;     /**<  deformation tensor */
    StdLargeVec<Matd> dF_dt_; /**<  deformation tensor change rate */
    ElasticSolid &elastic_solid_;
    //----------------------------------------------------------------------
    //		for fluid-structure interaction (FSI)
    //----------------------------------------------------------------------
    StdLargeVec<Vecd> vel_ave_; /**<  fluid time-step averaged particle velocity */
    StdLargeVec<Vecd> acc_ave_; /**<  fluid time-step averaged particle acceleration */

    /** Return the Lagrange strain. */
    Matd getGreenLagrangeStrain(size_t particle_i);
    /** Computing principal strain - returns the principal strains in descending order (starting from the largest) */
    Vecd getPrincipalStrains(size_t particle_i);
    /** Computing von Mises equivalent strain from a static (constant) formulation. */
    Real getVonMisesStrain(size_t particle_i);
    /** Computing von Mises equivalent strain from a "dynamic" formulation. This depends on the Poisson's ratio (from commercial FEM software Help). */
    Real getVonMisesStrainDynamic(size_t particle_i, Real poisson);
    /** Computing von Mises strain for all particles. - "static" or "dynamic"*/
    StdLargeVec<Real> getVonMisesStrainVector(std::string strain_measure = "static");
    /** Computing maximum von Mises strain from all particles. - "static" or "dynamic" */
    Real getVonMisesStrainMax(std::string strain_measure = "static");
    /** Return the max principal strain. */
    Real getPrincipalStrainMax();
    /** get the Cauchy stress. */
    Matd getStressCauchy(size_t particle_i);
    /** get the PK2 stress. */
    Matd getStressPK2(size_t particle_i);
    /** Computing principal_stresses - returns the principal stresses in descending order (starting from the largest) */
    Vecd getPrincipalStresses(size_t particle_i);
    /** Computing von_Mises_stress - "Cauchy" or "PK2" decided based on the stress_measure_ */
    Real getVonMisesStress(size_t particle_i);
    /** Computing von Mises stress for all particles. - "Cauchy" or "PK2" decided based on the stress_measure_ */
    StdLargeVec<Real> getVonMisesStressVector();
    /** Computing maximum von Mises stress from all particles. - "Cauchy" or "PK2" decided based on the stress_measure_ */
    Real getVonMisesStressMax();
    Real getPrincipalStressMax();

    /** Computing displacement. */
    Vecd displacement(size_t particle_i);
    /** Return the displacement. */
    StdLargeVec<Vecd> getDisplacement();
    /** get the max displacement. */
    Real getMaxDisplacement();

    /**< Computing normal vector. */
    Vecd normal(size_t particle_i);
    /** get the normal vector. */
    StdLargeVec<Vecd> getNormal();

    /** relevant stress measure */
    std::string stress_measure_;

    /** Get wall average velocity when interacting with fluid. */
    virtual StdLargeVec<Vecd> *AverageVelocity() override { return &vel_ave_; };
    /** Get wall average acceleration when interacting with fluid. */
    virtual StdLargeVec<Vecd> *AverageAcceleration() override { return &acc_ave_; };

    /** Initialize the variables for elastic particle. */
    virtual void initializeOtherVariables() override;
    /** Return this pointer. */
    virtual ElasticSolidParticles *ThisObjectPtr() override { return this; };
};

/**
 * @class ShellParticles
 * @brief A group of particles with shell particle data.
 */
class ShellParticles : public ElasticSolidParticles
{
  public:
    ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid);
    virtual ~ShellParticles(){};

    Real thickness_ref_;                      /**< Shell thickness. */
    StdLargeVec<Matd> transformation_matrix_; /**< initial transformation matrix from global to local coordinates */
    StdLargeVec<Real> thickness_;             /**< shell thickness */
    /**
     *	extra generalized coordinates in global coordinate
     */
    StdLargeVec<Vecd> pseudo_n_;      /**< current pseudo-normal vector */
    StdLargeVec<Vecd> dpseudo_n_dt_;  /**< pseudo-normal vector change rate */
    StdLargeVec<Vecd> dpseudo_n_d2t_; /**< pseudo-normal vector second order time derivation */
    /**
     *	extra generalized coordinate and velocity in local coordinate
     */
    StdLargeVec<Vecd> rotation_;        /**< rotation angle of the initial normal respective to each axis */
    StdLargeVec<Vecd> angular_vel_;     /**< angular velocity respective to each axis */
    StdLargeVec<Vecd> dangular_vel_dt_; /**< angular acceleration of respective to each axis*/
    /**
     *	extra deformation and deformation rate in local coordinate
     */
    StdLargeVec<Matd> F_bending_;     /**< bending deformation gradient	*/
    StdLargeVec<Matd> dF_bending_dt_; /**< bending deformation gradient change rate	*/
    /**
     *	extra stress for pair interaction in global coordinate
     */
    StdLargeVec<Vecd> global_shear_stress_; /**< global shear stress */
    StdLargeVec<Matd> global_stress_;       /**<  global stress for pair interaction */
    StdLargeVec<Matd> global_moment_;       /**<  global bending moment for pair interaction */
    /**
     *	extra stress for calculating von Mises stress of shell mid-surface
     */
    StdLargeVec<Matd> mid_surface_cauchy_stress_;
    /**
     *	extra scaling matrix fot numerical damping
     */
    StdLargeVec<Matd> numerical_damping_scaling_;

    /** get particle volume. */
    virtual Real ParticleVolume(size_t index_i) override { return Vol_[index_i] * thickness_[index_i]; }
    /** get particle mass. */
    virtual Real ParticleMass(size_t index_i) override { return mass_[index_i] * thickness_[index_i]; }
    /** Initialize variable for shell particles. */
    virtual void initializeOtherVariables() override;
    /** Return this pointer. */
    virtual ShellParticles *ThisObjectPtr() override { return this; };
};



}
#endif // SOLID_PARTICLES_H
