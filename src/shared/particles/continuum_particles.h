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
 * @file 	continuum_particles.h
 * @brief 	This is the derived class of base particles.
 * @author	Shuaihao Zhang and Xiangyu Hu
 */
#ifndef CONTINUUM_PARTICLES_H
#define CONTINUUM_PARTICLES_H

#include "base_particles.hpp"
#include "general_continuum.h"
namespace SPH
{
class GeneralContinuum;

class ContinuumParticles : public BaseParticles
{
  public:
    StdLargeVec<Matd> strain_tensor_;
    StdLargeVec<Matd> strain_tensor_rate_;
    StdLargeVec<Vecd> acc_shear_;

    StdLargeVec<Matd> shear_stress_;
    StdLargeVec<Matd> shear_stress_rate_;
    StdLargeVec<Matd> velocity_gradient_;

    StdLargeVec<Real> von_mises_stress_;
    StdLargeVec<Real> von_mises_strain_;

    StdLargeVec<Vecd> pos0_; /**< initial position */
    StdLargeVec<Vecd> n_;    /**<  current normal direction */
    StdLargeVec<Vecd> n0_;   /**<  initial normal direction */

    GeneralContinuum &continuum_;

    ContinuumParticles(SPHBody &sph_body, GeneralContinuum *continuum);
    virtual ~ContinuumParticles(){};

    virtual void initializeOtherVariables() override;
    virtual ContinuumParticles *ThisObjectPtr() override { return this; };
};

class PlasticContinuumParticles : public ContinuumParticles
{
  public:
    StdLargeVec<Mat3d> elastic_strain_tensor_3D_;
    StdLargeVec<Mat3d> elastic_strain_rate_3D_;

    StdLargeVec<Mat3d> strain_tensor_3D_;
    StdLargeVec<Mat3d> stress_tensor_3D_;
    StdLargeVec<Mat3d> strain_rate_3D_;
    StdLargeVec<Mat3d> stress_rate_3D_;

    StdLargeVec<Real> vertical_stress_;
    StdLargeVec<Real> acc_deviatoric_plastic_strain_;

    Real getDeviatoricPlasticStrain(Mat3d &strain_tensor);

    PlasticContinuum &plastic_continuum_;

    PlasticContinuumParticles(SPHBody &sph_body, PlasticContinuum *plastic_continuum);
    virtual ~PlasticContinuumParticles(){};

    virtual void initializeOtherVariables() override;
    virtual ContinuumParticles *ThisObjectPtr() override { return this; };
};
} // namespace SPH

#endif // CONTINUUM_PARTICLES_H