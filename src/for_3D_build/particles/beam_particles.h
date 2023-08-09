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
 * @file 	beam_particles.h
 * @brief 	This is the derived class of shell particles.
 * @author	Xipeng Lyu
 */

#ifndef BEAM_PARTICLES_H
#define BEAM_PARTICLES_H

#include "base_particles.h"
#include "base_particles.hpp"
#include "elastic_solid.h"
#include "solid_particles.h"
#include "particle_generator_lattice.h"
#include "base_body.h"
//#include "solid_particles.h"
#include "solid_particles_variable.h"

#include <iterator>
namespace SPH
{
class ShellParticles;
class Solid;
class ElasticSolid;
class SPHBody;
class BarParticles : public ShellParticles
        {
              public:
                BarParticles(SPHBody &sph_body, ElasticSolid *elastic_solid);
                virtual ~BarParticles(){};

                Real width_ref_;
                StdLargeVec<Vecd> b_n_;  /**< binormal direction */
                StdLargeVec<Real> width_;
                StdLargeVec<Vecd> b_n0_; /**< initial binormal direction */

                StdLargeVec<Vecd> pseudo_b_n_;            /**< current pseudo-binormal vector */
                StdLargeVec<Vecd> dpseudo_b_n_dt_;        /**< pseudo-binormal vector change rate */
                StdLargeVec<Vecd> dpseudo_b_n_d2t_;       /**< pseudo-binormal vector second order time derivation */
                StdLargeVec<Vecd> global_b_shear_stress_; /**< global b shear stress */
                StdLargeVec<Matd> global_b_stress_;       /**<  global b stress for pair interaction */
                StdLargeVec<Matd> global_b_moment_;       /**<  global b bending moment for pair interaction */
                                                          //----------------------------------------------------------------------
                //	extra generalized coordinate and velocity in local coordinate
                //----------------------------------------------------------------------
                StdLargeVec<Vecd> rotation_b_; /**< rotation angle of the initial binormal respective to each axis */
                StdLargeVec<Vecd> angular_b_vel_;		/**< angular velocity respective to each axis */
                StdLargeVec<Vecd> dangular_b_vel_dt_; /**< angular acceleration of respective to each axis*/
                                                      //	extra deformation and deformation rate in local coordinate
                //----------------------------------------------------------------------
                StdLargeVec<Matd> F_b_bending_;     /**< bending deformation gradient	*/
                StdLargeVec<Matd> dF_b_bending_dt_; /**< bending deformation gradient change rate	*/
                //----------------------------------------------------------------------
                /** get particle volume. */
                virtual Real ParticleVolume(size_t index_i) override { return Vol_[index_i] * thickness_[index_i]*width_[index_i]; }
                /** get particle mass. */
                virtual Real ParticleMass(size_t index_i) override { return mass_[index_i] * thickness_[index_i] * width_[index_i]; }
                /** Initialize variable for shell particles. */
                virtual void initializeOtherVariables() override;
                virtual BarParticles *ThisObjectPtr() override{return this;};
        };

}
#endif // SOLID_PARTICLES_H
