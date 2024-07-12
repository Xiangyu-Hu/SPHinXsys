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
 * @file 	inelastic_solid_dynamics.h
 * @brief 	Here, we define the algorithm classes for inelastic_solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "elastic_dynamics.h"
#include "inelastic_solid.h"

namespace SPH
{
namespace solid_dynamics
{
/**
 * @class DecomposedPlasticIntegration1stHalf
 * @brief Generalized essentially non-hourglass control formulation based on volumetric-deviatoric stress decomposition.
 */
class DecomposedPlasticIntegration1stHalf
    : public DecomposedIntegration1stHalf
{
  public:
    DecomposedPlasticIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~DecomposedPlasticIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        // including gravity and force from fluid
        Vecd force = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real r_ij = inner_neighborhood.r_ij_[n];
            Vecd e_ij = inner_neighborhood.e_ij_[n];
            Vecd pair_distance = pos_[index_i] - pos_[index_j];
            Matd pair_scaling = scaling_matrix_[index_i] + scaling_matrix_[index_j];
            Matd pair_inverse_F = 0.5 * (inverse_F_[index_i] + inverse_F_[index_j]);
            Vecd e_ij_difference = pair_inverse_F * pair_distance / r_ij - e_ij;
            Real e_ij_difference_norm = e_ij_difference.norm();

            Real limiter = SMIN(10.0 * SMAX(e_ij_difference_norm - 0.05, 0.0), 1.0);

            Vecd shear_force_ij = plastic_solid_.ShearModulus() * pair_scaling * (e_ij + limiter * e_ij_difference);
            force += mass_[index_i] * ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * e_ij + shear_force_ij) *
                     inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inv_rho0_;
        }

        force_[index_i] = force;
    };

  protected:
    PlasticSolid &plastic_solid_;
    Matd *scaling_matrix_, *inverse_F_;
};
} // namespace solid_dynamics
} // namespace SPH
