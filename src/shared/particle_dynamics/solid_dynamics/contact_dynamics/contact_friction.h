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
 * @file 	contact_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid contact dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONTACT_FRICTION_H
#define CONTACT_FRICTION_H

#include "base_contact_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{

/**
 * @class PairwiseFrictionFromWall
 * @brief Damping to wall by which the wall velocity is not updated
 * and the mass of wall particle is not considered.
 * Note that, currently, this class works only when the contact
 * bodies have the same resolution.
 */
class PairwiseFrictionFromWall : public LocalDynamics, public DataDelegateContact
{
  public:
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta);
    virtual ~PairwiseFrictionFromWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        std::array<Real, MaximumNeighborhoodSize> parameter_b;

        /** Contact interaction. */
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_k = wall_vel_n_[k];
            Vecd *n_k = wall_n_[k];
            Real *Vol_k = wall_Vol_n_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            // forward sweep
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd &e_ij = contact_neighborhood.e_ij_[n];

                parameter_b[n] = eta_ * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] *
                                 Vol_[index_i] * dt / contact_neighborhood.r_ij_[n];

                // only update particle i
                Vecd vel_derivative = (vel_[index_i] - vel_k[index_j]);
                Vecd n_j = e_ij.dot(n_k[index_j]) > 0.0 ? n_k[index_j] : -1.0 * n_k[index_j];
                vel_derivative -= SMAX(Real(0), vel_derivative.dot(n_j)) * n_j;
                vel_[index_i] += parameter_b[n] * vel_derivative / (mass_[index_i] - 2.0 * parameter_b[n]);
            }
            // backward sweep
            for (size_t n = contact_neighborhood.current_size_; n != 0; --n)
            {
                size_t index_j = contact_neighborhood.j_[n - 1];
                Vecd &e_ij = contact_neighborhood.e_ij_[n];

                // only update particle i
                Vecd vel_derivative = (vel_[index_i] - vel_k[index_j]);
                Vecd n_j = e_ij.dot(n_k[index_j]) > 0.0 ? n_k[index_j] : -1.0 * n_k[index_j];
                vel_derivative -= SMAX(Real(0), vel_derivative.dot(n_j)) * n_j;
                vel_[index_i] += parameter_b[n - 1] * vel_derivative / (mass_[index_i] - 2.0 * parameter_b[n - 1]);
            }
        }
    };

  protected:
    Real eta_; /**< friction coefficient */
    Real *Vol_, *mass_;
    Vecd *vel_;
    StdVec<Real *> wall_Vol_n_;
    StdVec<Vecd *> wall_vel_n_, wall_n_;
};

} // namespace solid_dynamics
} // namespace SPH
#endif // CONTACT_FRICTION_H
