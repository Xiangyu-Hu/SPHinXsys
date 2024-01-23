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
 * @file 	contact_repulsion.h
 * @brief 	TBD.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONTACT_REPULSION_H
#define CONTACT_REPULSION_H

#include "base_contact_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{

/**
 * @class SelfContactForce
 * @brief Computing the self-contact force.
 */
class SelfContactForce : public LocalDynamics, public SolidDataInner
{
  public:
    explicit SelfContactForce(SelfSurfaceContactRelation &self_contact_relation);
    virtual ~SelfContactForce(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        Vecd vel_i = vel_[index_i];
        Real p_i = self_repulsion_density_[index_i] * solid_.ContactStiffness();

        /** Inner interaction. */
        Vecd force = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            const Vecd &e_ij = inner_neighborhood.e_ij_[n];
            Real p_star = 0.5 * (p_i + self_repulsion_density_[index_j] * solid_.ContactStiffness());
            Real impedance_p = 0.5 * contact_impedance_ * (vel_i - vel_[index_j]).dot(-e_ij);
            // force to mimic pressure
            force -= 2.0 * (p_star + impedance_p) * e_ij * Vol_i * inner_neighborhood.dW_ijV_j_[n];
        }
        force_prior_[index_i] += force;
    };

  protected:
    Solid &solid_;
    StdLargeVec<Real> &mass_, &self_repulsion_density_, &Vol_;
    StdLargeVec<Vecd> &force_prior_, &vel_;
    Real contact_impedance_;
};

/**
 * @class ContactForce
 * @brief Computing the contact force.
 */
class ContactForce : public LocalDynamics, public ContactDynamicsData
{
  public:
    explicit ContactForce(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForce(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        Real p_i = repulsion_density_[index_i] * solid_.ContactStiffness();
        /** Contact interaction. */
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &contact_density_k = *(contact_contact_density_[k]);
            Solid *solid_k = contact_solids_[k];

            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];

                Real p_star = 0.5 * (p_i + contact_density_k[index_j] * solid_k->ContactStiffness());
                // force due to pressure
                force -= 2.0 * p_star * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
            }
        }
        force_prior_[index_i] += force;
    };

  protected:
    Solid &solid_;
    StdLargeVec<Real> &repulsion_density_, &Vol_, &mass_;
    StdLargeVec<Vecd> &force_prior_;
    StdVec<Solid *> contact_solids_;
    StdVec<StdLargeVec<Real> *> contact_contact_density_;
};

/**
 * @class ContactForceFromWall
 * @brief Computing the contact force from a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class ContactForceFromWall : public LocalDynamics, public ContactWithWallData
{
  public:
    explicit ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForceFromWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        Real p_i = repulsion_density_[index_i] * solid_.ContactStiffness();
        /** Contact interaction. */
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                Vecd e_ij = contact_neighborhood.e_ij_[n];

                // force due to pressure
                force -= 2.0 * p_i * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
            }
        }
        force_prior_[index_i] += force;
    };

  protected:
    Solid &solid_;
    StdLargeVec<Real> &repulsion_density_, &Vol_, &mass_;
    StdLargeVec<Vecd> &force_prior_;
};

/**
 * @class ContactForceToWall
 * @brief Computing contact force acting on a rigid wall.
 */
class ContactForceToWall : public LocalDynamics, public ContactDynamicsData
{
  public:
    explicit ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForceToWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        /** Contact interaction. */
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &contact_density_k = *(contact_contact_density_[k]);
            Solid *solid_k = contact_solids_[k];

            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];

                Real p_star = contact_density_k[index_j] * solid_k->ContactStiffness();
                // force due to pressure
                force -= 2.0 * p_star * e_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n];
            }
        }
        force_prior_[index_i] += force;
    };

  protected:
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<Vecd> &force_prior_;
    StdVec<Solid *> contact_solids_;
    StdVec<StdLargeVec<Real> *> contact_contact_density_;
};

/**
 * @class DynamicContactForceWithWall
 * @brief Computing the contact force with a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class DynamicContactForceWithWall : public LocalDynamics, public ContactDynamicsData
{
  public:
    explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0);
    virtual ~DynamicContactForceWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        Vecd vel_i = vel_[index_i];

        /** Contact interaction. */
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
            Real particle_spacing_ratio2 =
                1.0 / (sph_body_.sph_adaptation_->ReferenceSpacing() * particle_spacing_j1);
            particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

            StdLargeVec<Vecd> &n_k = *(contact_n_[k]);
            StdLargeVec<Vecd> &vel_n_k = *(contact_vel_[k]);

            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Vecd n_k_j = n_k[index_j];

                Real impedance_p = 0.5 * impedance_ * (vel_i - vel_n_k[index_j]).dot(-n_k_j);
                Real overlap = contact_neighborhood.r_ij_[n] * n_k_j.dot(e_ij);
                Real delta = 2.0 * overlap * particle_spacing_j1;
                Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
                Real penalty_p = penalty_strength_ * beta * fabs(overlap) * reference_pressure_;

                // force due to pressure
                force -= 2.0 * (impedance_p + penalty_p) * e_ij.dot(n_k_j) *
                         n_k_j * Vol_i * contact_neighborhood.dW_ijV_j_[n];
            }
        }

        force_prior_[index_i] += force;
    };

  protected:
    Solid &solid_;
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<Vecd> &vel_, &force_prior_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_, contact_n_;
    Real penalty_strength_;
    Real impedance_, reference_pressure_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // CONTACT_REPULSION_H
