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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	neighborhood_ck.h
 * @brief 	There are the classes for particle neighborhood.
 * It saves the information for carrying out inter-particle (or pair) interaction,
 * and also considered as the local configuration of the particles.
 * @author	Xiangyu Hu
 */

#ifndef NEIGHBORHOOD_CK_H
#define NEIGHBORHOOD_CK_H

#include "neighbor_method.hpp"

namespace SPH
{
template <class NeighborMethodType>
class Neighbor : public NeighborMethodType
{
  public:
    template <class SourceIdentifier, class TargetIdentifier>
    Neighbor(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier,
             DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos)
        : NeighborMethodType(source_identifier, contact_identifier),
          dv_source_pos_(dv_source_pos), dv_target_pos_(dv_target_pos){};

    class NeighborKernel : public NeighborMethodType::SmoothingKernel
    {
        using BaseKernel = typename NeighborMethodType::SmoothingKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseKernel(ex_policy, encloser),
              source_pos_(encloser.dv_source_pos_->DelegatedData(ex_policy)),
              target_pos_(encloser.dv_target_pos_->DelegatedData(ex_policy)){};

        inline Vecd vec_r_ij(UnsignedInt i, UnsignedInt j) const { return source_pos_[i] - target_pos_[j]; };
        inline Vecd e_ij(UnsignedInt i, UnsignedInt j) const { return vec_r_ij(i, j).normalized(); };
        inline Real W_ij(UnsignedInt i, UnsignedInt j) const { return BaseKernel::W(vec_r_ij(i, j)); };
        inline Real dW_ij(UnsignedInt i, UnsignedInt j) const { return BaseKernel::dW(vec_r_ij(i, j)); };

      protected:
        Vecd *source_pos_;
        Vecd *target_pos_;
    };

    class NeighborCriterion : public NeighborMethodType::NeighborCriterion
    {
        using BaseKernel = typename NeighborMethodType::NeighborCriterion;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseKernel(ex_policy, encloser, encloser.dv_source_pos_, encloser.dv_target_pos_){};

        inline bool operator()(UnsignedInt target_index, UnsignedInt source_index) const
        {
            return BaseKernel::operator()(source_index, target_index); // Note the order of indices
        };
    };

  protected:
    DiscreteVariable<Vecd> *dv_source_pos_;
    DiscreteVariable<Vecd> *dv_target_pos_;
};
} // namespace SPH
#endif // NEIGHBORHOOD_CK_H