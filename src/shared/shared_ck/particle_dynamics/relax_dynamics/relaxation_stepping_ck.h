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
 * @file relaxation_stepping_ck.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef RELAXATION_STEPPING_CK_H
#define RELAXATION_STEPPING_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
class RelaxationScalingCK : public LocalDynamicsReduce<ReduceMax>
{
  public:
    RelaxationScalingCK(SPHBody &sph_body);
    virtual ~RelaxationScalingCK() {};

    class FinishDynamics
    {
        Real h_ref_;

      public:
        using OutputType = Real;
        FinishDynamics(RelaxationScalingCK &encloser);
        Real Result(Real reduced_value);
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy>
        ReduceKernel(const ExecutionPolicy &ex_policy, RelaxationScalingCK &encloser)
            : residue_(encloser.dv_residue_->DelegatedData(ex_policy)),
              h_ref_(encloser.h_ref_){};

        Real reduce(size_t index_i, Real dt)
        {
            return residue_[index_i].norm();
        };

      protected:
        Vecd *residue_;
        Real h_ref_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_residue_;
    Real h_ref_;
};

class PositionRelaxationCK : public LocalDynamics
{
  public:
    explicit PositionRelaxationCK(SPHBody &sph_body);
    virtual ~PositionRelaxationCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, PositionRelaxationCK &encloser)
            : pos_(encloser.pos_->DelegatedData(ex_policy)),
              residue_(encloser.residue_->DelegatedData(ex_policy)){};

        void update(size_t index_i, Real dt_square)
        {
            pos_[index_i] += residue_[index_i] * dt_square * 0.5;
        };

      protected:
        Vecd *pos_, *residue_;
    };

  protected:
    DiscreteVariable<Vecd> *pos_, *residue_;
};

} // namespace SPH
#endif // RELAXATION_STEPPING_CK_H
