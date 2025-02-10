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
 * @file 	derived_solid_state.h
 * @brief	tbd.
 * @author Xiangyu Hu
 */

#ifndef DERIVED_SOLID_STATE_H
#define DERIVED_SOLID_STATE_H

#include "base_general_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
class DisplacementAndPosition : public LocalDynamics
{
  public:
    explicit DisplacementAndPosition(SPHBody &sph_body);
    virtual ~DisplacementAndPosition() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
              displacement_(encloser.dv_displacement_->DelegatedData(ex_policy)){};

      protected:
        Vecd *pos_, *pos0_, *displacement_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos0_, *dv_displacement_;
};

class UpdateDisplacementFromPosition : public DisplacementAndPosition
{
  public:
    explicit UpdateDisplacementFromPosition(SPHBody &sph_body);
    virtual ~UpdateDisplacementFromPosition() {};

    class UpdateKernel : public DisplacementAndPosition::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : DisplacementAndPosition::UpdateKernel(ex_policy, encloser){};
        void update(size_t index_i, Real dt = 0.0)
        {
            displacement_[index_i] = pos_[index_i] - pos0_[index_i];
        };
    };
};

class UpdatePositionFromDisplacement : public DisplacementAndPosition
{
  public:
    explicit UpdatePositionFromDisplacement(SPHBody &sph_body);
    virtual ~UpdatePositionFromDisplacement() {};

    class UpdateKernel : public DisplacementAndPosition::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : DisplacementAndPosition::UpdateKernel(ex_policy, encloser){};
        void update(size_t index_i, Real dt = 0.0)
        {
            pos_[index_i] = pos0_[index_i] + displacement_[index_i];
        };
    };
};
} // namespace solid_dynamics
} // namespace SPH
#endif // DERIVED_SOLID_STATE_H
