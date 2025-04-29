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
 * @file neighbor_method.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef NEIGHBOR_METHOD_H
#define NEIGHBOR_METHOD_H

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
class ConstantSmoothingLength
{
  public:
    template <class DynamicsIdentifier>
    ConstantSmoothingLength(DynamicsIdentifier &identifier)
        : inv_h_(1.0 / identifier.getSPHAdaptation().ReferenceSmoothingLength()){};

    template <class SourceIdentifier, class TargetIdentifier>
    ConstantSmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : inv_h_(1.0 / SMAX(source_identifier.getSPHAdaptation().ReferenceSmoothingLength(),
                            contact_identifier.getSPHAdaptation().ReferenceSmoothingLength())){};

    class ComputingKernel
    {
        Real inv_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, ConstantSmoothingLength &smoothing_length)
            : inv_h_(smoothing_length.inv_h_){};
        Real
        operator()(UnsignedInt i, UnsignedInt j) const { return inv_h_; };
    };

  protected:
    Real inv_h_;
};

class VariableSmoothingLength
{
  public:
    template <class DynamicsIdentifier>
    VariableSmoothingLength(DynamicsIdentifier &identifier)
        : dv_source_h_(identifier.getBaseParticle().template getVariableByName<Real>("SmoothingLength")),
          dv_target_h_(dv_source_h_){};

    template <class SourceIdentifier, class TargetIdentifier>
    VariableSmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : dv_source_h_(source_identifier.getBaseParticle().template getVariableByName<Real>("SmoothingLength")),
          dv_target_h_(contact_identifier.getBaseParticle().template getVariableByName<Real>("SmoothingLength")){};

    class ComputingKernel
    {
        Real *source_h_;
        Real *target_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, VariableSmoothingLength &smoothing_length)
            : source_h_(smoothing_length.dv_source_h_->DelegatedData(ex_policy)),
              target_h_(smoothing_length.dv_target_h_->DelegatedData(ex_policy)){};
        Real
        operator()(UnsignedInt i, UnsignedInt j) const { return 1.0 / SMAX(source_h_[i], target_h_[j]); };
    };

  protected:
    DiscreteVariable<Real> *dv_source_h_;
    DiscreteVariable<Real> *dv_target_h_;
};

} // namespace SPH
#endif // NEIGHBOR_METHOD_H