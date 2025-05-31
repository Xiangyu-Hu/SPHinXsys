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
 * @file    adapt_indication.h
 * @brief   to indicate the adaptation level
 * @author	Xiangyu Hu
 */

#ifndef ADAPT_INDICATION_H
#define ADAPT_INDICATION_H

#include "adapt_criterion.h"
#include "base_general_dynamics.h"

namespace SPH
{
template <typename IndicationCriterion>
class AdaptLevelIndication : public LocalDynamics
{
    using IndicationKernel = typename IndicationCriterion::ComputingKernel;

  public:
    explicit AdaptLevelIndication(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          indication_method_(sph_adaptation_, particles_),
          dv_adapt_level_(
              particles_->registerStateVariableOnly<int>("AdaptLevel"))
    {
        if (indication_method_.isFixedIndication())
        {
            particles_->addEvolvingVariable<int>(dv_adapt_level_);
        }
        particles_->addVariableToWrite<int>(dv_adapt_level_);
    };
    virtual ~AdaptLevelIndication() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : indication_(ex_policy, encloser.indication_method_),
              adapt_level_(
                  encloser.dv_adapt_level_->DelegatedData(ex_policy)){};

        void update(size_t index_i, Real dt = 0.0)
        {
            adapt_level_[index_i] = indication_(index_i);
        };

      protected:
        IndicationKernel indication_;
        int *adapt_level_;
    };

  protected:
    IndicationCriterion indication_method_;
    DiscreteVariable<int> *dv_adapt_level_;
};
} // namespace SPH
#endif // ADAPT_INDICATION_H
