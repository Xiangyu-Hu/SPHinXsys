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
class AdaptIndication : public LocalDynamics
{
    using IndicationKernel = typename IndicationCriterion::IndicationKernel;

  public:
    explicit AdaptIndication(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          indication_method_(sph_adaptation_, particles_),
          dv_adapt_indicator_(
              particles_->registerStateVariableOnly(
                  indication_method_.getName() + "Indicator"))
    {
        if (indication_method_.isFixedIndication())
        {
            particles_->addEvolvingVariable(dv_adapt_indicator_);
        }
        particles_.addVariableToWrite(dv_adapt_indicator_);
    };
    virtual ~AdaptIndication() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : indication_(ex_policy, encloser.indication_method_),
              adaptation_indicator_(
                  encloser.dv_adapt_indicator_->DelegatedData(ex_policy)){};

        void update(size_t index_i, Real dt = 0.0)
        {
            adaptation_indicator_[index_i] = indication_(index_i);
        };

      protected:
        IndicationKernel indication_;
        int *adaptation_indicator_;
    };

  protected:
    IndicationCriterion indication_method_;
    DiscreteVariable<int> *dv_adapt_indicator_;
};
} // namespace SPH
#endif // ADAPT_INDICATION_H
