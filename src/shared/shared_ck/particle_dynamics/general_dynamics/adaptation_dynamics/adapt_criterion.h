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
 * @file    adapt_criterion.h
 * @brief   to define the criterion of indication for adaptation
 * @author	Xiangyu Hu
 */

#ifndef ADAPT_CRITERION_H
#define ADAPT_CRITERION_H

#include "base_general_dynamics.h"

namespace SPH
{
template <typename...>
class Refinement;

template <>
class Refinement<Base>
{
  public:
    Refinement(const std::string &name, bool is_fixed_indication)
        : name_(name), is_fixed_indication_(is_fixed_indication) {};
    virtual ~Refinement() {};
    std::string getName() { return name_; };
    bool isFixedIndication() { return is_fixed_indication_; };

  protected:
    std::string name_;
    bool is_fixed_indication_;
};

template <typename T>
class Refinement<Continuous, T> : public Refinement<Base>
{
  public:
    Refinement(SPHAdaptation *sph_adaptation, BaseParticles *particles)
        : Refinement<Base>("Refinement", T::is_fixed),
          refinement_level_(sph_adaptation->LocalRefinementLevel()),
          h_ref_(sph_adaptation->ReferenceSmoothingLength()),
          dv_h_(particles->getVariableByName<Real>("SmoothingLength")) {};
    virtual ~Refinement() {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : refinement_level_(encloser.refinement_level_), h_ref_(encloser.h_ref_),
              h_(encloser.dv_h_->DelegatedData(ex_policy)){};
        int operator()(size_t index_i, Real dt = 0.0)
        {
            Real h_level = h_ref_;
            Real h_current = h_[index_i];
            for (int j = 0; j < refinement_level_; ++j)
            {
                h_level *= 0.5;
                if (h_current > h_level)
                {
                    return j;
                }
            }
            return refinement_level_;
        };

      protected:
        int refinement_level_;
        Real h_ref_;
        Real *h_;
    };

  protected:
    int refinement_level_;
    Real h_ref_;
    DiscreteVariable<Real> *dv_h_;
};
} // namespace SPH
#endif // ADAPT_CRITERION_H
