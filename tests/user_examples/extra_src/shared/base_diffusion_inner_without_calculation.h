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
 * @file    diffusion_dynamics.h
 * @brief   These are particle dynamics applicable for all type of particles.
 * @author  Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_DIFFUSION_INNER_WITHOUT_CALCULATION_H
#define BASE_DIFFUSION_INNER_WITHOUT_CALCULATION_H

#include "diffusion_dynamics.h"

namespace SPH
{
template <class ParticlesType, class KernelGradientType = KernelGradientInner>
class DiffusionRelaxationInnerWithoutCalculation
    : public DiffusionRelaxationInner<ParticlesType, KernelGradientType>
{
  public:
    typedef BaseInnerRelation BodyRelationType;
    explicit DiffusionRelaxationInnerWithoutCalculation(BaseInnerRelation &inner_relation) : DiffusionRelaxationInner<ParticlesType, KernelGradientType>(inner_relation){};
    virtual ~DiffusionRelaxationInnerWithoutCalculation(){};
    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        /* initialize Diffusion ChangeRate*/
        for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
        {
            (*this->diffusion_dt_[m])[index_i] = 0;
        }
    };
};

} // namespace SPH
#endif // BASE_DIFFUSION_INNER_WITHOUT_CALCULATION_H