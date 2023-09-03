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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	general_interaction.h
 * @brief 	This is the interaction dynamics applicable for all type bodies
 * @author	Yaru Ren and Xiangyu Hu
 */

#ifndef GENERAL_INTERACTION_H
#define GENERAL_INTERACTION_H

#include "general_dynamics.h"

namespace SPH
{
class KernelCorrectionMatrixInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    KernelCorrectionMatrixInner(BaseInnerRelation &inner_relation, int beta = 0, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrixInner(){};

  protected:
    int beta_;
    Real alpha_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class KernelCorrectionMatrixComplex : public KernelCorrectionMatrixInner, public GeneralDataDelegateContactOnly
{
  public:
    KernelCorrectionMatrixComplex(ComplexRelation &complex_relation, int beta = 0, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrixComplex(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;

    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class KernelGradientCorrectionInner
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient
 */
class KernelGradientCorrectionInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    KernelGradientCorrectionInner(KernelCorrectionMatrixInner &kernel_correction_inner);
    virtual ~KernelGradientCorrectionInner(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &B_;
    void correctKernelGradient(Neighborhood &neighborhood, size_t index_i);
};

/**
 * @class KernelGradientCorrectionComplex
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient in complex topology
 */
class KernelGradientCorrectionComplex : public KernelGradientCorrectionInner, public GeneralDataDelegateContactOnly
{
  public:
    KernelGradientCorrectionComplex(KernelCorrectionMatrixComplex &kernel_correction_complex);
    virtual ~KernelGradientCorrectionComplex(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_INTERACTION_H