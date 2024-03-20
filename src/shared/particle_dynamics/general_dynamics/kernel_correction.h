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
 * @file kernel_correction.h
 * @brief This the methods related on the corrections of SPH smoothing kernel.
 * @details The corrections aim to increase the numerical consistency
 * or accuracy for kernel approximations.
 * @author Yaru Ren and Xiangyu Hu
 */

#ifndef KERNEL_CORRECTION_H
#define KERNEL_CORRECTION_H

#include "base_general_dynamics.h"

namespace SPH
{
template <typename... InteractionTypes>
class LinearCorrectionMatrix;

template <class DataDelegationType>
class LinearCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit LinearCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~LinearCorrectionMatrix(){};

  protected:
    StdLargeVec<Matd> &B_;
};

template <>
class LinearCorrectionMatrix<Inner<>>
    : public LinearCorrectionMatrix<GeneralDataDelegateInner>
{
    Real alpha_;

  public:
    explicit LinearCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : LinearCorrectionMatrix<GeneralDataDelegateInner>(inner_relation), alpha_(alpha){};
    template <typename BodyRelationType, typename FirstArg>
    explicit LinearCorrectionMatrix(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : LinearCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~LinearCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using LinearCorrectionMatrixInner = LinearCorrectionMatrix<Inner<>>;

template <>
class LinearCorrectionMatrix<Contact<>>
    : public LinearCorrectionMatrix<GeneralDataDelegateContact>
{
  public:
    explicit LinearCorrectionMatrix(BaseContactRelation &contact_relation);
    virtual ~LinearCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
};

using LinearCorrectionMatrixComplex = ComplexInteraction<LinearCorrectionMatrix<Inner<>, Contact<>>>;

template <typename... InteractionTypes>
class KernelGradientCorrection;

template <class DataDelegationType>
class KernelGradientCorrection<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit KernelGradientCorrection(BaseRelationType &base_relation);
    virtual ~KernelGradientCorrection(){};

  protected:
    template <class PairAverageType>
    void correctKernelGradient(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i);
};

template <>
class KernelGradientCorrection<Inner<>>
    : public KernelGradientCorrection<GeneralDataDelegateInner>
{
    PairAverageVariable<Matd> average_correction_matrix_;

  public:
    explicit KernelGradientCorrection(BaseInnerRelation &inner_relation);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using KernelGradientCorrectionInner = KernelGradientCorrection<Inner<>>;

template <>
class KernelGradientCorrection<Contact<>>
    : public KernelGradientCorrection<GeneralDataDelegateContact>
{
    StdVec<PairAverageVariable<Matd>> contact_average_correction_matrix_;

  public:
    KernelGradientCorrection(BaseContactRelation &contact_relation);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

using KernelGradientCorrectionComplex = ComplexInteraction<KernelGradientCorrection<Inner<>, Contact<>>>;
} // namespace SPH
#endif // KERNEL_CORRECTION_H
