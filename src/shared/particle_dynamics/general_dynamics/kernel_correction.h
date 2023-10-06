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
class KernelCorrectionMatrix;

template <class DataDelegationType>
class KernelCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit KernelCorrectionMatrix(BaseRelationType &base_relation, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrix(){};

  protected:
    Real alpha_;
    StdLargeVec<Matd> &B_;
};

template <>
class KernelCorrectionMatrix<Inner>
    : public KernelCorrectionMatrix<GeneralDataDelegateInner>
{
  public:
    explicit KernelCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : KernelCorrectionMatrix<GeneralDataDelegateInner>(inner_relation){};
    explicit KernelCorrectionMatrix(ConstructorArgs<BaseInnerRelation, Real> parameters)
        : KernelCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~KernelCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using KernelCorrectionMatrixInner = KernelCorrectionMatrix<Inner>;

template <>
class KernelCorrectionMatrix<Contact>
    : public KernelCorrectionMatrix<GeneralDataDelegateContact>
{
  public:
    explicit KernelCorrectionMatrix(BaseContactRelation &contact_relation, Real alpha = Real(0));
    explicit KernelCorrectionMatrix(ConstructorArgs<BaseContactRelation, Real> parameters)
        : KernelCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~KernelCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
};

template <class InnerInteractionType, class ContactInteractionType>
class BaseKernelCorrectionMatrixComplex
    : public ComplexInteraction<KernelCorrectionMatrix<InnerInteractionType, ContactInteractionType>>
{
  public:
    explicit BaseKernelCorrectionMatrixComplex(ComplexRelation &complex_relation)
        : ComplexInteraction<KernelCorrectionMatrix<InnerInteractionType, ContactInteractionType>>(
              complex_relation.getInnerRelation(), complex_relation.getContactRelation()){};
};
using KernelCorrectionMatrixComplex = BaseKernelCorrectionMatrixComplex<Inner, Contact>;

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
class KernelGradientCorrection<Inner>
    : public KernelGradientCorrection<GeneralDataDelegateInner>
{
    PairAverageInner<Matd> average_correction_matrix_;

  public:
    explicit KernelGradientCorrection(BaseInnerRelation &inner_relation);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class KernelGradientCorrection<Contact>
    : public KernelGradientCorrection<GeneralDataDelegateContact>
{
    StdVec<PairAverageContact<Matd>> contact_average_correction_matrix_;

  public:
    KernelGradientCorrection(BaseContactRelation &contact_relation);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class InnerInteractionType, class ContactInteractionType>
class BaseKernelGradientCorrectionComplex
    : public ComplexInteraction<KernelGradientCorrection<InnerInteractionType, ContactInteractionType>>
{
  public:
    explicit BaseKernelGradientCorrectionComplex(ComplexRelation &complex_relation)
        : ComplexInteraction<KernelGradientCorrection<InnerInteractionType, ContactInteractionType>>(
              complex_relation.getInnerRelation(), complex_relation.getContactRelation()){};
};
using KernelGradientCorrectionComplex = BaseKernelGradientCorrectionComplex<Inner, Contact>;
} // namespace SPH
#endif // KERNEL_CORRECTION_H
