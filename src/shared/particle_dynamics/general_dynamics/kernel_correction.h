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
 * @file 	kernel_correction.h
 * @brief 	This the methods related on the corrections of SPH smoothing kernel.
 * @author	Yaru Ren and Xiangyu Hu
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
    explicit KernelCorrectionMatrix(ConstructorArgs<BaseContactRelation, Real> parameters)
        : KernelCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~KernelCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

template <>
class KernelCorrectionMatrix<Contact>
    : public KernelCorrectionMatrix<GeneralDataDelegateContact>
{
  public:
    explicit KernelCorrectionMatrix(BaseContactRelation &contact_relation, Real alpha = Real(0));
    virtual ~KernelCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
};

template <class InnerInteractionType, class ContactInteractionType>
class BaseKernelCorrectionMatrixComplex
    : public ComplexInteraction<eKernelGradientCorrection<InnerInteractionType, ContactInteractionType>>
{
  public:
    explicit BaseDensitySummationComplex(ComplexRelation &complex_relation)
        : ComplexInteraction<eKernelGradientCorrection<InnerInteractionType, ContactInteractionType>>(
              complex_relation.getInnerRelation(), complex_relation.getContactRelation()){};
};
using KernelGradientCorrectionComplex = BaseKernelGradientCorrectionComplex<Inner, Contact>;

template <typename... InteractionTypes>
class KernelGradientCorrection;

template <class DataDelegationType>
class KernelGradientCorrection<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class KernelCorrectionMatrixType>
    explicit KernelGradientCorrection(KernelCorrectionMatrixType &kernel_correction);
    virtual ~KernelGradientCorrection(){};

  protected:
    template <class PairAverageType>
    void correctKernelGradient(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i);
};

class KernelGradientCorrection<Inner>
    : public KernelGradientCorrection<GeneralDataDelegateInner>
{
    PairAverageInner<Matd> average_correction_matrix_;

  public:
    explicit KernelGradientCorrection(KernelCorrectionMatrix<Inner> &kernel_correction_inner);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

class KernelGradientCorrection<Contact>
    : public KernelGradientCorrection<GeneralDataDelegateContact>
{
    StdVec<PairAverageContact<Matd>> contact_average_correction_matrix_;

  public:
    KernelGradientCorrection(KernelCorrectionMatrix<Contact> &kernel_correction_contact);
    virtual ~KernelGradientCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class InnerInteractionType, class ContactInteractionType>
class BaseKernelGradientCorrectionComplex
    : public ComplexInteraction<eKernelGradientCorrection<InnerInteractionType, ContactInteractionType>>
{
  public:
    explicit BaseKernelGradientCorrectionComplex(KernelGradientCorrectionComplex &kernel_correction_complex)
        : ComplexInteraction<eKernelGradientCorrection<InnerInteractionType, ContactInteractionType>>(
              complex_relation.getInnerRelation(), complex_relation.getContactRelation()){};
};
using KernelGradientCorrectionComplex = BaseKernelGradientCorrectionComplex<Inner, Contact>;

/**
 * @class ParticleSmoothing
 * @brief computing smoothed variable field by averaging with neighbors
 */
template <typename VariableType>
class ParticleSmoothing : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    explicit ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name)
        : LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
          W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)),
          smoothed_(*particles_->template getVariableByName<VariableType>(variable_name))
    {
        particles_->registerVariable(temp_, variable_name + "_temp");
    }

    virtual ~ParticleSmoothing(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real weight = W0_;
        VariableType summation = W0_ * smoothed_[index_i];
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            summation += inner_neighborhood.W_ij_[n] * smoothed_[index_j];
            weight += inner_neighborhood.W_ij_[n];
        }
        temp_[index_i] = summation / (weight + TinyReal);
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        smoothed_[index_i] = temp_[index_i];
    };

  protected:
    const Real W0_;
    StdLargeVec<VariableType> &smoothed_, temp_;
};
} // namespace SPH
#endif // KERNEL_CORRECTION_H
