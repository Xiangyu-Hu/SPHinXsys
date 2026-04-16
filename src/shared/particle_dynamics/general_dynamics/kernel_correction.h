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
#define KERNEL_CORRECTION_H //防止头文件被多次包含

#include "base_general_dynamics.h"

namespace SPH
{
template <typename... InteractionTypes>//可变参数模板 ... 表示可以接受任意多个类型参数，比如 <Inner<>, Contact<>> 或 <Inner<>> 等
class LinearGradientCorrectionMatrix;//主模板声明（只有声明，没有定义）

// WKGC1 will be used for calculate the KGC matrix.
// The difference between WKGC1 and WKGC2 can refer to https://doi.org/10.1016/j.cma.2023.116460

template <class DataDelegationType>//部分特化：处理一种交互类型
class LinearGradientCorrectionMatrix<DataDelegationType>//当用户只给一个模板参数时，使用这个版本
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit LinearGradientCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~LinearGradientCorrectionMatrix() {};

  protected:
    Real *Vol_;
    Matd *B_;
};

template <>//表示这是一个完全特化，模板参数列表为空
class LinearGradientCorrectionMatrix<Inner<>>//只处理内部粒子
    : public LinearGradientCorrectionMatrix<DataDelegateInner>
{
    Real alpha_;

  public:
    explicit LinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : LinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation), alpha_(alpha) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit LinearGradientCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : LinearGradientCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradientCorrectionMatrix() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using LinearGradientCorrectionMatrixInner = LinearGradientCorrectionMatrix<Inner<>>;//定义别名，这样以后写 LinearGradientCorrectionMatrixInner 就等价于这个特化类

template <>
class LinearGradientCorrectionMatrix<Contact<>>
    : public LinearGradientCorrectionMatrix<DataDelegateContact>
{
  public:
    explicit LinearGradientCorrectionMatrix(BaseContactRelation &contact_relation);
    virtual ~LinearGradientCorrectionMatrix() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Real *> contact_mass_;
};

using LinearGradientCorrectionMatrixComplex = ComplexInteraction<LinearGradientCorrectionMatrix<Inner<>, Contact<>>>;//ComplexInteraction 是一个辅助模板，把内部和接触的修正组合起来，一次调用就能同时处理内部和边界邻居。

template <typename... InteractionTypes>
class KernelGradientCorrection;

template <class DataDelegationType>
class KernelGradientCorrection<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit KernelGradientCorrection(BaseRelationType &base_relation);
    virtual ~KernelGradientCorrection() {};

  protected:
    template <class PairAverageType>
    void correctKernelGradient(PairAverageType &average_correction_matrix, Neighborhood &neighborhood, size_t index_i);
};

template <>
class KernelGradientCorrection<Inner<>>
    : public KernelGradientCorrection<DataDelegateInner>
{
    PairAverageVariable<Matd> average_correction_matrix_;//它是每个粒子对的平均修正矩阵 (B_i + B_j)/2。

  public:
    explicit KernelGradientCorrection(BaseInnerRelation &inner_relation);
    virtual ~KernelGradientCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0);
};
using KernelGradientCorrectionInner = KernelGradientCorrection<Inner<>>;

template <>
class KernelGradientCorrection<Contact<>>
    : public KernelGradientCorrection<DataDelegateContact>
{
    StdVec<PairAverageVariable<Matd>> contact_average_correction_matrix_;

  public:
    KernelGradientCorrection(BaseContactRelation &contact_relation);
    virtual ~KernelGradientCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0);
};

using KernelGradientCorrectionComplex = ComplexInteraction<KernelGradientCorrection<Inner<>, Contact<>>>;
} // namespace SPH
#endif // KERNEL_CORRECTION_H
