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
 * @file general_gradient.h
 * @brief This the methods related on the corrections of SPH smoothing kernel.
 * @details The corrections aim to increase the numerical consistency
 * or accuracy for kernel approximations.
 * @author Xiangyu Hu
 */

#ifndef GENERAL_GRADIENT_H
#define GENERAL_GRADIENT_H

#include "base_general_dynamics.h"
#include "kernel_correction_ck.h"

namespace SPH
{
struct Divergence;

template <typename...>
struct DerivativeHelper;

template <int N>
struct DerivativeHelper<Divergence, Eigen::Matrix<Real, N, 1>>
{
    using type = Real;
};

template <int N>
struct DerivativeHelper<Divergence, Eigen::Matrix<Real, N, N>>
{
    using type = Eigen::Matrix<Real, N, 1>;
};

template <typename T>
using Div = typename DerivativeHelper<Divergence, T>::type;

struct Divergence
{
    Real operator()(const Vecd &difference, const Vecd &kernel_gradient) const
    {
        return difference.dot(kernel_gradient);
    }

    Vecd operator()(const Matd &difference, const Vecd &kernel_gradient) const
    {
        return difference * kernel_gradient;
    }
};

struct Gradient;

template <int N>
struct DerivativeHelper<Gradient, Eigen::Matrix<Real, N, 1>>
{
    using type = Eigen::Matrix<Real, N, N>;
};

template <>
struct DerivativeHelper<Gradient, Real>
{
    using type = Vecd;
};

template <typename T>
using Grad = typename DerivativeHelper<Gradient, T>::type;

struct Gradient
{
    Vecd operator()(const Real &difference, const Vecd &kernel_gradient) const
    {
        return difference * kernel_gradient;
    }

    Matd operator()(const Vecd &difference, const Vecd &kernel_gradient) const
    {
        return difference * kernel_gradient.transpose();
    }
};

template <typename...>
class Derivative;

template <typename OperatorType, typename DataType, typename KernelCorrectionType,
          template <typename...> class RelationType, typename... Parameters>
class Derivative<Base, OperatorType, DataType, KernelCorrectionType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseDynamicsType = Interaction<RelationType<Parameters...>>;

  public:
    explicit Derivative(RelationType<Parameters...> &relation, const std::string &variable_name);
    template <typename FirstArg>
    explicit Derivative(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : Derivative(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~Derivative() {};

  protected:
    using DerivativeType = typename DerivativeHelper<OperatorType, DataType>::type;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

    std::string variable_name_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<DerivativeType> *dv_derivative_;
    KernelCorrectionType kernel_correction_;

    std::string getDerivativeName(const Gradient &) const { return "Gradient"; }
    std::string getDerivativeName(const Divergence &) const { return "Divergence"; }
};

template <typename OperatorType, typename DataType, typename KernelCorrectionType, typename... Parameters>
class Derivative<OperatorType, DataType, KernelCorrectionType, Inner<Parameters...>>
    : public Derivative<Base, OperatorType, DataType, KernelCorrectionType, Inner<Parameters...>>
{

    using BaseDynamicsType = Derivative<Base, OperatorType, DataType, KernelCorrectionType, Inner<Parameters...>>;
    using DerivativeType = typename BaseDynamicsType::DerivativeType;
    using CorrectionKernel = typename BaseDynamicsType::CorrectionKernel;

  public:
    explicit Derivative(Inner<Parameters...> &relation, const std::string &variable_name);
    virtual ~Derivative() {};
    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        OperatorType operator_;
        DerivativeType zero_derivative_;
        Real *Vol_;
        DataType *variable_;
        DerivativeType *derivative_;
        CorrectionKernel correction_;
    };
};

template <typename OperatorType, typename DataType, typename KernelCorrectionType, typename... Parameters>
class Derivative<OperatorType, DataType, KernelCorrectionType, Contact<Parameters...>>
    : public Derivative<Base, OperatorType, DataType, KernelCorrectionType, Contact<Parameters...>>
{

    using BaseDynamicsType = Derivative<Base, OperatorType, DataType, KernelCorrectionType, Contact<Parameters...>>;
    using DerivativeType = typename BaseDynamicsType::DerivativeType;
    using CorrectionKernel = typename BaseDynamicsType::CorrectionKernel;

  public:
    explicit Derivative(Contact<Parameters...> &relation, const std::string &variable_name);
    virtual ~Derivative() {};
    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        OperatorType operator_;
        DerivativeType zero_derivative_;
        Real *contact_Vol_;
        DataType *variable_, *contact_variable_;
        DerivativeType *derivative_;
        CorrectionKernel correction_;
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> dv_contact_variable_;
};

template <typename...>
class LinearGradient;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class LinearGradient<RelationType<DataType, Parameters...>>
    : public Derivative<Gradient, DataType, LinearCorrectionCK, RelationType<Parameters...>>
{
  public:
    explicit LinearGradient(RelationType<Parameters...> &relation, const std::string &variable_name)
        : Derivative<Gradient, DataType, LinearCorrectionCK, RelationType<Parameters...>>(relation, variable_name) {};
    template <typename FirstArg>
    explicit LinearGradient(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : LinearGradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradient() {};
};

template <typename...>
class LinearDivergence;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class LinearDivergence<RelationType<DataType, Parameters...>>
    : public Derivative<Divergence, DataType, LinearCorrectionCK, RelationType<Parameters...>>
{
  public:
    explicit LinearDivergence(RelationType<Parameters...> &relation, const std::string &variable_name)
        : Derivative<Divergence, DataType, LinearCorrectionCK, RelationType<Parameters...>>(relation, variable_name) {};
    template <typename FirstArg>
    explicit LinearDivergence(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : LinearDivergence(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearDivergence() {};
};

template <typename...>
class PlainDivergence;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class PlainDivergence<RelationType<DataType, Parameters...>>
    : public Derivative<Divergence, DataType, NoKernelCorrectionCK, RelationType<Parameters...>>
{
  public:
    explicit PlainDivergence(RelationType<Parameters...> &relation, const std::string &variable_name)
        : Derivative<Divergence, DataType, NoKernelCorrectionCK, RelationType<Parameters...>>(relation, variable_name) {};
    template <typename FirstArg>
    explicit PlainDivergence(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : PlainDivergence(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~PlainDivergence() {};
};

template <typename...>
struct HessHelper;

template <>
struct HessHelper<Vecd>
{
    using type = VecMatGrad;
};

template <>
struct HessHelper<Real>
{
    using type = VecMatd;
};

template <typename T>
using Hess = typename HessHelper<T>::type;

template <typename...>
class Hessian;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class Hessian<Base, DataType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseDynamicsType = Interaction<RelationType<Parameters...>>;

  public:
    explicit Hessian(RelationType<Parameters...> &relation, const std::string &variable_name);
    template <typename FirstArg>
    explicit Hessian(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : Hessian(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        Real *Vol_;
        DataType *variable_;
        Grad<DataType> *gradient_;
        Matd *B_;
        MatTend *M_;
        Hess<DataType> *hessian_;
    };

  protected:
    std::string variable_name_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Grad<DataType>> *dv_gradient_;
    DiscreteVariable<Matd> *dv_B_;
    DiscreteVariable<MatTend> *dv_M_;
    DiscreteVariable<Hess<DataType>> *dv_hessian_;
};

template <typename DataType, typename... Parameters>
class Hessian<Inner<DataType, Parameters...>>
    : public Hessian<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Inner<Parameters...>>;

  public:
    template <typename... Args>
    explicit Hessian(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };
};

template <typename DataType, typename... Parameters>
class Hessian<Contact<DataType, Parameters...>>
    : public Hessian<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Contact<Parameters...>>;

  public:
    template <typename... Args>
    explicit Hessian(Args &&...args);
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> dv_contact_variable_;
};

template <typename... RelationTypes>
class SecondOrderGradient;

template <typename DataType, typename... Parameters>
class SecondOrderGradient<Inner<DataType, Parameters...>>
    : public Hessian<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Inner<Parameters...>>;

  public:
    template <typename... Args>
    explicit SecondOrderGradient(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~SecondOrderGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };
};

template <typename DataType, typename... Parameters>
class SecondOrderGradient<Contact<DataType, Parameters...>>
    : public Hessian<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Contact<Parameters...>>;

  public:
    template <typename... Args>
    explicit SecondOrderGradient(Args &&...args);
    virtual ~SecondOrderGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> dv_contact_variable_;
};
} // namespace SPH
#endif // GENERAL_GRADIENT_H
