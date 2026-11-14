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

#include "base_local_dynamics.h"

#include <string>
#include <tuple>
#include <utility>

namespace SPH
{
template <typename...>
struct GradHelper;

template <int N>
struct GradHelper<Eigen::Matrix<Real, N, 1>>
{
    using type = Eigen::Matrix<Real, N, N>;
};

template <>
struct GradHelper<Real>
{
    using type = Vecd;
};

template <typename T>
using Grad = typename GradHelper<T>::type;

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

template <typename... RelationTypes>
class Gradient;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class Gradient<Base, DataType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseDynamicsType = Interaction<RelationType<Parameters...>>;

  public:
    explicit Gradient(RelationType<Parameters...> &relation, const std::string &variable_name);
    template <typename FirstArg>
    explicit Gradient(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : Gradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~Gradient(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        Real *Vol_;
        DataType *variable_;
        Grad<DataType> *gradient_;
    };

  protected:
    std::string variable_name_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Grad<DataType>> *dv_gradient_;
};

template <typename... RelationTypes>
class LinearGradient;

template <typename DataType, typename... Parameters>
class LinearGradient<Inner<DataType, Parameters...>>
    : public Gradient<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit LinearGradient(DynamicsIdentifier &identifier, const std::string &variable_name)
        : BaseDynamicsType(identifier, variable_name){};
    template <typename FirstArg>
    explicit LinearGradient(DynamicsArgs<Inner<Parameters...>, FirstArg> parameters)
        : LinearGradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradient(){};

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
class LinearGradient<Contact<DataType, Parameters...>>
    : public Gradient<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit LinearGradient(DynamicsIdentifier &identifier, const std::string &variable_name);
    template <typename FirstArg>
    explicit LinearGradient(DynamicsArgs<Contact<Parameters...>, FirstArg> parameters)
        : LinearGradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradient(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
};

template <typename... RelationTypes>
class Hessian;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class Hessian<Base, DataType, RelationType<Parameters...>>
    : public Gradient<Base, DataType, RelationType<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, RelationType<Parameters...>>;

  public:
    template <typename... Args>
    explicit Hessian(Args &&...args);
    virtual ~Hessian(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        Hess<DataType> *hessian_;
    };

  protected:
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
    virtual ~Hessian(){};

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
    virtual ~Hessian(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
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
    virtual ~SecondOrderGradient(){};

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
    virtual ~SecondOrderGradient(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
};

//=================================================================================================//
template <typename... RelationTypes>
class GradientApplyCorrection;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class GradientApplyCorrection<Base, DataType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseDynamicsType = Interaction<RelationType<Parameters...>>;

  public:
    explicit GradientApplyCorrection(RelationType<Parameters...> &relation, const std::string &variable_name);
    template <typename FirstArg>
    explicit GradientApplyCorrection(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : GradientApplyCorrection(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~GradientApplyCorrection(){};

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
    };

  protected:
    std::string variable_name_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Grad<DataType>> *dv_gradient_;
    DiscreteVariable<Matd> *dv_B_;
};

template <typename... RelationTypes>
class LinearGradientApplyCorrection;

template <typename DataType, typename... Parameters>
class LinearGradientApplyCorrection<Inner<DataType, Parameters...>>
    : public GradientApplyCorrection<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = GradientApplyCorrection<Base, DataType, Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit LinearGradientApplyCorrection(DynamicsIdentifier &identifier, const std::string &variable_name)
        : BaseDynamicsType(identifier, variable_name){};
    template <typename FirstArg>
    explicit LinearGradientApplyCorrection(DynamicsArgs<Inner<Parameters...>, FirstArg> parameters)
        : LinearGradientApplyCorrection(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradientApplyCorrection(){};

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
class LinearGradientApplyCorrection<Contact<DataType, Parameters...>>
    : public GradientApplyCorrection<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = GradientApplyCorrection<Base, DataType, Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit LinearGradientApplyCorrection(DynamicsIdentifier &identifier, const std::string &variable_name);
    template <typename FirstArg>
    explicit LinearGradientApplyCorrection(DynamicsArgs<Contact<Parameters...>, FirstArg> parameters)
        : LinearGradientApplyCorrection(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearGradientApplyCorrection(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
};

template <typename... RelationTypes>
class HessianApplyCorrection;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class HessianApplyCorrection<Base, DataType, RelationType<Parameters...>>
    : public GradientApplyCorrection<Base, DataType, RelationType<Parameters...>>
{
    using BaseDynamicsType = GradientApplyCorrection<Base, DataType, RelationType<Parameters...>>;

  public:
    template <typename... Args>
    explicit HessianApplyCorrection(Args &&...args);
    virtual ~HessianApplyCorrection(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        MatTend *M_;
        Hess<DataType> *hessian_;
    };

  protected:
    DiscreteVariable<MatTend> *dv_M_;
    DiscreteVariable<Hess<DataType>> *dv_hessian_;
};

template <typename DataType, typename... Parameters>
class HessianApplyCorrection<Inner<DataType, Parameters...>>
    : public HessianApplyCorrection<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianApplyCorrection<Base, DataType, Inner<Parameters...>>;

  public:
    template <typename... Args>
    explicit HessianApplyCorrection(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~HessianApplyCorrection(){};

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
class HessianApplyCorrection<Contact<DataType, Parameters...>>
    : public HessianApplyCorrection<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianApplyCorrection<Base, DataType, Contact<Parameters...>>;

  public:
    template <typename... Args>
    explicit HessianApplyCorrection(Args &&...args);
    virtual ~HessianApplyCorrection(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
};

template <typename... RelationTypes>
class SecondOrderGradientApplyCorrection;

template <typename DataType, typename... Parameters>
class SecondOrderGradientApplyCorrection<Inner<DataType, Parameters...>>
    : public HessianApplyCorrection<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianApplyCorrection<Base, DataType, Inner<Parameters...>>;

  public:
    template <typename... Args>
    explicit SecondOrderGradientApplyCorrection(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~SecondOrderGradientApplyCorrection(){};

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
class SecondOrderGradientApplyCorrection<Contact<DataType, Parameters...>>
    : public HessianApplyCorrection<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianApplyCorrection<Base, DataType, Contact<Parameters...>>;

  public:
    template <typename... Args>
    explicit SecondOrderGradientApplyCorrection(Args &&...args);
    virtual ~SecondOrderGradientApplyCorrection(){};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_contact_variable_;
};
} // namespace SPH
#endif // GENERAL_GRADIENT_H
