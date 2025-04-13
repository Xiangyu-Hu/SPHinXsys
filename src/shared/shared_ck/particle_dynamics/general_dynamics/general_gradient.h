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
#include "scalar_numerics.h"

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
    template <class DynamicsIdentifier>
    explicit Gradient(DynamicsIdentifier &identifier, std::string &variable_name);
    virtual ~Gradient() {};

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
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Grad<DataType>> *dv_gradient_;
    DiscreteVariable<Matd> *dv_B_;
};

template <typename... RelationTypes>
class LinearGradient;

template <typename DataType, typename... Parameters>
class LinearGradient<Inner<DataType, Parameters...>>
    : public Gradient<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, Inner<Parameters...>>;

  public:
    explicit LinearGradient(Relation<Inner<Parameters...>> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
    virtual ~LinearGradient() {};

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
    explicit LinearGradient(Relation<Contact<Parameters...>> &contact_relation)
        : BaseDynamicsType(contact_relation) {};
    virtual ~LinearGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
              contact_Vol_(encloser.dv_contact_Vol_->DataDelegate(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

template <typename... RelationTypes>
class Hessian;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class Hessian<Base, DataType, RelationType<Parameters...>>
    : public Gradient<Base, DataType, RelationType<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, RelationType<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit Hessian(DynamicsIdentifier &identifier, std::string &variable_name);
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        MatTend *M_;
        Hess<DataType> *hessian_;

        Eigen::Matrix<Real, 1, 1> transferToMatrix(Real value)
        {
            return Eigen::Matrix<Real, 1, 1>::Identity() * value;
        };

        Vecd transferToMatrix(const Vecd &value)
        {
            return value;
        };
    };

  protected:
    DiscreteVariable<MatTend> *dv_M_;
    DiscreteVariable<Hess<DataType>> *dv_hessian_;
};

template <typename DataType, typename... Parameters>
class Hessian<Inner<DataType, Parameters...>>
    : public Hessian<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Inner<Parameters...>>;

  public:
    explicit Hessian(Relation<Inner<Parameters...>> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
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
    explicit Hessian(Relation<Contact<Parameters...>> &contact_relation)
        : BaseDynamicsType(contact_relation) {};
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
              contact_Vol_(encloser.dv_contact_Vol_->DataDelegate(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

template <typename... RelationTypes>
class SecondOrderGradient;

template <typename DataType, typename... Parameters>
class SecondOrderGradient<Inner<DataType, Parameters...>>
    : public Hessian<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Hessian<Base, DataType, Inner<Parameters...>>;

  public:
    explicit SecondOrderGradient(Relation<Inner<Parameters...>> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
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
    explicit SecondOrderGradient(Relation<Contact<Parameters...>> &contact_relation)
        : BaseDynamicsType(contact_relation) {};
    virtual ~SecondOrderGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
              contact_Vol_(encloser.dv_contact_Vol_->DataDelegate(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};
} // namespace SPH
#endif // GENERAL_GRADIENT_H
