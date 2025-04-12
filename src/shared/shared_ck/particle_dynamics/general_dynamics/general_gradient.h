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
struct Grad;

template <int N>
struct Grad<Eigen::Matrix<Real, N, 1>>
{
    typedef Eigen::Matrix<Real, N, N> GradType;
};

template <>
struct Grad<Real>
{
    typedef Vecd GradType;
};

template <typename...>
struct Hessian;

template <>
struct Hessian<Vecd>
{
    typedef VecMatGrad HessianType;
};

template <>
struct Hessian<Real>
{
    typedef VecMatd HessianType;
};

template <typename... RelationTypes>
class Gradient;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class Gradient<Base, DataType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseDynamicsType = Interaction<RelationType<Parameters...>>;
    using GradType = typename Grad<DataType>::GradType;

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
        GradType *gradient_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<GradType> *dv_gradient_;
};

template <typename... RelationTypes>
class HessianMatrix;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class HessianMatrix<Base, DataType, RelationType<Parameters...>>
    : public Gradient<Base, DataType, RelationType<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, RelationType<Parameters...>>;
    using HessianType = typename Hessian<DataType>::HessianType;

  public:
    template <class DynamicsIdentifier>
    explicit HessianMatrix(DynamicsIdentifier &identifier, std::string &variable_name);
    virtual ~HessianMatrix() {};

    class InteractKernel : public BaseDynamicsTypeInteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        Matd *B_;
        MatTend *M_;
        HessianType *hessian_;

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
    DiscreteVariable<Matd> *dv_B_;
    DiscreteVariable<MatTend> *dv_M_;
    DiscreteVariable<HessianType> *dv_hessian_;
};

template <typename DataType, typename... Parameters>
class HessianMatrix<Inner<DataType, Parameters...>>
    : public HessianMatrix<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianMatrix<Base, DataType, Inner<Parameters...>>;

  public:
    explicit HessianMatrix(Relation<Inner<Parameters...>> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
    virtual ~HessianMatrix() {};

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
class HessianMatrix<Contact<DataType, Parameters...>>
    : public HessianMatrix<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianMatrix<Base, DataType, Contact<Parameters...>>;

  public:
    explicit HessianMatrix(Relation<Contact<Parameters...>> &contact_relation)
        : BaseDynamicsType(contact_relation) {};
    virtual ~HessianMatrix() {};

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
    : public HessianMatrix<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianMatrix<Base, DataType, Inner<Parameters...>>;

  public:
    explicit SecondOrderGradient(Relation<Inner<Parameters...>> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
    virtual ~HessianMatrix() {};

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
    : public HessianMatrix<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianMatrix<Base, DataType, Contact<Parameters...>>;

  public:
    explicit HessianMatrix(Relation<Contact<Parameters...>> &contact_relation)
        : BaseDynamicsType(contact_relation) {};
    virtual ~HessianMatrix() {};

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
