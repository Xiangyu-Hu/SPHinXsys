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
 * or accuracy for kernel approximations. Note that the formulation is
 * not conservative (for divergence) but to achieve better accuracy.
 * @author Xiangyu Hu
 */

#ifndef GENERAL_GRADIENT_H
#define GENERAL_GRADIENT_H

#include "base_local_dynamics.h"

#include <string>
#include <tuple>
#include <type_traits>
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
struct HessHelper<Vec3d>
{
    using type = VecMat3d;
};

template <>
struct HessHelper<Real>
{
    using type = VecMatd;
};

template <typename T>
using Hess = typename HessHelper<T>::type;

template <typename... RelationTypes>
class Hessian;

template <typename... RelationTypes>
class Gradient;

template <typename... RelationTypes>
using LinearGradient = Gradient<RelationTypes...>;

template <typename T, typename = void>
struct TransportDataType
{
    using type = T;
};

template <typename T>
struct TransportDataType<T, std::void_t<typename T::TransportDataType>>
{
    using type = typename T::TransportDataType;
};

template <typename TransportType, typename DataType, typename = void>
struct InterParticleCoeffType
{
    struct type
    {
        DataType operator()(size_t, size_t) const { return DataType(1); }
    };
};

template <typename TransportType, typename DataType>
struct InterParticleCoeffType<TransportType, DataType, std::void_t<typename TransportType::InterParticleCoeff>>
{
    using type = typename TransportType::InterParticleCoeff;
};

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
    std::string variable_name_;
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Grad<DataType>> *dv_gradient_;
    DiscreteVariable<Matd> *dv_B_;
};

template <typename DataType, typename... Parameters>
class Gradient<Inner<DataType, Parameters...>>
    : public Gradient<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, Inner<Parameters...>>;

  public:
    explicit Gradient(Inner<Parameters...> &inner_relation, const std::string &variable_name)
        : BaseDynamicsType(inner_relation, variable_name) {};
    template <typename FirstArg>
    explicit Gradient(DynamicsArgs<Inner<Parameters...>, FirstArg> parameters)
        : Gradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~Gradient() {};

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
class Gradient<Contact<DataType, Parameters...>>
    : public Gradient<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, DataType, Contact<Parameters...>>;

  public:
    explicit Gradient(Contact<Parameters...> &contact_relation, const std::string &variable_name);
    template <typename FirstArg>
    explicit Gradient(DynamicsArgs<Contact<Parameters...>, FirstArg> parameters)
        : Gradient(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~Gradient() {};

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

template <typename TransportType, template <typename...> class RelationType, typename... Parameters>
class Hessian<Base, TransportType, RelationType<Parameters...>>
    : public Gradient<Base, typename TransportDataType<TransportType>::type, RelationType<Parameters...>>
{
    using DataType = typename TransportDataType<TransportType>::type;
    using BaseDynamicsType = Gradient<Base, DataType, RelationType<Parameters...>>;

  public:
    template <typename... Args>
    explicit Hessian(Args &&...args);
    virtual ~Hessian() {};

  protected:
    TransportType *transport_;
    DiscreteVariable<MatTend> *dv_M_;
    DiscreteVariable<Hess<DataType>> *dv_hessian_;
};

template <typename TransportType, typename... Parameters>
class Hessian<Inner<TransportType, Parameters...>>
    : public Hessian<Base, TransportType, Inner<Parameters...>>
{
    using DataType = typename TransportDataType<TransportType>::type;
    using BaseDynamicsType = Hessian<Base, TransportType, Inner<Parameters...>>;
    using InterParticleCoeff = typename InterParticleCoeffType<TransportType, DataType>::type;

  public:
    template <typename... Args>
    explicit Hessian(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InterParticleCoeff inter_particle_coeff_;
        MatTend *M_;
        Hess<DataType> *hessian_;
    };
};

template <typename TransportType, typename... OtherTransportType,
          template <typename...> class InterfaceType, typename... Parameters>
class Hessian<Contact<InterfaceType<TransportType, OtherTransportType...>, Parameters...>>
    : public Hessian<Base, TransportType, Contact<Parameters...>>
{
    using DataType = typename TransportDataType<TransportType>::type;
    using BaseDynamicsType = Hessian<Base, TransportType, Contact<Parameters...>>;
    using InterfaceModel = InterfaceType<TransportType, OtherTransportType...>;
    using InterfaceCoeff = typename InterfaceModel::Coefficient;

  public:
    template <typename FirstArg>
    explicit Hessian(DynamicsArgs<Contact<Parameters...>, FirstArg,
                                  StdVec<InterfaceType<TransportType, OtherTransportType...> *>> parameters);
    Hessian(Contact<Parameters...> &contact_relation, std::string &variable_name,
            StdVec<InterfaceType<TransportType, OtherTransportType...> *> interface_models);
    virtual ~Hessian() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InterfaceCoeff interface_coeff_;
        MatTend *M_;
        Hess<DataType> *hessian_;
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    StdVec<InterfaceModel *> interface_models_;
    StdVec<DiscreteVariable<DataType> *> dv_contact_variable_;
};

template <typename... RelationTypes>
class SecondOrderGradient;

template <typename DataType, typename... Parameters>
class SecondOrderGradient<Inner<DataType, Parameters...>>
    : public Hessian<Inner<DataType, Parameters...>>
{
    using BaseDynamicsType = Hessian<Inner<DataType, Parameters...>>;

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

template <typename...>
struct CurlHelper;

template <>
struct CurlHelper<Vec2d>
{
    using type = Real;
};

template <>
struct CurlHelper<Vec3d>
{
    using type = Vec3d;
};

using Curl = typename CurlHelper<Vecd>::type;

template <typename...>
class LinearCurl;

template <typename... Parameters>
class LinearCurl<Inner<WithUpdate, Parameters...>>
    : public Gradient<Inner<Vecd, Parameters...>>
{
    using BaseDynamicsType = Gradient<Inner<Vecd, Parameters...>>;

  public:
    LinearCurl(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~LinearCurl() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Grad<Vecd> *gradient_;
        Curl *curl_;
    };

  protected:
    DiscreteVariable<Curl> *dv_curl_;
};

template <typename... Parameters>
class LinearCurl<Contact<Parameters...>>
    : public Gradient<Contact<Vecd, Parameters...>>
{
  public:
    LinearCurl(Contact<Parameters...> &contact_relation, const std::string &variable_name)
        : Gradient<Contact<Vecd, Parameters...>>(contact_relation, variable_name) {};
    virtual ~LinearCurl() {};
};

template <typename...>
class DoubleCurl;

template <typename TransportType, typename... Parameters>
class DoubleCurl<Inner<TransportType, Parameters...>>
    : public Hessian<Inner<TransportType, Parameters...>>
{
    using DataType = typename TransportDataType<TransportType>::type;
    using BaseDynamicsType = Hessian<Inner<TransportType, Parameters...>>;

  public:
    DoubleCurl(Inner<Parameters...> &inner_relation, const std::string &variable_name);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Hess<DataType> *hessian_;
        Vecd *double_curl_;

        Vec2d computeDoubleCurl(Hess<Vec2d> &hessian);
        Vec3d computeDoubleCurl(Hess<Vec3d> &hessian);
    };

  protected:
    DiscreteVariable<Vecd> *dv_double_curl_;
};

} // namespace SPH
#endif // GENERAL_GRADIENT_H
