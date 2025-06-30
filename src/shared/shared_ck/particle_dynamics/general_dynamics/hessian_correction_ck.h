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
 * @file hessian_correction_ck.h
 * @brief This the methods related on the corrections of SPH smoothing kernel.
 * @details The corrections aim to increase the numerical consistency
 * or accuracy for kernel approximations.
 * @author Xiangyu Hu
 */

#ifndef HESSIAN_CORRECTION_CK_H
#define HESSIAN_CORRECTION_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
template <typename... RelationTypes>
class HessianCorrectionMatrix;

template <template <typename...> class RelationType, typename... Parameters>
class HessianCorrectionMatrix<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit HessianCorrectionMatrix(DynamicsIdentifier &identifier);
    virtual ~HessianCorrectionMatrix() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       HessianCorrectionMatrix<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

      protected:
        Real *Vol_;
        Matd *B_;
        VecMatGrad *displacement_matrix_grad_;
        MatTend *M_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Matd> *dv_B_;
    DiscreteVariable<VecMatGrad> *dv_displacement_matrix_grad_;
    DiscreteVariable<MatTend> *dv_M_;
};

template <typename... RelationTypes>
class DisplacementMatrixGradient; // preparation for hessian correction matrix

template <typename... Parameters>
class DisplacementMatrixGradient<Inner<Parameters...>>
    : public HessianCorrectionMatrix<Base, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianCorrectionMatrix<Base, Inner<Parameters...>>;

  public:
    explicit DisplacementMatrixGradient(Inner<Parameters...> &inner_relation)
        : BaseDynamicsType(inner_relation) {};
    virtual ~DisplacementMatrixGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };
};

template <typename... Parameters>
class DisplacementMatrixGradient<Contact<Parameters...>>
    : public HessianCorrectionMatrix<Base, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianCorrectionMatrix<Base, Contact<Parameters...>>;

  public:
    explicit DisplacementMatrixGradient(Contact<Parameters...> &contact_relation);
    virtual ~DisplacementMatrixGradient() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
              contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

template <typename... Parameters>
class HessianCorrectionMatrix<Inner<WithUpdate, Parameters...>>
    : public HessianCorrectionMatrix<Base, Inner<Parameters...>>
{
    using BaseDynamicsType = HessianCorrectionMatrix<Base, Inner<Parameters...>>;

  public:
    explicit HessianCorrectionMatrix(Inner<Parameters...> &inner_relation, Real alpha = Real(0))
        : BaseDynamicsType(inner_relation), alpha_(alpha) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit HessianCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : HessianCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~HessianCorrectionMatrix() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };

    class UpdateKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser), alpha_(encloser.alpha_){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real alpha_;
    };

  protected:
    Real alpha_;
};
using HessianCorrectionMatrixInner = HessianCorrectionMatrix<Inner<WithUpdate>>;

template <typename... Parameters>
class HessianCorrectionMatrix<Contact<Parameters...>>
    : public HessianCorrectionMatrix<Base, Contact<Parameters...>>
{
    using BaseDynamicsType = HessianCorrectionMatrix<Base, Contact<Parameters...>>;

  public:
    explicit HessianCorrectionMatrix(Contact<Parameters...> &contact_relation);
    virtual ~HessianCorrectionMatrix() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
              contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};
using HessianCorrectionMatrixComplex = HessianCorrectionMatrix<Inner<WithUpdate>, Contact<>>;
} // namespace SPH
#endif // HESSIAN_CORRECTION_CK_H
