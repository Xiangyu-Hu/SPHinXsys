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
 * @file kernel_correction_ck.h
 * @brief This the methods related on the corrections of SPH smoothing kernel.
 * @details The corrections aim to increase the numerical consistency
 * or accuracy for kernel approximations.
 * @author Xiangyu Hu
 */

#ifndef KERNEL_CORRECTION_CK_H
#define KERNEL_CORRECTION_CK_H

#include "base_general_dynamics.h"
#include "interaction_ck.h"

namespace SPH
{
template <typename... RelationTypes>
class LinearCorrectionMatrix;

template <template <typename...> class RelationType, typename... Parameters>
class LinearCorrectionMatrix<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{

  public:
    template <class DynamicsIdentifier>
    explicit LinearCorrectionMatrix(DynamicsIdentifier &identifier);
    virtual ~LinearCorrectionMatrix() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       LinearCorrectionMatrix<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

      protected:
        Real *Vol_;
        Matd *B_;
    };

  protected:
    DiscreteVariable<Matd> *dv_B_;
};

template <typename... Parameters>
class LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>>
    : public LinearCorrectionMatrix<Base, Inner<Parameters...>>
{

  public:
    explicit LinearCorrectionMatrix(Inner<Parameters...> &inner_relation, Real alpha = Real(0))
        : LinearCorrectionMatrix<Base, Inner<Parameters...>>(inner_relation), alpha_(alpha) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit LinearCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : LinearCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~LinearCorrectionMatrix() {};

    class InteractKernel
        : public LinearCorrectionMatrix<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>> &encloser);
        void interact(size_t index_i, Real dt = 0.0);
    };

    class UpdateKernel
        : public LinearCorrectionMatrix<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     LinearCorrectionMatrix<Inner<WithUpdate, Parameters...>> &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real alpha_;
    };

  protected:
    Real alpha_;
};
using LinearCorrectionMatrixInner = LinearCorrectionMatrix<Inner<WithUpdate>>;

template <typename... Parameters>
class LinearCorrectionMatrix<Contact<Parameters...>>
    : public LinearCorrectionMatrix<Base, Contact<Parameters...>>
{
  public:
    explicit LinearCorrectionMatrix(Contact<Parameters...> &contact_relation);
    virtual ~LinearCorrectionMatrix() {};

    class InteractKernel
        : public LinearCorrectionMatrix<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       LinearCorrectionMatrix<Contact<Parameters...>> &encloser,
                       size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_k_;
    };
};

using LinearCorrectionMatrixComplex = LinearCorrectionMatrix<Inner<WithUpdate>, Contact<>>;

class NoKernelCorrectionCK : public KernelCorrection
{
  public:
    typedef Real CorrectionDataType;
    NoKernelCorrectionCK(BaseParticles *particles) : KernelCorrection() {};

    class ComputingKernel : public ParameterFixed<Real>
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        NoKernelCorrectionCK &encloser)
            : ParameterFixed<Real>(1.0){};
    };
};

class LinearCorrectionCK : public KernelCorrection
{
  public:
    typedef Matd CorrectionDataType;
    LinearCorrectionCK(BaseParticles *particles)
        : KernelCorrection(),
          dv_B_(particles->getVariableByName<Matd>("LinearCorrectionMatrix")) {};

    class ComputingKernel : public ParameterVariable<Matd>
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, LinearCorrectionCK &encloser)
            : ParameterVariable<Matd>(encloser.dv_B_->DelegatedData(ex_policy)){};
    };

  protected:
    DiscreteVariable<Matd> *dv_B_;
};

} // namespace SPH
#endif // KERNEL_CORRECTION_CK_H
