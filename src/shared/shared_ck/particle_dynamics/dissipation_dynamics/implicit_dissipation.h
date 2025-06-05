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
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file implicit_dissipation.h
 * @brief TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef IMPLICIT_DISSIPATION_H
#define IMPLICIT_DISSIPATION_H

#include "base_dissipation.h"
#include "general_reduce_ck.h"

namespace SPH
{
// Note: DissipativeTransform method is for obtaining the matrix multiplication
// between the dissipation matrix and a variable, which can be the dissipation
// variable or residue.
template <typename...>
class DissipativeTransform;

template <typename DissipationType, typename... Parameters>
class DissipativeTransform<Inner<DissipationType, Parameters...>>
    : public Dissipation<Base, DissipationType, Inner<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseDissipationType = Dissipation<Base, DissipationType, Inner<Parameters...>>;

  public:
    DissipativeTransform(Inner<Parameters...> &inner_relation, const std::string &variable_name);
    virtual ~DissipativeTransform() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        DataType *transformed_;
    };

  protected:
    DiscreteVariable<DataType> *dv_transformed_;
};

template <typename DataType, class DynamicsIdentifier>
class DissipationRHS : BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    DissipationRHS(DynamicsIdentifier &identifier, const std::string &variable_name);
    virtual ~DissipationRHS() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        DataType *variable_;
        DataType *old_state_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<DataType> *dv_old_state_;
};

template <typename DataType, class DynamicsIdentifier>
class FullDissipationResidue : BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    FullDissipationResidue(DynamicsIdentifier &identifier, const std::string &variable_name);
    virtual ~FullDissipationResidue() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        DataType *residue_;
        DataType *old_state_;
        DataType *transformed_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_old_state_;
    DiscreteVariable<DataType> *dv_transformed_;
};

template <typename DataType, typename TensorProductType, class DynamicsIdentifier>
class UpdateDissipationResidue : BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    UpdateDissipationResidue(DynamicsIdentifier &identifier,
                             const std::string &variable_name,
                             SingularVariable<TensorProductType> *sv_search_depth);
    virtual ~UpdateDissipationResidue() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        DataType *residue_;
        DataType *transformed_;
        TensorProductType *search_depth_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_transformed_;
    SingularVariable<TensorProductType> *sv_search_depth_;
};

template <typename DataType, typename TensorProductType, class DynamicsIdentifier>
class UpdateDissipationSolution : BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    UpdateDissipationSolution(DynamicsIdentifier &identifier,
                              const std::string &variable_name,
                              SingularVariable<TensorProductType> *sv_search_depth);
    virtual ~UpdateDissipationSolution() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real search_depth);

      protected:
        DataType *residue_;
        DataType *variable_;
        TensorProductType *search_depth_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_variable_;
    SingularVariable<TensorProductType> *sv_search_depth_;
};

template <typename DataType, typename TensorProductType, class DynamicsIdentifier>
class DissipationResidueSum
    : public QuantityTensorProductSum<DataType, TensorProductType, DynamicsIdentifier>
{
  public:
    DissipationResidueSum(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductSum<DataType, TensorProductType, DynamicsIdentifier>(
              identifier, "Residue" + variable_name, "Residue" + variable_name) {};
    virtual ~DissipationResidueSum() {};
};

template <typename DataType, typename TensorProductType, class DynamicsIdentifier>
class TransformedDissipationResidueSum
    : public QuantityTensorProductSum<DataType, TensorProductType, DynamicsIdentifier>
{
  public:
    TransformedDissipationResidueSum(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductSum<DataType, TensorProductType, DynamicsIdentifier>(
              identifier, "Residue" + variable_name, "TransformedResidue" + variable_name) {};
    virtual ~TransformedDissipationResidueSum() {};
};

template <typename...>
class ImplicitDissipation;

template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename TensorProductType, typename... Parameters>
class ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, TensorProductType, Parameters...>>
    : public BaseDynamics<void>
{
    using DataType = typename DissipationType::DataType;
    using DynamicsIdentifier = typename RelationType<Parameters...>::Identifier;

  public:
    ImplicitDissipation(RelationType<Parameters...> &first_relation,
                        const std::string &variable_name,
                        Real sqr_norm_criteria);
    template <typename... ControlParameters, typename... RelationParameters, typename... Args>
    auto &addContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args);
    virtual void exec(Real dt = 0.0) override;

  protected:
    Real sqr_norm_criteria_;
    DynamicsIdentifier &dynamics_identifier_;
    SingularVariable<TensorProductType> sv_search_depth_;
    DissipationRHS<DataType, DynamicsIdentifier> dissipation_rhs_;
    InteractionDynamicsCK<
        ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>>
        transformed_variable_;
    InteractionDynamicsCK<
        ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>>
        transformed_residue_;
    StateDynamics<ExecutionPolicy, FullDissipationResidue<DataType, DynamicsIdentifier>>
        full_dissipation_residue_;
    StateDynamics<
        ExecutionPolicy, UpdateDissipationResidue<DataType, TensorProductType, DynamicsIdentifier>>
        update_dissipation_residue_;
    StateDynamics<
        ExecutionPolicy, DissipationResidueSum<DataType, TensorProductType, DynamicsIdentifier>>
        dissipation_residue_sum_;
    StateDynamics<
        ExecutionPolicy, TransformedDissipationResidueSum<DataType, TensorProductType, DynamicsIdentifier>>
        transformed_dissipation_residue_sum_;
    StateDynamics<
        ExecutionPolicy, UpdateDissipationSolution<DataType, TensorProductType, DynamicsIdentifier>>
        update_dissipation_solution_;
};
} // namespace SPH
#endif // IMPLICIT_DISSIPATION_H