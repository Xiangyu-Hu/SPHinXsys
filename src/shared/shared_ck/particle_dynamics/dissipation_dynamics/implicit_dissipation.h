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
template <typename...>
struct TensorProductReturnType;

template <>
struct TensorProductReturnType<Real>
{
    typedef Real type;
};

template <int N>
struct TensorProductReturnType<Eigen::Matrix<Real, N, 1>>
{
    typedef Eigen::Matrix<Real, N, N> type;
};

template <typename DataType>
using TensorProductType = typename TensorProductReturnType<DataType>::type;

template <typename DataType, class DynamicsIdentifier>
class QuantityTensorProductSum
    : public BaseLocalDynamicsReduce<
          ReduceSum<TensorProductType<DataType>>, DynamicsIdentifier>
{
    using BaseReduceDynamics = BaseLocalDynamicsReduce<
        ReduceSum<TensorProductType<DataType>>, DynamicsIdentifier>;

  public:
    QuantityTensorProductSum(DynamicsIdentifier &identifier,
                             const std::string &variable_name1,
                             const std::string &variable_name2);
    virtual ~QuantityTensorProductSum() {};

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        TensorProductType<DataType> reduce(size_t index_i, Real dt = 0.0)
        {
            return tensorProduct(variable1_[index_i], variable2_[index_i]);
        };

      protected:
        DataType *variable1_;
        DataType *variable2_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable1_;
    DiscreteVariable<DataType> *dv_variable2_;
};

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
class DissipationRHS : public BaseLocalDynamics<DynamicsIdentifier>
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
class FullDissipationResidue : public BaseLocalDynamics<DynamicsIdentifier>
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

template <typename DataType, class DynamicsIdentifier>
class UpdateDissipationResidue : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    UpdateDissipationResidue(DynamicsIdentifier &identifier,
                             const std::string &variable_name,
                             SingularVariable<TensorProductType<DataType>> *sv_search_depth);
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
        TensorProductType<DataType> *search_depth_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_transformed_;
    SingularVariable<TensorProductType<DataType>> *sv_search_depth_;
};

template <typename DataType, class DynamicsIdentifier>
class UpdateDissipationSolution : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    UpdateDissipationSolution(DynamicsIdentifier &identifier,
                              const std::string &variable_name,
                              SingularVariable<TensorProductType<DataType>> *sv_search_depth);
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
        TensorProductType<DataType> *search_depth_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_variable_;
    SingularVariable<TensorProductType<DataType>> *sv_search_depth_;
};

template <typename DataType, class DynamicsIdentifier>
class DissipationResidueSum : public QuantityTensorProductSum<DataType, DynamicsIdentifier>
{
  public:
    DissipationResidueSum(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductSum<DataType, DynamicsIdentifier>(
              identifier, "Residue" + variable_name, "Residue" + variable_name) {};
    virtual ~DissipationResidueSum() {};
};

template <typename DataType, class DynamicsIdentifier>
class TransformedDissipationResidueSum
    : public QuantityTensorProductSum<DataType, DynamicsIdentifier>
{
  public:
    TransformedDissipationResidueSum(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductSum<DataType, DynamicsIdentifier>(
              identifier, "Residue" + variable_name, "TransformedResidue" + variable_name) {};
    virtual ~TransformedDissipationResidueSum() {};
};

template <typename...>
class ImplicitDissipation;

template <class ExecutionPolicy, template <typename...> class RelationType,
          typename DissipationType, typename... Parameters>
class ImplicitDissipation<ExecutionPolicy, RelationType<DissipationType, Parameters...>>
    : public BaseDynamics<void>
{
    using DataType = typename DissipationType::DataType;
    using DynamicsIdentifier = typename RelationType<Parameters...>::SourceType;

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
    SingularVariable<TensorProductType<DataType>> sv_search_depth_;
    SingularVariable<UnsignedInt> *sv_total_real_particles_;
    StateDynamics<ExecutionPolicy, DissipationRHS<DataType, DynamicsIdentifier>> dissipation_rhs_;
    InteractionDynamicsCK<ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>> transformed_variable_;
    StateDynamics<ExecutionPolicy, FullDissipationResidue<DataType, DynamicsIdentifier>> full_residue_;
    InteractionDynamicsCK<ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>> transformed_residue_;
    StateDynamics<ExecutionPolicy, UpdateDissipationResidue<DataType, DynamicsIdentifier>> update_residue_;
    ReduceDynamicsCK<ExecutionPolicy, DissipationResidueSum<DataType, DynamicsIdentifier>> dissipation_residue_sum_;
    ReduceDynamicsCK<ExecutionPolicy, TransformedDissipationResidueSum<DataType, DynamicsIdentifier>> transformed_dissipation_residue_sum_;
    StateDynamics<ExecutionPolicy, UpdateDissipationSolution<DataType, DynamicsIdentifier>> update_dissipation_solution_;
};
} // namespace SPH
#endif // IMPLICIT_DISSIPATION_H