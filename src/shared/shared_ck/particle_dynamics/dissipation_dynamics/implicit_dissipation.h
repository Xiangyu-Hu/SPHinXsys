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
 * @brief Implicit dissipation methods are for obtaining solution for both for steady and transient problems
 * of dissipation system. their numerical stability is dependent on the condition number of the dissipative
 * operators. Generally, an implicit method is less accurate than explicit methods.
 * Therefore, it is used when explicit time-step size is too small to be practical,
 * or when the problem is too stiff for explicit methods..
 * @details Here, we use conjugate gradient (CG) method as the starting point of implicit solver.
 * Other methods, especially Chebyshev iteration methods should be considered as it is said has better
 * performance on GPUs.
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
class QuantityTensorProductAverage
    : public BaseLocalDynamicsReduce<
          ReduceSum<std::pair<TensorProductType<DataType>, Real>>, DynamicsIdentifier>
{
    using ReduceReturnType = std::pair<TensorProductType<DataType>, Real>;
    using BaseReduceDynamics = BaseLocalDynamicsReduce<
        ReduceSum<std::pair<TensorProductType<DataType>, Real>>, DynamicsIdentifier>;

  public:
    QuantityTensorProductAverage(DynamicsIdentifier &identifier,
                                 const std::string &variable_name1,
                                 const std::string &variable_name2);
    virtual ~QuantityTensorProductAverage() {};

    class FinishDynamics
    {
      public:
        using OutputType = TensorProductType<DataType>;
        template <class EncloserType>
        FinishDynamics(EncloserType &encloser){};
        OutputType Result(const ReduceReturnType &reduced_value)
        {
            return reduced_value.first / reduced_value.second;
        }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            return ReduceReturnType(
                tensorProduct(variable1_[index_i], variable2_[index_i]), Real(1));
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

template <typename DissipationType, typename... Parameters>
class DissipativeTransform<Contact<Dirichlet<DissipationType>, Parameters...>>
    : public Dissipation<Base, DissipationType, Contact<Parameters...>>
{
    using DataType = typename DissipationType::DataType;
    using BaseDissipationType = Dissipation<Base, DissipationType, Contact<Parameters...>>;

  public:
    DissipativeTransform(Contact<Parameters...> &contact_relation, const std::string &variable_name);
    virtual ~DissipativeTransform() {};

    class InteractKernel : public BaseDissipationType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        DataType *transformed_;
        Real *contact_Vol_;
        DataType *contact_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_transformed_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariable<DataType> *> contact_dv_variable_;
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
        void update(size_t index_i, Real dt = 0.0)
        {
            rhs_[index_i] = variable_[index_i];
        };

      protected:
        DataType *variable_;
        DataType *rhs_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<DataType> *dv_rhs_;
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
        void update(size_t index_i, Real dt = 0.0)
        {
            residue_[index_i] = rhs_[index_i] - transformed_[index_i];
        };

      protected:
        DataType *residue_;
        DataType *rhs_;
        DataType *transformed_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_rhs_;
    DiscreteVariable<DataType> *dv_transformed_;
};

template <typename DataType, class DynamicsIdentifier>
class UpdateDissipationSearch : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    UpdateDissipationSearch(DynamicsIdentifier &identifier,
                            const std::string &variable_name,
                            SingularVariable<TensorProductType<DataType>> *sv_residue_ratio);
    virtual ~UpdateDissipationSearch() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0)
        {
            search_[index_i] = residue_[index_i] + (*residue_ratio_) * search_[index_i];
        };

      protected:
        DataType *residue_;
        DataType *search_;
        TensorProductType<DataType> *residue_ratio_;
    };

  protected:
    DiscreteVariable<DataType> *dv_residue_;
    DiscreteVariable<DataType> *dv_search_;
    SingularVariable<TensorProductType<DataType>> *sv_residue_ratio_;
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
        void update(size_t index_i, Real search_depth)
        {
            variable_[index_i] += (*search_depth_) * search_[index_i];
        };

      protected:
        DataType *search_;
        DataType *variable_;
        TensorProductType<DataType> *search_depth_;
    };

  protected:
    DiscreteVariable<DataType> *dv_search_;
    DiscreteVariable<DataType> *dv_variable_;
    SingularVariable<TensorProductType<DataType>> *sv_search_depth_;
};

template <typename DataType, class DynamicsIdentifier>
class DissipationResidueAverage
    : public QuantityTensorProductAverage<DataType, DynamicsIdentifier>
{
  public:
    DissipationResidueAverage(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductAverage<DataType, DynamicsIdentifier>(
              identifier, "Residue" + variable_name, "Residue" + variable_name) {};
    virtual ~DissipationResidueAverage() {};
};

template <typename DataType, class DynamicsIdentifier>
class TransformedDissipationSearchAverage
    : public QuantityTensorProductAverage<DataType, DynamicsIdentifier>
{
  public:
    TransformedDissipationSearchAverage(DynamicsIdentifier &identifier, const std::string &variable_name)
        : QuantityTensorProductAverage<DataType, DynamicsIdentifier>(
              identifier, "Search" + variable_name, "TransformedSearch" + variable_name) {};
    virtual ~TransformedDissipationSearchAverage() {};
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
    UniquePtrsKeeper<BaseDynamics<void>> initialization_ptrs_keeper_;

  public:
    ImplicitDissipation(RelationType<Parameters...> &first_relation,
                        const std::string &variable_name,
                        Real convergence_criteria);
    template <typename... ControlParameters, typename... RelationParameters, typename... Args>
    auto &addPostContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args);
    void initializeImplicitDissipation();
    virtual void exec(Real dt = 0.0) override;

  protected:
    std::string variable_name_;
    Real convergence_criteria_;
    DynamicsIdentifier &identifier_;
    StdVec<BaseDynamics<void> *> initialization_methods_;
    TensorProductType<DataType> zero_residue_ratio_;
    SingularVariable<TensorProductType<DataType>> sv_search_depth_, sv_residue_ratio_;
    StateDynamics<ExecutionPolicy, DissipationRHS<DataType, DynamicsIdentifier>> dissipation_rhs_;
    InteractionDynamicsCK<ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>> transformed_variable_;
    StateDynamics<ExecutionPolicy, FullDissipationResidue<DataType, DynamicsIdentifier>> full_residue_;
    StateDynamics<ExecutionPolicy, VariableAssignment<DynamicsIdentifier, CopyVariable<DataType>>> initial_search_;
    InteractionDynamicsCK<ExecutionPolicy, DissipativeTransform<RelationType<DissipationType, Parameters...>>> transformed_search_;
    StateDynamics<ExecutionPolicy, UpdateDissipationSearch<DataType, DynamicsIdentifier>> update_search_;
    ReduceDynamicsCK<ExecutionPolicy, DissipationResidueAverage<DataType, DynamicsIdentifier>> residue_average_;
    ReduceDynamicsCK<ExecutionPolicy, TransformedDissipationSearchAverage<DataType, DynamicsIdentifier>> transformed_search_average_;
    StateDynamics<ExecutionPolicy, UpdateDissipationSolution<DataType, DynamicsIdentifier>> update_dissipation_solution_;
};
} // namespace SPH
#endif // IMPLICIT_DISSIPATION_H