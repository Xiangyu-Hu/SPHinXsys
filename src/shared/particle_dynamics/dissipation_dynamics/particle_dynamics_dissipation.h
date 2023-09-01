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
 * @file 	particle_dynamics_dissipation.h
 * @brief 	Here are the classes for damping the magnitude of
 * 			any variables.
 * 			Note that, currently, these classes works only in single resolution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DISSIPATION_H
#define PARTICLE_DYNAMICS_DISSIPATION_H

#include "all_particle_dynamics.h"

namespace SPH
{
class SolidParticles;

typedef DataDelegateInner<BaseParticles> DissipationDataInner;
typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
    DissipationDataContact;
typedef DataDelegateContact<BaseParticles, SolidParticles, DataDelegateEmptyBase>
    DissipationDataWithWall;

template <typename VariableType>
struct ErrorAndParameters
{
    VariableType error_;
    Real a_, c_;
    ErrorAndParameters() : error_(ZeroData<VariableType>::value), a_(0), c_(0){};
};

/**
 * @class DampingBySplittingAlgorithm
 * @brief A quantity damping by splitting scheme
 * this method modifies the quantity directly.
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */
template <typename VariableType>
class DampingBySplittingInner : public LocalDynamics, public DissipationDataInner
{
  protected:
  public:
    DampingBySplittingInner(BaseInnerRelation &inner_relation, const std::string &variable_name, Real eta);
    virtual ~DampingBySplittingInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real eta_; /**< damping coefficient */
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<VariableType> &variable_;

    virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters);
};

template <typename VariableType>
class DampingBySplittingComplex : public DampingBySplittingInner<VariableType>, public DissipationDataContact
{
  public:
    DampingBySplittingComplex(ComplexRelation &complex_relation, const std::string &variable_name, Real eta);
    virtual ~DampingBySplittingComplex(){};

  protected:
    virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters) override;

  private:
    StdVec<StdLargeVec<Real> *> contact_Vol_, contact_mass_;
    StdVec<StdLargeVec<VariableType> *> contact_variable_;
};

template <typename VariableType,
          template <typename BaseVariableType>
          class BaseDampingBySplittingType>
class DampingBySplittingWithWall : public BaseDampingBySplittingType<VariableType>, public DissipationDataWithWall
{
  public:
    DampingBySplittingWithWall(ComplexRelation &complex_wall_relation, const std::string &variable_name, Real eta);
    virtual ~DampingBySplittingWithWall(){};

  protected:
    virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;

  private:
    StdVec<StdLargeVec<Real> *> wall_Vol_;
    StdVec<StdLargeVec<VariableType> *> wall_variable_;
};

/**
 * @class DampingPairwiseInner
 * @brief A quantity damping by a pairwise splitting scheme
 * this method modifies the quantity directly
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */
template <typename VariableType>
class DampingPairwiseInner : public LocalDynamics, public DissipationDataInner
{
  public:
    DampingPairwiseInner(BaseInnerRelation &inner_relation, const std::string &variable_name, Real eta);
    virtual ~DampingPairwiseInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<VariableType> &variable_;
    Real eta_; /**< damping coefficient */
};

template <typename VariableType>
class DampingPairwiseComplex : public DampingPairwiseInner<VariableType>, public DissipationDataContact
{
  public:
    DampingPairwiseComplex(BaseInnerRelation &inner_relation,
                           BaseContactRelation &contact_relation, const std::string &variable_name, Real eta);
    DampingPairwiseComplex(ComplexRelation &complex_relation, const std::string &variable_name, Real eta);
    virtual ~DampingPairwiseComplex(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  private:
    StdVec<StdLargeVec<Real> *> contact_Vol_, contact_mass_;
    StdVec<StdLargeVec<VariableType> *> contact_variable_;
};

/**
 * @class DampingPairwiseWithWall
 * @brief Damping with wall by which the wall velocity is not updated
 * and the mass of wall particle is not considered.
 */
template <typename VariableType,
          template <typename BaseVariableType> class BaseDampingPairwiseType>
class DampingPairwiseWithWall : public BaseDampingPairwiseType<VariableType>,
                                public DissipationDataWithWall
{
  public:
    DampingPairwiseWithWall(BaseInnerRelation &inner_relation,
                            BaseContactRelation &contact_relation, const std::string &variable_name, Real eta);
    DampingPairwiseWithWall(ComplexRelation &complex_wall_relation, const std::string &variable_name, Real eta);
    virtual ~DampingPairwiseWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  private:
    StdVec<StdLargeVec<Real> *> wall_Vol_;
    StdVec<StdLargeVec<VariableType> *> wall_variable_;
};

/**
 * @class DampingPairwiseFromWall
 * @brief Damping to wall by which the wall velocity is not updated
 * and the mass of wall particle is not considered.
 */
template <typename VariableType>
class DampingPairwiseFromWall : public LocalDynamics,
                                public DataDelegateContact<BaseParticles, SolidParticles>
{
  public:
    DampingPairwiseFromWall(BaseContactRelation &contact_relation, const std::string &variable_name, Real eta);
    virtual ~DampingPairwiseFromWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  private:
    Real eta_; /**< damping coefficient */
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<VariableType> &variable_;
    StdVec<StdLargeVec<Real> *> wall_Vol_;
    StdVec<StdLargeVec<VariableType> *> wall_variable_;
};

/**
 * @class DampingWithRandomChoice
 * @brief A random choice method for obtaining static equilibrium state
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */
template <class DampingAlgorithmType>
class DampingWithRandomChoice : public DampingAlgorithmType
{
  protected:
    Real random_ratio_;
    bool RandomChoice();

  public:
    template <typename... ConstructorArgs>
    DampingWithRandomChoice(Real random_ratio, ConstructorArgs &&...args);
    virtual ~DampingWithRandomChoice(){};

    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // PARTICLE_DYNAMICS_DISSIPATION_H