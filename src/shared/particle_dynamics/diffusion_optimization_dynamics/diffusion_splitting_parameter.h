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
 * @file 	diffusion_splitting_parameter.h
 * @brief    This is the splitting method for solving parameter field in optimization problem.
 * Note that here inner interaction and that with boundary are derived from the inner interaction.
 * This is because the error and parameters are computed based on both.
 * @author   Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_SPLITTING_PARAMETER_H
#define DIFFUSION_SPLITTING_PARAMETER_H

#include "diffusion_splitting_base.h"
#include "diffusion_splitting_state.h"

namespace SPH
{
/**
 * @class ParameterSplittingByPDEInner
 * @brief Modify the parameter inner by splitting operator based on PDEs.
 */
template <typename VariableType>
class ParameterSplittingByPDEInner
    : public OptimizationBySplittingAlgorithmBase<VariableType>
{
  public:
    ParameterSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &variable_name);
    virtual ~ParameterSplittingByPDEInner(){};

  protected:
    virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters);
    virtual void interaction(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class ParameterSplittingByPDEWithBoundary
 * @brief Modify the parameter contact with the boundary by splitting operator based on PDEs.
 */
template <typename VariableType>
class ParameterSplittingByPDEWithBoundary
    : public ParameterSplittingByPDEInner<VariableType>,
      public DataDelegateContactOnly
{
  public:
    ParameterSplittingByPDEWithBoundary(BaseInnerRelation &inner_relation,
                                        BaseContactRelation &contact_relation, const std::string &variable_name);
    virtual ~ParameterSplittingByPDEWithBoundary(){};

  protected:
    StdVec<StdLargeVec<Vecd> *> boundary_normal_vector_;
    StdVec<StdLargeVec<Real> *> boundary_heat_flux_, boundary_Vol_;
    StdVec<StdLargeVec<Real> *> boundary_species_;
    virtual ErrorAndParameters<VariableType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class UpdateParameterPDEResidual
 * @brief Update the global residual from parameter after one splitting loop finished.
 */
template <typename ParameterSplittingType>
class UpdateParameterPDEResidual : public ParameterSplittingType
{
  public:
    template <typename... Args>
    UpdateParameterPDEResidual(Args &&...args);
    virtual ~UpdateParameterPDEResidual(){};

  protected:
    virtual void interaction(size_t index_i, Real dt = 0.0) override;
};
} // namespace SPH

#endif // DIFFUSION_SPLITTING_PARAMETER_H
