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
 * @file 	diffusion_splitting_state.h
 * @brief 	This is the splitting method for solving state field in optimization problem.
 * Note that here inner interaction and that with boundary are derived from the inner interaction.
 * This is because the error and parameters are computed based on both.
 * @author   Bo Zhang and Xiangyu H
 */

#ifndef DIFFUSION_SPLITTING_STATE_H
#define DIFFUSION_SPLITTING_STATE_H

#include "diffusion_splitting_base.h"
#include "diffusion_splitting_parameter.h"

namespace SPH
{
/**
 * @class TemperatureSplittingByPDEInner
 * @brief The temperature on each particle will be modified innerly to satisfy the PDEs.
 */
template <typename DataType>
class TemperatureSplittingByPDEInner
    : public OptimizationBySplittingAlgorithmBase<DataType>
{
  public:
    TemperatureSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &variable_name);
    virtual ~TemperatureSplittingByPDEInner(){};

  protected:
    virtual ErrorAndParameters<DataType> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters);
    virtual void interaction(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class TemperatureSplittingByPDEWithBoundary
 * @brief The temperature on each particle will be modified with boundary to satisfy the PDEs.
 */
template <typename DataType>
class TemperatureSplittingByPDEWithBoundary
    : public TemperatureSplittingByPDEInner<DataType>,
      public DataDelegateContact
{
  public:
    TemperatureSplittingByPDEWithBoundary(BaseInnerRelation &inner_relation,
                                          BaseContactRelation &contact_relation,
                                          const std::string &variable_name);
    virtual ~TemperatureSplittingByPDEWithBoundary(){};

  protected:
    StdVec<DataType *> boundary_variable_;
    StdVec<Real *> boundary_heat_flux_, boundary_Vol_;
    StdVec<Vecd *> boundary_normal_vector_;
    virtual ErrorAndParameters<DataType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class UpdateTemperaturePDEResidual
 * @brief Update the global residual from temperature after one splitting loop finished.
 */
template <typename TemperatureSplittingType>
class UpdateTemperaturePDEResidual : public TemperatureSplittingType
{
  public:
    template <typename... Args>
    UpdateTemperaturePDEResidual(Args &&...args);
    virtual ~UpdateTemperaturePDEResidual(){};

  protected:
    virtual void interaction(size_t index_i, Real dt = 0.0) override;
};
} // namespace SPH

#endif // DIFFUSION_SPLITTING_STATE_H