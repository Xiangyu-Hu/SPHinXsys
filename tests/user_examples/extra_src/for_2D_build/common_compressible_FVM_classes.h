/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2023 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	common_compressible_FVM_classes.h
 * @brief 	Here, we define the common compressible classes for fluid dynamics in FVM.
 * @author	Zhentong Wang and Xiangyu Hu
 */

#ifndef COMMON_COMPRESSIBLE_FVM_CLASSES_H
#define COMMON_COMPRESSIBLE_FVM_CLASSES_H
#include "common_compressible_eulerian_classes.hpp"
#include "common_shared_FVM_classes.h"
namespace SPH
{
/**
 * @class CompressibleAcousticTimeStepSizeInFVM
 * @brief Computing the acoustic time step size
 */
class CompressibleAcousticTimeStepSizeInFVM : public fluid_dynamics::AcousticTimeStepSize
{
  protected:
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_;
    Real max_distance_between_nodes_;

  public:
    explicit CompressibleAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real max_distance_between_nodes, Real acousticCFL = 0.6);
    virtual ~CompressibleAcousticTimeStepSizeInFVM(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    CompressibleFluid compressible_fluid_;
    Real acousticCFL_;
};
} // namespace SPH
#endif // COMMON_COMPRESSIBLE_FVM_CLASSES_H