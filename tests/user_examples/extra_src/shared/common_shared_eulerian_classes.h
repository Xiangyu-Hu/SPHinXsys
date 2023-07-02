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
 * @file 	common_shared_eulerian_classes.h
 * @brief 	Here, we define the common shared eulerian classes for compressible and weakly compressible fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef COMMON_SHARED_EULERIAN_CLASSES_H
#define COMMON_SHARED_EULERIAN_CLASSES_H

#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
namespace SPH
{
/**
 * @class EulerianFluidBody
 * @brief Eulerian Fluid body uses smoothing length to particle spacing 1.3
 */
class EulerianFluidBody : public FluidBody
{
  public:
    explicit EulerianFluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr) : FluidBody(system, shape_ptr)
    {
        defineAdaptation<SPHAdaptation>(1.3);
    };
    virtual ~EulerianFluidBody(){};
    virtual EulerianFluidBody *ThisObjectPtr() override { return this; };
};

/**
 * @class KernelGradientWithCorrectionInner
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient
 */
class KernelGradientWithCorrectionInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    KernelGradientWithCorrectionInner(BaseInnerRelation &inner_relation);
    virtual ~KernelGradientWithCorrectionInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> B_, local_configuration_inner_;
};

/**
 * @class KernelGradientWithCorrectionComplex
 * @brief obtain the corrected initial configuration in strong form and correct kernel gradient in complex topology
 */
class KernelGradientWithCorrectionComplex : public BaseInteractionComplex<KernelGradientWithCorrectionInner, GeneralDataDelegateContactOnly>
{
  public:
    template <typename... Args>
    KernelGradientWithCorrectionComplex(Args &&...args)
        : BaseInteractionComplex<KernelGradientWithCorrectionInner, GeneralDataDelegateContactOnly>(std::forward<Args>(args)...){};
    virtual ~KernelGradientWithCorrectionComplex(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // COMMON_SHARED_EULERIAN_CLASSES_H