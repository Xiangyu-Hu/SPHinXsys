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
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
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
#include "common_compressible_eulerian_classes.h"
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

      public:
        explicit CompressibleAcousticTimeStepSizeInFVM(SPHBody &sph_body);
        virtual ~CompressibleAcousticTimeStepSizeInFVM(){};
        Real reduce(size_t index_i, Real dt = 0.0);
        virtual Real outputResult(Real reduced_value) override;
        CompressibleFluid compressible_fluid_;
    };

    /**
	* @class BaseIntegrationInCompressibleFVM
	* @brief Pure abstract base class for all fluid relaxation schemes in compressible flows
	*/
    class BaseIntegrationInCompressibleFVM : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
    {
      public:
        explicit BaseIntegrationInCompressibleFVM(BaseInnerRelationInFVM &inner_relation);
        virtual ~BaseIntegrationInCompressibleFVM(){};

      protected:
        CompressibleFluid compressible_fluid_;
        StdLargeVec<Real> &E_, &dE_dt_, &dE_dt_prior_, &rho_, &drho_dt_, &p_;
        StdLargeVec<Vecd> &mom_, &dmom_dt_, &dmom_dt_prior_, &vel_, &pos_;
    };

    /**
	* @class BaseIntegration1stHalfInFVM
	* @brief Template class for pressure relaxation scheme with the Riemann solver In FVM
	* as template variable
	*/
	template <class RiemannSolverType>
    class BaseIntegration1stHalfInFVM : public BaseIntegrationInCompressibleFVM
	{
	public:
		explicit BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter);
		virtual ~BaseIntegration1stHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void initialization(size_t index_i, Real dt);
		void interaction(size_t index_i, Real dt);
		void update(size_t index_i, Real dt);
	};
	using Integration1stHalfHLLCRiemannInFVM = BaseIntegration1stHalfInFVM<HLLCRiemannSolver>;
	using Integration1stHalfHLLCWithLimiterRiemannInFVM = BaseIntegration1stHalfInFVM<HLLCWithLimiterRiemannSolver>;

	/**
	 * @class BaseIntegration2ndHalfInFVM
	 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter In FVM
	 */
	template <class RiemannSolverType>
	class BaseIntegration2ndHalfInFVM : public BaseIntegrationInCompressibleFVM
	{
	public:
		explicit BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter);
		virtual ~BaseIntegration2ndHalfInFVM() {};
		RiemannSolverType riemann_solver_;
		void interaction(size_t index_i, Real dt);
		void update(size_t index_i, Real dt);
	};
	using Integration2ndHalfHLLCRiemannInFVM = BaseIntegration2ndHalfInFVM<HLLCRiemannSolver>;
	using Integration2ndHalfHLLCWithLimiterRiemannInFVM = BaseIntegration2ndHalfInFVM<HLLCWithLimiterRiemannSolver>;
}
#endif // COMMON_COMPRESSIBLE_FVM_CLASSES_H