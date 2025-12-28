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
 * @file eulerian_compressible_fluid_integration.h
 * @brief Here, we define the common compressible eulerian classes for fluid dynamics.
 * @author Zhentong Wang and Xiangyu Hu
 */

#ifndef EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H
#define EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H

#include "base_general_dynamics.h"
#include "compressible_fluid.h"
#include "eulerian_riemann_solver.h"
#include "fluid_body.h"
#include "fluid_integration.hpp"
#include "fluid_time_step.h"
#include "viscous_dynamics.hpp"
#include "muscl_hllc_integration.h"
#include <tuple>

namespace SPH
{
namespace fluid_dynamics
{
template <class DataDelegationType>
class BaseIntegrationInCompressibleType : public BaseIntegration<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegrationInCompressibleType(BaseRelationType &relation)
        : BaseIntegration<DataDelegationType>(relation),
          compressible_fluid_(CompressibleFluid(1.0, 1.4)),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          E_(this->particles_->template registerStateVariable<Real>("TotalEnergy")),
          dE_dt_(this->particles_->template registerStateVariable<Real>("TotalEnergyChangeRate")),
          dmass_dt_(this->particles_->template registerStateVariable<Real>("MassChangeRate")),
          mom_(this->particles_->template registerStateVariable<Vecd>("Momentum")),
          force_(this->particles_->template registerStateVariable<Vecd>("Force")),
          force_prior_(this->particles_->template registerStateVariable<Vecd>("ForcePrior")){};
    virtual ~BaseIntegrationInCompressibleType(){};

  protected:
    CompressibleFluid compressible_fluid_;
    Real *Vol_, *E_, *dE_dt_, *dmass_dt_;
    Vecd *mom_, *force_, *force_prior_;
};
using BaseIntegrationInCompressible = BaseIntegrationInCompressibleType<DataDelegateInner>;

template <class RiemannSolverType>
class EulerianCompressibleIntegration1stHalf : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~EulerianCompressibleIntegration1stHalf(){};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using EulerianCompressibleIntegration1stHalfNoRiemann = EulerianCompressibleIntegration1stHalf<NoRiemannSolverInCompressibleEulerianMethod>;
using EulerianCompressibleIntegration1stHalfHLLCRiemann = EulerianCompressibleIntegration1stHalf<HLLCRiemannSolver>;
using EulerianCompressibleIntegration1stHalfHLLCWithLimiterRiemann = EulerianCompressibleIntegration1stHalf<HLLCWithLimiterRiemannSolver>;

/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
 */
template <class RiemannSolverType>
class EulerianCompressibleIntegration2ndHalf : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~EulerianCompressibleIntegration2ndHalf(){};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using EulerianCompressibleIntegration2ndHalfNoRiemann = EulerianCompressibleIntegration2ndHalf<NoRiemannSolverInCompressibleEulerianMethod>;
using EulerianCompressibleIntegration2ndHalfHLLCRiemann = EulerianCompressibleIntegration2ndHalf<HLLCRiemannSolver>;
using EulerianCompressibleIntegration2ndHalfHLLCWithLimiterRiemann = EulerianCompressibleIntegration2ndHalf<HLLCWithLimiterRiemannSolver>;

/**
 * @brief MUSCL-HLLC variant: reconstruct interface states with MUSCL and call HLLC.
 */
template <typename... InteractionTypes>
class EulerianCompressibleIntegration1stHalfMUSCL;

template <>
class EulerianCompressibleIntegration1stHalfMUSCL<Inner<>> : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration1stHalfMUSCL(BaseInnerRelation &inner_relation, const MUSCLHLLCBridgeConfig &bridge_cfg = MUSCLHLLCBridgeConfig());
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianCompressibleIntegration1stHalfMUSCL(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianCompressibleIntegration1stHalfMUSCL(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianCompressibleIntegration1stHalfMUSCL(){};
    // Expose Riemann solver type for compatibility with interfaces expecting it.
    using RiemannSolver = HLLCFSIRiemannSolver;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    MUSCL_HLLC_Bridge bridge_;
    MUSCLHLLCBridgeConfig bridge_cfg_;
    Vecd *rho_grad_;
    Matd *vel_grad_;
    Vecd *p_grad_;
};

template <>
class EulerianCompressibleIntegration1stHalfMUSCL<Contact<Wall>> : public InteractionWithWall<BaseIntegrationInCompressibleType>
{
  public:
    explicit EulerianCompressibleIntegration1stHalfMUSCL(BaseContactRelation &contact_relation, const MUSCLHLLCBridgeConfig &bridge_cfg = MUSCLHLLCBridgeConfig());
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianCompressibleIntegration1stHalfMUSCL(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianCompressibleIntegration1stHalfMUSCL(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianCompressibleIntegration1stHalfMUSCL(){};
    using RiemannSolver = HLLCFSIRiemannSolver;
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    MUSCL_HLLC_Bridge bridge_;
    MUSCLHLLCBridgeConfig bridge_cfg_;
    Vecd *rho_grad_;
    Matd *vel_grad_;
    Vecd *p_grad_;
};

template <typename... InteractionTypes>
class EulerianCompressibleIntegration2ndHalfMUSCL;

template <>
class EulerianCompressibleIntegration2ndHalfMUSCL<Inner<>> : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration2ndHalfMUSCL(BaseInnerRelation &inner_relation, const MUSCLHLLCBridgeConfig &bridge_cfg = MUSCLHLLCBridgeConfig());
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianCompressibleIntegration2ndHalfMUSCL(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianCompressibleIntegration2ndHalfMUSCL(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianCompressibleIntegration2ndHalfMUSCL(){};
    // Expose Riemann solver type for interfaces such as PressureForceFromFluid.
    using RiemannSolver = HLLCFSIRiemannSolver;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    MUSCL_HLLC_Bridge bridge_;
    MUSCLHLLCBridgeConfig bridge_cfg_;
    Vecd *rho_grad_;
    Matd *vel_grad_;
    Vecd *p_grad_;
};

template <>
class EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>> : public InteractionWithWall<BaseIntegrationInCompressibleType>
{
  public:
    explicit EulerianCompressibleIntegration2ndHalfMUSCL(BaseContactRelation &contact_relation, const MUSCLHLLCBridgeConfig &bridge_cfg = MUSCLHLLCBridgeConfig());
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianCompressibleIntegration2ndHalfMUSCL(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianCompressibleIntegration2ndHalfMUSCL(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianCompressibleIntegration2ndHalfMUSCL(){};
    using RiemannSolver = HLLCFSIRiemannSolver;
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    MUSCL_HLLC_Bridge bridge_;
    MUSCLHLLCBridgeConfig bridge_cfg_;
    Vecd *rho_grad_;
    Matd *vel_grad_;
    Vecd *p_grad_;
};

class CompressibleFluidInitialCondition : public FluidInitialCondition
{
  public:
    explicit CompressibleFluidInitialCondition(SPHBody &sph_body);

  protected:
    Vecd *mom_;
    Real *rho_, *Vol_, *mass_, *p_, *E_;
};

class EulerianCompressibleAcousticTimeStepSize : public AcousticTimeStep
{
  protected:
    Real *rho_, *p_;
    Vecd *vel_;
    Real smoothing_length_;
    Fluid &fluid_;

  public:
    explicit EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~EulerianCompressibleAcousticTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
};

// With-wall aliases for MUSCL (inner + wall contact)
using EulerianCompressibleIntegration1stHalfMUSCLWithWall =
    ComplexInteraction<EulerianCompressibleIntegration1stHalfMUSCL<Inner<>, Contact<Wall>>>;
using EulerianCompressibleIntegration2ndHalfMUSCLWithWall =
    ComplexInteraction<EulerianCompressibleIntegration2ndHalfMUSCL<Inner<>, Contact<Wall>>>;
    
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H