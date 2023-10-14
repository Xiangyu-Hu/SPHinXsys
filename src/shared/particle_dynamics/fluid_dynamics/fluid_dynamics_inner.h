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
 * @file 	fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_INNER_H
#define FLUID_DYNAMICS_INNER_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "fluid_body.h"
#include "riemann_solver.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateSimple<BaseParticles> FluidDataSimple;
typedef DataDelegateInner<BaseParticles> FluidDataInner;

/**
 * @class FluidInitialCondition
 * @brief  Set initial condition for a fluid body.
 * This is a abstract class to be override for case specific initial conditions
 */
class FluidInitialCondition : public LocalDynamics, public FluidDataSimple
{
  public:
    explicit FluidInitialCondition(SPHBody &sph_body);
    virtual ~FluidInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
};

class BaseDensitySummationInnerKernel {
  public:
    BaseDensitySummationInnerKernel(NeighborhoodDevice* inner_configuration, BaseParticles* particles,
                                    DeviceReal rho0, DeviceReal invSigma0) :
        inner_configuration_(inner_configuration), rho_(particles->getDeviceVariableByName<DeviceReal>("Density")),
        rho_sum_(particles->registerDeviceVariable<DeviceReal>("DensitySummation", particles->total_real_particles_)),
        mass_(particles->getDeviceVariableByName<DeviceReal>("Mass")), rho0_(rho0), inv_sigma0_(invSigma0) {}
  protected:
    NeighborhoodDevice* inner_configuration_;
    DeviceReal *rho_, *rho_sum_, *mass_;
    DeviceReal rho0_, inv_sigma0_;
};

/**
 * @class BaseDensitySummationInner
 * @brief Base class for computing density by summation
 */
class BaseDensitySummationInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseDensitySummationInner(BaseInnerRelation &inner_relation);
    virtual ~BaseDensitySummationInner(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_, rho_sum_, &mass_;
    Real rho0_, inv_sigma0_;
};

class DensitySummationInnerKernel : public BaseDensitySummationInnerKernel {
  public:
    template<class ...Args>
    DensitySummationInnerKernel(DeviceReal W0, Args ...baseArgs) :
        BaseDensitySummationInnerKernel(std::forward<Args>(baseArgs)...),  W0_(W0) {}

    template<class RealType, class NeighborhoodType>
    static void interaction(size_t index_i, Real dt, NeighborhoodType* inner_configuration, RealType W0,
                            RealType* rho_sum, RealType rho0, RealType inv_sigma0) {
        RealType sigma = W0;
        const auto& inner_neighborhood = inner_configuration[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size(); ++n)
            sigma += inner_neighborhood.W_ij_[n];
        rho_sum[index_i] = sigma * rho0 * inv_sigma0;
    }

    void interaction(size_t index_i, Real dt = 0.0) {
        interaction(index_i, dt, inner_configuration_, W0_, rho_sum_, rho0_, inv_sigma0_);
    }

  private:
    DeviceReal W0_;
};

/**
 * @class DensitySummationInner
 * @brief  computing density by summation
 */
class DensitySummationInner : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInner(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    auto& getDeviceProxy() {
        return device_proxy;
    }

  protected:
    Real W0_;
    ExecutionProxy<DensitySummationInner, DensitySummationInnerKernel> device_proxy;
};

/**
 * @class DensitySummationInnerAdaptive
 * @brief  computing density by summation with variable smoothing length
 */
class DensitySummationInnerAdaptive : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInnerAdaptive(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class BaseViscousAccelerationInner
 * @brief Base class for the viscosity force induced acceleration
 */
class BaseViscousAccelerationInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseViscousAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~BaseViscousAccelerationInner(){};

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &vel_, &acc_prior_;
    Real mu_;
    Real smoothing_length_;
};

/**
 * @class ViscousAccelerationInner
 * @brief  the viscosity force induced acceleration
 */
class ViscousAccelerationInner : public BaseViscousAccelerationInner
{
  public:
    explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAccelerationInner(inner_relation){};
    virtual ~ViscousAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class AngularConservativeViscousAccelerationInner
 * @brief the viscosity force induced acceleration, a formulation for conserving
 * angular momentum, to be tested for its practical applications.
 */
class AngularConservativeViscousAccelerationInner : public BaseViscousAccelerationInner
{
  public:
    explicit AngularConservativeViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAccelerationInner(inner_relation){};
    virtual ~AngularConservativeViscousAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class TransportVelocityCorrectionInner
 * @brief transport velocity correction
 */
class TransportVelocityCorrectionInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInner(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &surface_indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
};

/**
 * @class TransportVelocityCorrectionInner
 * @brief transport velocity correction
 */
class TransportVelocityCorrectionInnerAdaptive : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionInnerAdaptive(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &surface_indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
};

template<class FluidT>
class AcousticTimeStepSizeKernel {
  public:
    explicit AcousticTimeStepSizeKernel(BaseParticles* particles) :
        rho_(particles->getDeviceVariableByName<DeviceReal>("Density")),
        p_(particles->registerDeviceVariable<DeviceReal>("Pressure", particles->total_real_particles_)),
        vel_(particles->getDeviceVariableByName<DeviceVecd>("Velocity")),
        fluid_(DynamicCast<FluidT>(this, particles->getBaseMaterial())) {}

    template<class RealType, class VecType, class FluidType, class SoundSpeedFunc, class NormVecdFunc>
    static RealType reduce(size_t index_i, Real dt, FluidType&& fluid, RealType* p, RealType* rho, VecType* vel,
                        SoundSpeedFunc&& getSoundSpeed, NormVecdFunc&& norm) {
        return getSoundSpeed(fluid, p[index_i], rho[index_i]) + norm(vel[index_i]);
    }

    DeviceReal reduce(size_t index_i, Real dt = 0.0) const {
        return reduce(index_i, dt, fluid_, p_, rho_, vel_,
                      [](const FluidT& fluid, DeviceReal p_i, DeviceReal rho_i) {
                          return fluid.getSoundSpeed_Device(p_i, rho_i);
                      },
                      [](const DeviceVecd& vel) {
                          return sycl::length(vel);
                      });
    }

  private:
    DeviceReal *rho_, *p_;
    DeviceVecd *vel_;
    FluidT fluid_;
};

/**
 * @class AcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
{
  public:
    explicit AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

    auto& getDeviceProxy() {
        return device_proxy;
    }

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real acousticCFL_;

  private:
    ExecutionProxy<AcousticTimeStepSize, AcousticTimeStepSizeKernel<WeaklyCompressibleFluid>> device_proxy;
};

using namespace execution;

class AdvectionTimeStepSizeForImplicitViscosityKernel {
  public:
    AdvectionTimeStepSizeForImplicitViscosityKernel(BaseParticles* particles):
        vel_(particles->getDeviceVariableByName<DeviceVecd>("Velocity")) {}

    template<class VecType, class SquareNormFunction>
    static Real reduce(size_t index_i, Real dt, VecType* vel, SquareNormFunction&& squareNorm) {
        return squareNorm(vel[index_i]);
    }

    Real reduce(size_t index_i, Real dt = 0.0) const {
        return reduce(index_i, dt, vel_, [](const DeviceVecd& vel){ return sycl::dot(vel, vel); });
    }

  private:
    DeviceVecd* vel_;
};

/**
 * @class AdvectionTimeStepSizeForImplicitViscosity
 * @brief Computing the advection time step size when viscosity is handled implicitly
 */
class AdvectionTimeStepSizeForImplicitViscosity
    : public LocalDynamicsReduce<Real, ReduceMax>,
      public FluidDataSimple
{
  public:
    explicit AdvectionTimeStepSizeForImplicitViscosity(
        SPHBody &sph_body, Real U_max, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSizeForImplicitViscosity(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real advectionCFL_;

    ExecutionProxy<AdvectionTimeStepSizeForImplicitViscosity,
                   AdvectionTimeStepSizeForImplicitViscosityKernel> device_proxy;

  public:
    auto& getDeviceProxy() {
        return device_proxy;
    }

    using Proxy = decltype(device_proxy);
};

/**
 * @class AdvectionTimeStepSize
 * @brief Computing the advection time step size
 */
class AdvectionTimeStepSize : public AdvectionTimeStepSizeForImplicitViscosity
{
  public:
    explicit AdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
};

/**
 * @class VorticityInner
 * @brief  compute vorticity in the fluid field
 */
class VorticityInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit VorticityInner(BaseInnerRelation &inner_relation);
    virtual ~VorticityInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<AngularVecd> vorticity_;
};

template<class FluidT>
class BaseIntegrationKernel {
  public:
    BaseIntegrationKernel(BaseParticles *particles) :
        fluid_(DynamicCast<FluidT>(this, particles->getBaseMaterial())),
        rho_(particles->getDeviceVariableByName<DeviceReal>("Density")),
        p_(particles->registerDeviceVariable<DeviceReal>("Pressure", particles->total_real_particles_)),
        drho_dt_(particles->registerDeviceVariable<DeviceReal>("DensityChangeRate", particles->total_real_particles_)),
        Vol_(particles->getDeviceVariableByName<DeviceReal>("Volume")),
        mass_(particles->getDeviceVariableByName<DeviceReal>("Mass")),
        pos_(particles->getDeviceVariableByName<DeviceVecd>("Position")),
        vel_(particles->getDeviceVariableByName<DeviceVecd>("Velocity")),
        acc_(particles->getDeviceVariableByName<DeviceVecd>("Acceleration")),
        acc_prior_(particles->getDeviceVariableByName<DeviceVecd>("AccelerationPrior")) {
        particles->registerSortableVariable<Real>("Pressure");
    }
        
  protected:
    FluidT fluid_;
    DeviceReal *rho_, *p_, *drho_dt_, *Vol_, *mass_;
    DeviceVecd *pos_, *vel_, *acc_, *acc_prior_;
};

/**
 * @class BaseIntegration
 * @brief Pure abstract base class for all fluid relaxation schemes
 */
class BaseIntegration : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseIntegration(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
};

template<class RiemannSolverType>
class BaseIntegration1stHalfKernel : public BaseIntegrationKernel<WeaklyCompressibleFluid> {
  public:
    BaseIntegration1stHalfKernel(BaseParticles* particles, NeighborhoodDevice* inner_configuration,
                                 const RiemannSolverType& riemannSolver) :
                                 BaseIntegrationKernel<WeaklyCompressibleFluid>(particles),
                                 riemann_solver_(this->fluid_, this->fluid_),
                                 inner_configuration_(inner_configuration){}

    template<class RealType, class VecType, class FluidType, class PressureFunc>
    static void initialization(size_t index_i, Real dt, RealType* rho, const RealType *drho_dt, RealType *p, VecType* pos,
                               VecType *vel, FluidType&& fluid, PressureFunc&& getPressure) {
        rho[index_i] += drho_dt[index_i] * dt * 0.5;
        p[index_i] = getPressure(fluid, rho[index_i]);
        pos[index_i] += vel[index_i] * dt * 0.5;
    }
    
    template<class Vec>
    static void update(size_t index_i, Real dt, Vec *vel, const Vec *acc_prior, const Vec *acc) {
        vel[index_i] += (acc_prior[index_i] + acc[index_i]) * dt;
    }

    template<class RealType, class VecType, class NeighborhoodType, class RiemannSolver>
    static void interaction(size_t index_i, Real dt, RealType *p, RealType *rho, RealType *drho_dt, VecType* acc,
                            NeighborhoodType* inner_configuration, RiemannSolver&& riemann_solver) {
        auto acceleration = VecdZero<VecType>();
        RealType rho_dissipation(0);
        const auto &inner_neighborhood = inner_configuration[index_i];
        for (size_t n = 0; n < inner_neighborhood.current_size(); ++n)
        {
            const auto& index_j = inner_neighborhood.j_[n];
            const auto& dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
            const auto &e_ij = inner_neighborhood.e_ij_[n];
            
            acceleration -= (p[index_i] + p[index_j]) * dW_ijV_j * e_ij;
            rho_dissipation += riemann_solver.DissipativeUJump(p[index_i] - p[index_j]) * dW_ijV_j;
        }
        acc[index_i] += acceleration / rho[index_i];
        drho_dt[index_i] = rho_dissipation * rho[index_i];
    }
    
    void initialization(size_t index_i, Real dt = 0.0) {
        initialization(index_i, dt, rho_, drho_dt_, p_, pos_, vel_, fluid_, [](const auto& fluid, DeviceReal rho) {
            return fluid.getPressure_Device(rho);
        });
    }

    void update(size_t index_i, Real dt = 0.0) {
        update(index_i, dt, vel_, acc_prior_, acc_);
    }
    
    void interaction(size_t index_i, Real dt = 0.0) {
        interaction(index_i, dt, p_, rho_, drho_dt_, acc_, inner_configuration_, riemann_solver_);
    }
    
  protected:
    RiemannSolverType riemann_solver_;
    NeighborhoodDevice* inner_configuration_;
};

/**
 * @class BaseIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class BaseIntegration1stHalf : public BaseIntegration
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalf(){};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

    using DeviceKernel = BaseIntegration1stHalfKernel<RiemannSolverType>;

  protected:
    virtual Vecd computeNonConservativeAcceleration(size_t index_i);
};
using Integration1stHalf = BaseIntegration1stHalf<NoRiemannSolver>;
/** define the mostly used pressure relaxation scheme using Riemann solver */
using Integration1stHalfRiemann = BaseIntegration1stHalf<AcousticRiemannSolver>;
using Integration1stHalfDissipativeRiemann = BaseIntegration1stHalf<DissipativeRiemannSolver>;


template<class RiemannSolverType>
class BaseIntegration2ndHalfKernel : public BaseIntegrationKernel<WeaklyCompressibleFluid> {
  public:
    BaseIntegration2ndHalfKernel(BaseParticles* particles, NeighborhoodDevice* inner_configuration,
                                 const RiemannSolverType& riemannSolver) :
            BaseIntegrationKernel<WeaklyCompressibleFluid>(particles),
            riemann_solver_(this->fluid_, this->fluid_),
            inner_configuration_(inner_configuration){}

    template<class VecType>
    static inline void initialization(size_t index_i, Real dt, VecType* pos, VecType *vel) {
        pos[index_i] += vel[index_i] * dt * 0.5;
    }

    template<class RealType>
    static void update(size_t index_i, Real dt, RealType* rho, RealType* drho_dt, RealType* Vol, RealType* mass) {
        rho[index_i] += drho_dt[index_i] * dt * 0.5;
        Vol[index_i] = mass[index_i] / rho[index_i];
    }

    template<class RealType, class VecType, class NeighborhoodType, class RiemannSolver, class DotFunc>
    static void interaction(size_t index_i, Real dt, RealType *rho, RealType *drho_dt, VecType* vel, VecType* acc,
                            NeighborhoodType* inner_configuration, RiemannSolver&& riemann_solver, DotFunc&& dot) {
        RealType density_change_rate(0);
        auto p_dissipation = VecdZero<VecType>();
        const auto &inner_neighborhood = inner_configuration[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size(); ++n)
        {
            const auto &index_j = inner_neighborhood.j_[n];
            const auto &e_ij = inner_neighborhood.e_ij_[n];
            const auto &dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

            const RealType u_jump = dot(vel[index_i] - vel[index_j], e_ij);
            density_change_rate += u_jump * dW_ijV_j;
            p_dissipation += static_cast<RealType>(riemann_solver.DissipativePJump(u_jump)) * dW_ijV_j * e_ij;
        }
        drho_dt[index_i] += density_change_rate * rho[index_i];
        acc[index_i] = p_dissipation / rho[index_i];
    }

    void initialization(size_t index_i, Real dt = 0.0) {
        initialization(index_i, dt, pos_, vel_);
    }

    void update(size_t index_i, Real dt = 0.0) {
        update(index_i, dt, rho_, drho_dt_, Vol_, mass_);
    }

    void interaction(size_t index_i, Real dt = 0.0) {
        interaction(index_i, dt, rho_, drho_dt_, vel_, acc_, inner_configuration_, riemann_solver_,
                    [](const DeviceVecd& v1, const DeviceVecd& v2) { return sycl::dot(v1, v2); });
    }

  protected:
    RiemannSolverType riemann_solver_;
    NeighborhoodDevice* inner_configuration_;
};


/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalf : public BaseIntegration
{
  public:
    explicit BaseIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration2ndHalf(){};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

    using DeviceKernel = BaseIntegration2ndHalfKernel<RiemannSolverType>;

  protected:
    StdLargeVec<Real> &Vol_, &mass_;
};
using Integration2ndHalf = BaseIntegration2ndHalf<NoRiemannSolver>;
/** define the mostly used density relaxation scheme using Riemann solver */
using Integration2ndHalfRiemann = BaseIntegration2ndHalf<AcousticRiemannSolver>;
using Integration2ndHalfDissipativeRiemann = BaseIntegration2ndHalf<DissipativeRiemannSolver>;

/**
 * @class Oldroyd_BIntegration1stHalf
 * @brief Pressure relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration1stHalf : public Integration1stHalfDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> tau_, dtau_dt_;
};

/**
 * @class Oldroyd_BIntegration2ndHalf
 * @brief Density relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration2ndHalf : public Integration2ndHalfDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_INNER_H
