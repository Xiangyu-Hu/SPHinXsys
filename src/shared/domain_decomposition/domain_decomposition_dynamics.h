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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    domain_decomposition_dynamics.h
 * @brief   The time loop entries of a decomposed run.
 * @details Each dynamics forwards to one operation of a BodyDecomposition. Under a
 *          non-decomposed execution policy the BodyDecomposition is a no-op, and so
 *          are these, which lets a case place them unconditionally in its loop:
 *
 *            migrate  ->  (sort)  ->  update halo  ->  cell linked list  ->  ...
 *
 *          The interaction algorithms refresh the halo of their interact variables
 *          themselves; SyncHaloStateCK is for variables refreshed once per advection
 *          step by the case, such as the volume or the kernel correction matrix.
 * @author  Niki Loppi, Xiangyu Hu
 */

#ifndef DOMAIN_DECOMPOSITION_DYNAMICS_H
#define DOMAIN_DECOMPOSITION_DYNAMICS_H

#include "base_dynamics.h"
#include "body_decomposition.h"

namespace SPH
{
/**
 * @class UpdateHaloCK
 * @brief Rebuild the halo plan and refresh the full state of the halo particles.
 */
template <class ExecutionPolicy>
class UpdateHaloCK : public BaseDynamics<void>
{
  public:
    explicit UpdateHaloCK(BodyDecomposition<ExecutionPolicy> &decomposition)
        : BaseDynamics<void>(), decomposition_(decomposition) {};
    virtual ~UpdateHaloCK() {};

    virtual void exec(Real dt = 0.0) override { decomposition_.updateHaloPlan(); };

  protected:
    BodyDecomposition<ExecutionPolicy> &decomposition_;
};

/**
 * @class SyncHaloStateCK
 * @brief Refresh a subset of the state variables of the halo particles.
 * @details Sending only what the next interaction reads is what keeps the per-stage
 *          exchange affordable; sending everything would move several times the data.
 */
template <class ExecutionPolicy>
class SyncHaloStateCK : public BaseDynamics<void>
{
  public:
    explicit SyncHaloStateCK(BodyDecomposition<ExecutionPolicy> &decomposition)
        : BaseDynamics<void>(), decomposition_(decomposition),
          particles_(decomposition.getParticles()) {};
    virtual ~SyncHaloStateCK() {};

    /** Usage: sync.addVariable<Vecd>("Velocity").addVariable<Real>("Pressure"); */
    template <typename DataType>
    SyncHaloStateCK &addVariable(const std::string &name)
    {
        particles_.addDiscreteVariableToList<DataType>(variables_, name);
        return *this;
    };

    virtual void exec(Real dt = 0.0) override { decomposition_.refreshHalo(variables_); };

  protected:
    BodyDecomposition<ExecutionPolicy> &decomposition_;
    BaseParticles &particles_;
    DiscreteVariables variables_;
};

/**
 * @class MigrateParticlesCK
 * @brief Transfer ownership of the particles that crossed a cut plane.
 */
template <class ExecutionPolicy>
class MigrateParticlesCK : public BaseDynamics<void>
{
  public:
    explicit MigrateParticlesCK(BodyDecomposition<ExecutionPolicy> &decomposition)
        : BaseDynamics<void>(), decomposition_(decomposition) {};
    virtual ~MigrateParticlesCK() {};

    virtual void exec(Real dt = 0.0) override { decomposition_.migrateParticles(); };

  protected:
    BodyDecomposition<ExecutionPolicy> &decomposition_;
};

/**
 * @class RebalanceSubdomainsCK
 * @brief Move the cut planes towards an equal particle count per subdomain.
 * @details Run rarely, for instance every few hundred advection steps: each move
 *          triggers a migration of everything between the old and the new plane.
 */
template <class ExecutionPolicy>
class RebalanceSubdomainsCK : public BaseDynamics<void>
{
  public:
    explicit RebalanceSubdomainsCK(BodyDecomposition<ExecutionPolicy> &decomposition,
                                   Real relaxation = Real(0.5))
        : BaseDynamics<void>(), decomposition_(decomposition), relaxation_(relaxation) {};
    virtual ~RebalanceSubdomainsCK() {};

    virtual void exec(Real dt = 0.0) override { decomposition_.rebalance(relaxation_); };

  protected:
    BodyDecomposition<ExecutionPolicy> &decomposition_;
    Real relaxation_;
};
} // namespace SPH
#endif // DOMAIN_DECOMPOSITION_DYNAMICS_H
