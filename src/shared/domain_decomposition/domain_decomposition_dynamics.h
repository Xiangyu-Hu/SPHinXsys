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
 * @file 	domain_decomposition_dynamics.h
 * @brief 	Time loop entries of a domain decomposed run.
 * @details These are ordinary BaseDynamics objects, so a case file inserts them into
 *          the time loop next to the cell linked list and relation updates. The order
 *          within an advection step is:
 *
 *              update particle positions
 *              MigrateParticlesCK          ownership follows the positions
 *              particle sort               optional, per subdomain
 *              UpdateHaloCK                new plan, and a full state refresh
 *              update cell linked list     over owned plus halo
 *              update body relations       neighbor lists of the owned particles
 *
 *          The halo must be rebuilt before the cell linked list: the list is built
 *          over the local particles, and it is the halo plan that publishes how many
 *          of those there are.
 *
 *          and within an acoustic step, after every stage that writes a state variable
 *          read by the next interaction:
 *
 *              SyncHaloStateCK({"Velocity", "Pressure", ...})
 *
 *          The last point is the one that costs: a halo refresh per acoustic stage.
 *          It is the price of a halo one cut-off deep. A deeper halo would let several
 *          stages run between exchanges, at the cost of redundant computation on the
 *          halo particles; which of the two wins depends on the particle count per
 *          device and on the interconnect, so both should be measured.
 * @author	Niki Loppi
 */

#ifndef DOMAIN_DECOMPOSITION_DYNAMICS_H
#define DOMAIN_DECOMPOSITION_DYNAMICS_H

#include "base_dynamics.h"
#include "subdomain_exchange.hpp"

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
    explicit UpdateHaloCK(SubdomainExchange<ExecutionPolicy> &exchange)
        : BaseDynamics<void>(), exchange_(exchange) {};
    virtual ~UpdateHaloCK() {};

    virtual void exec(Real dt = 0.0) override { exchange_.updateHaloPlan(); };

  protected:
    SubdomainExchange<ExecutionPolicy> &exchange_;
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
    SyncHaloStateCK(SubdomainExchange<ExecutionPolicy> &exchange, BaseParticles &particles)
        : BaseDynamics<void>(), exchange_(exchange), particles_(particles) {};
    virtual ~SyncHaloStateCK() {};

    /** Usage: sync.addVariable<Vecd>("Velocity").addVariable<Real>("Pressure"); */
    template <typename DataType>
    SyncHaloStateCK &addVariable(const std::string &name)
    {
        particles_.addDiscreteVariableToList<DataType>(variables_, name);
        return *this;
    };

    virtual void exec(Real dt = 0.0) override { exchange_.refreshHalo(variables_); };

  protected:
    SubdomainExchange<ExecutionPolicy> &exchange_;
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
    explicit MigrateParticlesCK(SubdomainExchange<ExecutionPolicy> &exchange)
        : BaseDynamics<void>(), exchange_(exchange) {};
    virtual ~MigrateParticlesCK() {};

    virtual void exec(Real dt = 0.0) override { exchange_.migrateParticles(); };

  protected:
    SubdomainExchange<ExecutionPolicy> &exchange_;
};

/**
 * @class RebalanceSubdomainsCK
 * @brief Move the cut planes towards an equal particle count per device.
 * @details Run rarely, for instance every few hundred advection steps: each move
 *          triggers a migration of everything between the old and the new plane.
 */
template <class ExecutionPolicy>
class RebalanceSubdomainsCK : public BaseDynamics<void>
{
  public:
    RebalanceSubdomainsCK(SlabDecomposition &decomposition, SubdomainExchange<ExecutionPolicy> &exchange,
                          Real relaxation = Real(0.5))
        : BaseDynamics<void>(), decomposition_(decomposition),
          exchange_(exchange), relaxation_(relaxation) {};
    virtual ~RebalanceSubdomainsCK() {};

    virtual void exec(Real dt = 0.0) override
    {
        if (decomposition_.rebalance(exchange_.OwnedParticlesPerSubdomain(), relaxation_))
        {
            exchange_.migrateParticles();
        }
    };

  protected:
    SlabDecomposition &decomposition_;
    SubdomainExchange<ExecutionPolicy> &exchange_;
    Real relaxation_;
};
} // namespace SPH
#endif // DOMAIN_DECOMPOSITION_DYNAMICS_H
