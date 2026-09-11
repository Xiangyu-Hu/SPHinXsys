/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-sph_body dynamics and beyond with SPH   *
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
 * @file    body_decomposition.h
 * @brief   The decomposition of one sph_body over the subdomains of a run, selected by
 *          the execution policy.
 * @details A case file drives the decomposition through one object of this class
 *          and a few dynamics built on it (domain_decomposition_dynamics.h). The
 *          primary template is the non-decomposed case: every operation is a no-op
 *          and the sph_body keeps its single global particle set. The partial
 *          specialization on DecomposedExecution<> owns the cut planes and the halo
 *          and migration exchange. Since MainExecutionPolicy is a decomposed policy
 *          exactly when SPHINXSYS_DECOMPOSITION is on, a case written against this
 *          class runs decomposed or not without any conditional compilation.
 *
 *          The exchange set, that is the variables a particle carries when it
 *          migrates and the variables staged to the host for output, is fixed when
 *          the object is constructed: the evolving variables of the sph_body and the
 *          variables registered for output at that time, plus whatever the case adds
 *          with addExchangeVariable(). Construct it after every dynamics and after
 *          the output variables are registered.
 * @author  Niki Loppi, Xiangyu Hu
 */

#ifndef BODY_DECOMPOSITION_H
#define BODY_DECOMPOSITION_H

#include "adaptation.h"
#include "base_body.h"
#include "base_particles.hpp"
#include "domain_decomposition.h"
#include "execution_policy.h"
#include "sph_system.h"
#include "subdomain_exchange.hpp"
#include "subdomain_runner.h"

#include <memory>
#include <string>

namespace SPH
{
class BaseDecomposition
{
  public:
    virtual ~BaseDecomposition() {};
};

/**
 * @class BodyDecomposition
 * @brief Non-decomposed run: the single global particle set is the only "subdomain".
 * @details The recorders gather a decomposed sph_body to the host themselves, through the
 *          SubdomainExchangeInterface the decomposed specialization installs on the
 *          particles; a case therefore only places scatterFromHost() and the loop
 *          dynamics. The gather and finish methods stay public for host side work a
 *          case does on its own.
 */
template <class ExecutionPolicy>
class BodyDecomposition : public BaseDecomposition
{
  public:
    explicit BodyDecomposition(SPHBody &sph_body, int split_axis = -1)
        : body_(sph_body), particles_(sph_body.getBaseParticles()) {};
    virtual ~BodyDecomposition() {};

    BaseParticles &getParticles() { return particles_; };
    template <typename DataType>
    BodyDecomposition &addExchangeVariable(const std::string &name) { return *this; };
    template <class RecorderType>
    void addSubdomainIDToWrite(RecorderType &recorder) {};

    void scatterFromHost() {};
    void gatherToHost() {};
    void finishHostAccess() {};
    void updateHaloPlan() {};
    void refreshHalo(DiscreteVariables &variables) {};
    void migrateParticles() {};
    bool rebalance(Real relaxation) { return false; };

    int NumberOfSubdomains() const { return 1; };
    StdVec<UnsignedInt> OwnedParticlesPerSubdomain() { return {particles_.TotalRealParticles()}; };
    UnsignedInt TotalOwnedParticles() { return particles_.TotalRealParticles(); };
    Real HaloLoadFactor() const { return Real(0); };
    std::string checkConsistency() const { return std::string(); };
    std::string describe() const { return "Body " + body_.Name() + " is not decomposed\n"; };

  protected:
    SPHBody &body_;
    BaseParticles &particles_;
};

/** Append every variable of one data assemble to an exchange set, skipping duplicates. */
struct AddVariablesToExchangeSet
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    BaseParticles &particles, DiscreteVariables &exchange_set)
    {
        for (std::size_t k = 0; k != variables.size(); ++k)
        {
            particles.template addDiscreteVariableToList<DataType>(exchange_set, variables[k]);
        }
    }
};

/**
 * @class BodyDecomposition<DecomposedExecution<PolicyType>>
 * @brief Decomposed run: slab cut planes plus the halo and migration exchange of the sph_body.
 */
template <class PolicyType>
class BodyDecomposition<DecomposedExecution<PolicyType>>
    : public BaseDecomposition, public SubdomainExchangeInterface
{
    using ExecutionPolicy = DecomposedExecution<PolicyType>;
    std::unique_ptr<SubdomainExchange<ExecutionPolicy>> exchange_;

  public:
    /**
     * @param sph_body        the sph_body to decompose; its particles must be generated
     * @param split_axis  axis to cut along; by default the longest one
     *
     * The cut planes are balanced on the initial particle positions, so that a sph_body
     * occupying only part of the domain does not leave a subdomain empty.
     */
    explicit BodyDecomposition(SPHBody &sph_body, int split_axis = -1)
        : body_(sph_body), particles_(sph_body.getBaseParticles()),
          decomposition_(sph_body.getSPHSystem().getSystemDomainBounds(),
                         sph_body.getSPHAdaptation().getKernel()->CutOffRadius(),
                         execution::subdomain_runner.NumberOfSubdomains(), split_axis),
          dv_subdomain_id_(nullptr)
    {
        balanceOnInitialPositions();
        std::cout << decomposition_.describe();
        particles_.setSubdomainExchange(this);
    };
    virtual ~BodyDecomposition() {};

    BaseParticles &getParticles() { return particles_; };
    SlabDecomposition &getSlabDecomposition() { return decomposition_; };
    SubdomainExchange<ExecutionPolicy> &getExchange()
    {
        if (!exchange_)
        {
            OperationOnDataAssemble<DiscreteVariables, AddVariablesToExchangeSet> add_variables;
            add_variables(particles_.EvolvingVariables(), particles_, exchange_variables_);
            add_variables(particles_.VariablesToWrite(), particles_, exchange_variables_);
            exchange_ = std::make_unique<SubdomainExchange<ExecutionPolicy>>(
                decomposition_, particles_, exchange_variables_);
        }
        return *exchange_;
    };

    /** Add a variable to the exchange set, for instance one that is read at the
     *  neighbors by an interaction but is neither evolving nor written out. Must be
     *  called before scatterFromHost(). */
    template <typename DataType>
    BodyDecomposition &addExchangeVariable(const std::string &name)
    {
        auto &exchange = getExchange();
        exchange.template addExchangeVariable<DataType>(
            particles_.template getVariableByName<DataType>(name));
        return *this;
    };

    /** Write the owning subdomain of every particle as the variable "SubdomainID".
     *  The tag is filled on the host by gatherToHost(), which lays the owned
     *  particles out subdomain by subdomain, so it costs no exchange. */
    template <class RecorderType>
    void addSubdomainIDToWrite(RecorderType &recorder)
    {
        dv_subdomain_id_ = particles_.template registerStateVariable<int>("SubdomainID");
        recorder.template addToWrite<int>(body_, "SubdomainID");
    };

    void scatterFromHost()
    {
        auto &exchange = getExchange();
        exchange.scatterFromHost();
    };

    virtual void gatherToHost() override
    {
        auto &exchange = getExchange();
        exchange.gatherToHost();

        if (dv_subdomain_id_ != nullptr)
        {
            int *subdomain_id = dv_subdomain_id_->Data();
            UnsignedInt offset = 0;
            const StdVec<UnsignedInt> owned = exchange.OwnedParticlesPerSubdomain();
            for (int s = 0; s < decomposition_.NumberOfSubdomains(); ++s)
            {
                for (UnsignedInt i = offset; i < offset + owned[s]; ++i)
                    subdomain_id[i] = s;
                offset += owned[s];
            }
        }
    };

    virtual void finishHostAccess() override
    {
        auto &exchange = getExchange();
        exchange.finishHostAccess();
    };

    virtual int subdomainOf(const Vecd &position) const override
    {
        return decomposition_.getSubdomainMap().subdomainOf(position);
    };

    void updateHaloPlan()
    {
        auto &exchange = getExchange();
        exchange.updateHaloPlan();
    };

    virtual void refreshHalo(DiscreteVariables &variables) override
    {
        auto &exchange = getExchange();
        exchange.refreshHalo(variables);
    };

    void migrateParticles()
    {
        auto &exchange = getExchange();
        exchange.migrateParticles();
    };

    /** Move the cut planes towards equal owned counts; migrates when a plane moved. */
    bool rebalance(Real relaxation)
    {
        auto &exchange = getExchange();
        if (decomposition_.rebalance(exchange.OwnedParticlesPerSubdomain(), relaxation))
        {
            exchange.migrateParticles();
            return true;
        }
        return false;
    };

    int NumberOfSubdomains() const { return decomposition_.NumberOfSubdomains(); };
    StdVec<UnsignedInt> OwnedParticlesPerSubdomain()
    {
        auto &exchange = getExchange();
        return exchange.OwnedParticlesPerSubdomain();
    };

    UnsignedInt TotalOwnedParticles()
    {
        auto &exchange = getExchange();
        return exchange.TotalOwnedParticles();
    };

    Real HaloLoadFactor() const
    {
        auto &exchange = getExchange();
        return exchange.HaloLoadFactor();
    };

    std::string checkConsistency() const
    {
        auto &exchange = getExchange();
        return exchange.checkConsistency();
    };

    std::string describe() const { return decomposition_.describe(); };

  protected:
    /** Iterate the one dimensional rebalance on the initial positions until the cut
     *  planes settle. A full correction per iteration is fine here: no migration is
     *  triggered before the scatter. */
    void balanceOnInitialPositions()
    {
        const int number_of_subdomains = decomposition_.NumberOfSubdomains();
        Vecd *position = particles_.dvParticlePosition()->Data();
        const UnsignedInt total_particles = particles_.TotalRealParticles();
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            StdVec<UnsignedInt> counts(number_of_subdomains, 0);
            for (UnsignedInt i = 0; i < total_particles; ++i)
                counts[decomposition_.getSubdomainMap().subdomainOf(position[i])]++;
            if (!decomposition_.rebalance(counts, Real(1)))
                break;
        }
    };

    SPHBody &body_;
    BaseParticles &particles_;
    SlabDecomposition decomposition_;
    DiscreteVariables exchange_variables_;
    DiscreteVariable<int> *dv_subdomain_id_;
};
} // namespace SPH
#endif // BODY_DECOMPOSITION_H
