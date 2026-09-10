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
 * @file 	subdomain_exchange.h
 * @brief 	Halo exchange and particle migration between subdomains.
 * @details Policy generic: the same code drives a multi-GPU run and a host side
 *          decomposition. Everything it needs is already policy generic --
 *          particle_for over an IndexRange, exclusive_scan, DelegatedData -- except
 *          the copy between two subdomains' replicas, which is
 *          execution::copyBetweenSubdomains and is the only backend specific call here.
 *
 *          Layout of the particle arrays of one subdomain:
 *
 *              [0, n_owned)                    particles owned by this subdomain
 *              [n_owned, n_local)              halo copies received from neighbors
 *              [n_local, particles_bound)      spare capacity
 *
 *          Physics is integrated over [0, n_owned): a halo particle is a read-only copy
 *          whose value is authored by its owner. The cell linked list, and hence the
 *          neighbor search, covers [0, n_local), so a particle next to a cut plane
 *          keeps a complete support.
 *
 *          Two operations maintain that layout:
 *          - the halo plan, rebuilt whenever particles have moved between cells,
 *            determines which owned particles fall into a neighbor's halo band and
 *            where the incoming copies land;
 *          - the halo refresh, run after every step that writes a state variable used
 *            by the next interaction, re-sends the values along the existing plan.
 *          Splitting the two matters: the plan involves a scan and a reallocation
 *          check, the refresh is a pack plus a copy.
 *
 *          Migration transfers ownership of the particles that crossed a cut plane.
 *          It is only needed at the advection step, where positions change by more
 *          than a fraction of a cell.
 *
 *          Ordering, per exchange, with two barriers over the subdomains:
 *            1. every subdomain packs into its own send buffers;
 *            2. barrier;
 *            3. every subdomain pulls from its neighbors' send buffers into its own
 *               arrays;
 *            4. barrier.
 *          The pull direction means a subdomain only ever writes memory it owns, so the
 *          exchange needs no cross-subdomain atomics or locks: the barriers alone order
 *          it. Under the sequential host runner the barriers are trivially satisfied,
 *          which is the configuration to debug the logic in; under the threaded runner
 *          and on the device they are real, which is the configuration to debug the
 *          synchronization in.
 * @author	Niki Loppi
 */

#ifndef SUBDOMAIN_EXCHANGE_H
#define SUBDOMAIN_EXCHANGE_H

#include "base_particles.h"
#include "domain_decomposition.h"
#include "sphinxsys_variable.h"
#include "subdomain_fan_out.h"

namespace SPH
{
/** Lower and upper neighbor of a slab. */
constexpr int NumberOfSides = 2;

/**
 * @class VariableExchangeBuffer
 * @brief Packing buffers of one particle variable, replicated on every subdomain.
 * @details The buffers are plain DiscreteVariables, which gives the per-subdomain
 *          replication for free: a DelegatedData() call inside SubdomainScope(d) yields
 *          the buffer belonging to subdomain d.
 */
template <class ExecutionPolicy, typename DataType>
class VariableExchangeBuffer
{
    UniquePtrsKeeper<DiscreteVariable<DataType>> buffer_keeper_;

  public:
    VariableExchangeBuffer(DiscreteVariable<DataType> *variable, UnsignedInt capacity);

    DiscreteVariable<DataType> *Variable() { return variable_; };
    /** Send buffer of `subdomain_id` for the given side. Callable from another
     *  subdomain's thread; the replicas are pre-touched at construction so that no
     *  allocation can race here. */
    DataType *SendBuffer(int subdomain_id, int side);
    void reserve(UnsignedInt capacity);
    UnsignedInt Capacity() const { return capacity_; };

  protected:
    DiscreteVariable<DataType> *variable_;
    DiscreteVariable<DataType> *send_buffer_[NumberOfSides];
    UnsignedInt capacity_;
};

/**
 * @class SubdomainExchange
 * @brief Halo and migration exchange of the particles of one body.
 */
template <class ExecutionPolicy>
class SubdomainExchange
{
    /** The buffer assembles are keyed on the policy too, so that a translation unit
     *  which instantiates both policies keeps their buffers apart. */
    template <typename DataType>
    using BufferType = VariableExchangeBuffer<ExecutionPolicy, DataType>;

  public:
    using ExchangeBuffers = DataContainerAddressAssemble<BufferType>;
    using ExchangeBufferPtrs = DataContainerUniquePtrAssemble<BufferType>;

    SubdomainExchange(SlabDecomposition &decomposition, BaseParticles &particles,
                      DiscreteVariables &variables_to_exchange);
    virtual ~SubdomainExchange() {};

    /** Add a variable to the exchange set, with its send buffers. Skips a variable
     *  already in the set. Must precede scatterFromHost(), which permutes the host
     *  arrays of the set only. */
    template <typename DataType>
    void addExchangeVariable(DiscreteVariable<DataType> *variable);

    /** Distribute an initially global particle set over the subdomains. Called once,
     *  after particle generation, from the host thread. */
    void scatterFromHost();
    /** Assemble the owned particles of all subdomains back into the host arrays, in
     *  subdomain order. Used for output, restart and host side dynamics. While the
     *  host arrays are in use, the particle counter seen from the host thread is the
     *  global one; call finishHostAccess() afterwards. */
    void gatherToHost();
    /** Undo the counter publication of gatherToHost(), so that subdomain 0 sees its
     *  own owned count again. Must be called before the next fan-out. */
    void finishHostAccess();

    /** Recompute which particles are exchanged and how many arrive, and publish the
     *  resulting n_local. Must run after every migration or sort, and before the cell
     *  linked list is rebuilt, since that list covers the local particles. */
    void updateHaloPlan();
    /** Re-send the values of the given variables along the current plan. Also reached
     *  through BaseParticles::refreshHalo(), which the interaction algorithms call with
     *  their interact variables before every interaction step. */
    void refreshHalo(DiscreteVariables &variables);
    void refreshHalo() { refreshHalo(variables_to_exchange_); };
    /** Transfer ownership of the particles that crossed a cut plane.
     *
     *  The departing particles are packed into the send buffers first, then removed
     *  by filling their slots from the tail: the k-th departing slot below the new
     *  owned count receives the k-th staying particle at or above it, so that the
     *  staying particles end up contiguous in [0, n_new) with O(departing) copies.
     *  This is the "swap with the last real particle" of particle deletion, but
     *  driven by the same flag-and-scan lists as the packing, which keeps it free of
     *  atomics and deterministic. The arrivals are then appended after n_new, which is
     *  the particle generation side: their state comes from the neighbor's buffer. */
    void migrateParticles();

    UnsignedInt OwnedParticles(int subdomain_id) const { return owned_count_[subdomain_id]; };
    StdVec<UnsignedInt> OwnedParticlesPerSubdomain() const { return owned_count_; };
    UnsignedInt TotalOwnedParticles() const;
    /** Largest fill ratio of any send buffer, for diagnosing a halo that outgrew the
     *  reserved capacity. */
    Real HaloLoadFactor() const;

    /** Invariant check, cheap enough to leave on in a debug run: the owned counts sum
     *  to the global particle number, every owned particle really is inside its slab,
     *  and no subdomain exceeded its particle bound. Returns an empty string when the
     *  state is consistent, a description of the first violation otherwise. This is
     *  the main reason the host path exists. */
    std::string checkConsistency() const;

  protected:
    /** Steps of an exchange, each executed inside the subdomain fan-out. */
    void buildSendListsOnCurrentSubdomain();
    void publishLocalCountOnCurrentSubdomain();
    void packOnCurrentSubdomain(DiscreteVariables &variables);
    void pullFromNeighborsOnCurrentSubdomain(DiscreteVariables &variables);
    void reserveBuffers(UnsignedInt capacity);
    UnsignedInt largestSendCount() const;

    SlabDecomposition &decomposition_;
    BaseParticles &particles_;
    DiscreteVariables &variables_to_exchange_;
    ExchangeBufferPtrs exchange_buffer_ptrs_;
    ExchangeBuffers exchange_buffers_;

    /** Per-subdomain work arrays; replicated by DiscreteVariable itself. */
    DiscreteVariable<UnsignedInt> *dv_send_flag_;
    DiscreteVariable<UnsignedInt> *dv_send_scan_;
    DiscreteVariable<UnsignedInt> *dv_send_index_[NumberOfSides];
    /** Migration: slots of departing particles below the new owned count, and the
     *  staying particles at or above it which move into those slots. */
    DiscreteVariable<UnsignedInt> *dv_hole_index_;
    DiscreteVariable<UnsignedInt> *dv_donor_index_;

    /** Host side bookkeeping, written by the subdomain threads at disjoint indices. */
    StdVec<std::array<UnsignedInt, NumberOfSides>> send_count_;
    StdVec<std::array<UnsignedInt, NumberOfSides>> recv_count_;
    StdVec<std::array<UnsignedInt, NumberOfSides>> halo_offset_;
    StdVec<UnsignedInt> owned_count_;
    StdVec<UnsignedInt> keep_count_; /**< owned count after the departures, before arrivals */
    StdVec<UnsignedInt> fill_count_; /**< number of holes filled from the tail */
    UnsignedInt buffer_capacity_;
};
} // namespace SPH
#endif // SUBDOMAIN_EXCHANGE_H
