/**
 * @file    particle_migration_ck.h
 * @brief   Hand particles off between neighbouring MPI ranks when they cross a
 *          subdomain boundary.
 * @details Runs once per step, before HaloExchangeCK::exec() and before the cell
 *          linked list rebuild (see MULTI_GPU_DESIGN.md section 5):
 *
 *              water_update_particle_position.exec();
 *              particle_migration.exec();      // this class
 *              halo_exchange.exec();
 *              water_cell_linked_list.exec();
 *
 *          Unlike the halo exchange, a migrated particle stops existing on the
 *          sending rank and starts existing on the receiving rank, so its *entire*
 *          registered state must move, not just the variables read by neighbour
 *          interactions. The case therefore registers the same or a larger variable
 *          set than it gives to HaloExchangeCK via addToExchange().
 *
 *          IDs are not part of that payload: ParticleOriginalIds()/ParticleSortedIds()
 *          are rank-local bookkeeping (same convention as
 *          BaseParticles::createRealParticleFrom), so the receiving rank assigns a
 *          fresh original id equal to the new slot index rather than copying the
 *          sender's.
 *
 *          Reuses HaloVariablePacker/TypedHaloVariablePacker from halo_exchange_ck.h:
 *          "pack N particles' worth of one variable into a buffer" is exactly the
 *          same operation for a migration payload as for a halo payload.
 * @author  SPHinXsys multi-GPU draft
 */

#ifndef PARTICLE_MIGRATION_CK_H
#define PARTICLE_MIGRATION_CK_H

#include "base_body.h"
#include "base_particles.h"
#include "domain_partition.h"
#include "halo_exchange_ck.h"
#include "mpi_environment.h"
#include "sphinxsys_variable.h"

#include <memory>
#include <vector>

namespace SPH
{
/** MPI message tags used by particle migration. Distinct from the halo exchange tags. */
constexpr int tag_migration_count = 7201;
constexpr int tag_migration_state = 7202;

/**
 * @class ParticleMigrationCK
 * @brief Per-body particle hand-off across the neighbour ranks of a DomainPartition.
 */
template <class ExecutionPolicy>
class ParticleMigrationCK
{
  public:
    /**
     * @param sv_total_owned_particles  The same SingleVariable a paired HaloExchangeCK
     *   reads via svTotalOwnedParticles(). Migration updates it after changing how many
     *   particles this rank owns, so the halo exchange that runs afterwards in the same
     *   step sees the post-migration count rather than a stale one.
     */
    ParticleMigrationCK(RealBody &real_body, const DomainPartition &partition,
                        const MPIEnvironment &mpi_env,
                        SingleVariable<UnsignedInt> *sv_total_owned_particles);
    virtual ~ParticleMigrationCK() {};

    /** Register a variable to be handed off with a migrating particle. Call at setup. */
    template <typename DataType>
    void addToMigrate(const std::string &variable_name);

    /**
     * Move particles that left this rank's subdomain to their new owner, and receive
     * particles handed off by neighbours. Must be called with
     * particles_.TotalRealParticles() == owned particles only, i.e. before any halo
     * particles have been appended for this step.
     */
    void exec(Real dt = 0.0);

    /** Particles received this call. For diagnostics/logging. */
    UnsignedInt LastReceivedCount() const { return last_received_count_; };
    /** Particles sent away this call. For diagnostics/logging. */
    UnsignedInt LastSentCount() const { return last_sent_count_; };

  protected:
    RealBody &real_body_;
    BaseParticles &particles_;
    const DomainPartition &partition_;
    const MPIEnvironment &mpi_env_;

    DiscreteVariable<Vecd> *dv_pos_;
    SingleVariable<UnsignedInt> *sv_total_owned_particles_;

    StdVec<int> neighbour_ranks_;
    /** Per neighbour: indices of owned particles that now belong to that neighbour. */
    StdVec<StdVec<UnsignedInt>> send_indices_;
    StdVec<UnsignedInt> send_counts_;
    StdVec<UnsignedInt> recv_counts_;

    StdVec<std::unique_ptr<HaloVariablePacker>> packers_;
    size_t bytes_per_particle_ = 0;
    StdVec<std::vector<char>> send_buffers_;
    StdVec<std::vector<char>> recv_buffers_;

    UnsignedInt last_sent_count_ = 0;
    UnsignedInt last_received_count_ = 0;

    /** Scan owned particles, sort each one leaving this subdomain into a send list. */
    void buildMigrationLists();
    /** Exchange send_counts_ <-> recv_counts_ with each neighbour. */
    void negotiateCounts();
    /** Pack the registered variables for every particle in the send lists. */
    void packSendBuffers();
    /**
     * Remove particles queued in send_indices_ from the owned range via swap-with-last,
     * processed in descending index order so an already-removed slot's swap partner is
     * never itself pending removal.
     */
    void removeMigratedParticles();
    /** Post the MPI transfer, then append received particles as new real particles. */
    void exchangeAndAppend();
};
} // namespace SPH
#endif // PARTICLE_MIGRATION_CK_H
