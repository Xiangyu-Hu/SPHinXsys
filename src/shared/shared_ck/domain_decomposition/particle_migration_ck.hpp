#ifndef PARTICLE_MIGRATION_CK_HPP
#define PARTICLE_MIGRATION_CK_HPP

#include "particle_migration_ck.h"

#include <algorithm>
#include <iostream>

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
ParticleMigrationCK<ExecutionPolicy>::ParticleMigrationCK(
    RealBody &real_body, const DomainPartition &partition, const MPIEnvironment &mpi_env,
    SingleVariable<UnsignedInt> *sv_total_owned_particles)
    : real_body_(real_body), particles_(real_body.getBaseParticles()),
      partition_(partition), mpi_env_(mpi_env),
      dv_pos_(particles_.getVariableByName<Vecd>("Position")),
      sv_total_owned_particles_(sv_total_owned_particles),
      neighbour_ranks_(partition.NeighbourRanks())
{
    const size_t number_of_neighbours = neighbour_ranks_.size();
    send_indices_.resize(number_of_neighbours);
    send_counts_.resize(number_of_neighbours, 0);
    recv_counts_.resize(number_of_neighbours, 0);
    send_buffers_.resize(number_of_neighbours);
    recv_buffers_.resize(number_of_neighbours);
}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename DataType>
void ParticleMigrationCK<ExecutionPolicy>::addToMigrate(const std::string &variable_name)
{
    DiscreteVariable<DataType> *dv_variable =
        particles_.getVariableByName<DataType>(variable_name);
    if (dv_variable == nullptr)
    {
        std::cout << "\n ERROR: ParticleMigrationCK cannot migrate unknown variable "
                  << variable_name << "!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    packers_.push_back(std::make_unique<TypedHaloVariablePacker<DataType>>(dv_variable));
    bytes_per_particle_ += sizeof(DataType);
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::buildMigrationLists()
{
    for (size_t k = 0; k != send_indices_.size(); ++k)
    {
        send_indices_[k].clear();
    }

    // TODO(device+parallel): serial host scan, same caveat as HaloExchangeCK::buildSendLists.
    Vecd *pos = dv_pos_->Data();
    const UnsignedInt total_owned = particles_.TotalRealParticles();
    for (UnsignedInt index_i = 0; index_i != total_owned; ++index_i)
    {
        const int owner = partition_.ownerRank(pos[index_i]);
        if (owner == partition_.Rank())
        {
            continue;
        }
        if (owner < 0)
        {
            std::cout << "\n ERROR: ParticleMigrationCK: particle " << index_i
                      << " left the system bounds entirely; open boundaries across "
                      << "rank boundaries are not handled." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        bool found = false;
        for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
        {
            if (neighbour_ranks_[k] == owner)
            {
                send_indices_[k].push_back(index_i);
                found = true;
                break;
            }
        }
        if (!found)
        {
            // Owner is a real rank but not one of the ~26 precomputed neighbours: the
            // particle crossed more than one subdomain width this step. That means dt
            // (or the halo width) is too large relative to the subdomain size for this
            // one-hop migration scheme.
            std::cout << "\n ERROR: ParticleMigrationCK: particle " << index_i
                      << " moved into rank " << owner << ", which is not a neighbour of "
                      << "rank " << partition_.Rank() << ". The time step is too large "
                      << "relative to the subdomain size for single-hop migration."
                      << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }

    for (size_t k = 0; k != send_indices_.size(); ++k)
    {
        send_counts_[k] = UnsignedInt(send_indices_[k].size());
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::negotiateCounts()
{
#if SPHINXSYS_USE_MPI
    const int number_of_neighbours = int(neighbour_ranks_.size());
    std::vector<MPI_Request> requests(2 * number_of_neighbours);
    for (int k = 0; k != number_of_neighbours; ++k)
    {
        MPI_Irecv(&recv_counts_[k], sizeof(UnsignedInt), MPI_BYTE, neighbour_ranks_[k],
                  tag_migration_count, MPI_COMM_WORLD, &requests[2 * k]);
        MPI_Isend(&send_counts_[k], sizeof(UnsignedInt), MPI_BYTE, neighbour_ranks_[k],
                  tag_migration_count, MPI_COMM_WORLD, &requests[2 * k + 1]);
    }
    MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
#else
    // Single-rank build: there are no neighbours, so nothing can migrate.
    for (size_t k = 0; k != recv_counts_.size(); ++k)
    {
        recv_counts_[k] = 0;
    }
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::packSendBuffers()
{
    for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
    {
        send_buffers_[k].resize(size_t(send_counts_[k]) * bytes_per_particle_);
        recv_buffers_[k].resize(size_t(recv_counts_[k]) * bytes_per_particle_);

        size_t offset = 0;
        for (size_t v = 0; v != packers_.size(); ++v)
        {
            packers_[v]->pack(send_indices_[k].data(), send_counts_[k],
                              send_buffers_[k].data() + offset);
            offset += size_t(send_counts_[k]) * packers_[v]->bytesPerParticle();
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::removeMigratedParticles()
{
    StdVec<UnsignedInt> departing;
    for (size_t k = 0; k != send_indices_.size(); ++k)
    {
        departing.insert(departing.end(), send_indices_[k].begin(), send_indices_[k].end());
    }
    last_sent_count_ = UnsignedInt(departing.size());

    // A particle is sorted into at most one neighbour's send list in buildMigrationLists,
    // so indices here are already unique; only the order needs fixing. Descending order
    // means every swap-with-last pulls in a particle that was never itself queued for
    // removal, because all indices greater than the current one have already been
    // removed (see particle_migration_ck.h for the argument).
    std::sort(departing.begin(), departing.end(), std::greater<UnsignedInt>());
    for (UnsignedInt index : departing)
    {
        particles_.switchToBufferParticle(index);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::exchangeAndAppend()
{
#if SPHINXSYS_USE_MPI
    const int number_of_neighbours = int(neighbour_ranks_.size());
    std::vector<MPI_Request> requests(2 * number_of_neighbours);
    for (int k = 0; k != number_of_neighbours; ++k)
    {
        MPI_Irecv(recv_buffers_[k].data(), int(recv_buffers_[k].size()), MPI_BYTE,
                  neighbour_ranks_[k], tag_migration_state, MPI_COMM_WORLD, &requests[2 * k]);
        MPI_Isend(send_buffers_[k].data(), int(send_buffers_[k].size()), MPI_BYTE,
                  neighbour_ranks_[k], tag_migration_state, MPI_COMM_WORLD, &requests[2 * k + 1]);
    }
    MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
#endif // SPHINXSYS_USE_MPI

    last_received_count_ = 0;
    for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
    {
        if (recv_counts_[k] == 0)
        {
            continue;
        }

        if (particles_.TotalRealParticles() + recv_counts_[k] > particles_.ParticlesBound())
        {
            std::cout << "\n ERROR: ParticleMigrationCK: rank " << mpi_env_.Rank()
                      << " received " << recv_counts_[k] << " particles from rank "
                      << neighbour_ranks_[k] << " but only "
                      << (particles_.ParticlesBound() - particles_.TotalRealParticles())
                      << " buffer slots remain. Increase the particle buffer reserve."
                      << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        const UnsignedInt first_slot = particles_.TotalRealParticles();
        for (UnsignedInt p = 0; p != recv_counts_[k]; ++p)
        {
            // Allocates a slot with correctly-assigned original/sorted id (same
            // bookkeeping as any other new real particle). The state copied from
            // particle 0 is a throwaway placeholder, overwritten by unpack() below.
            particles_.createRealParticleFrom(0);
        }

        size_t offset = 0;
        for (size_t v = 0; v != packers_.size(); ++v)
        {
            packers_[v]->unpack(first_slot, recv_counts_[k],
                                recv_buffers_[k].data() + offset);
            offset += size_t(recv_counts_[k]) * packers_[v]->bytesPerParticle();
        }

        last_received_count_ += recv_counts_[k];
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void ParticleMigrationCK<ExecutionPolicy>::exec(Real dt)
{
    (void)dt;
    if (packers_.empty())
    {
        std::cout << "\n ERROR: ParticleMigrationCK has no registered variables; "
                  << "call addToMigrate() at setup!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    buildMigrationLists();
    negotiateCounts();
    packSendBuffers();
    // Data is captured in the send buffers before any indices move, so removal can now
    // freely swap particles around.
    removeMigratedParticles();
    exchangeAndAppend();

    // HaloExchangeCK::buildSendLists() reads this at the start of its next exec(), so it
    // must reflect the post-migration owned count before that call happens.
    sv_total_owned_particles_->setValue(particles_.TotalRealParticles());
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_MIGRATION_CK_HPP
