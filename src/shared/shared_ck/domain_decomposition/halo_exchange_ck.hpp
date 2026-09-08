#ifndef HALO_EXCHANGE_CK_HPP
#define HALO_EXCHANGE_CK_HPP

#include "halo_exchange_ck.h"

#include <iostream>

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void TypedHaloVariablePacker<DataType>::pack(const UnsignedInt *indices, UnsignedInt count, void *buffer)
{
    // TODO(device): Data() is the host allocation. On ParallelDevicePolicy the authoritative
    // copy lives in device memory, so this needs either
    //   (a) a device-side gather kernel into a device staging buffer + GPU-aware MPI, or
    //   (b) dv_variable_->synchronizeWithDevice() before packing (correct but a full-array
    //       device-to-host copy every step, which defeats the purpose).
    // Option (a) is the one worth having; this host path exists so the CPU policy works first.
    DataType *data = dv_variable_->Data();
    DataType *out = reinterpret_cast<DataType *>(buffer);
    for (UnsignedInt i = 0; i != count; ++i)
    {
        out[i] = data[indices[i]];
    }
}
//=================================================================================================//
template <typename DataType>
void TypedHaloVariablePacker<DataType>::unpack(UnsignedInt first_slot, UnsignedInt count, const void *buffer)
{
    DataType *data = dv_variable_->Data();
    const DataType *in = reinterpret_cast<const DataType *>(buffer);
    for (UnsignedInt i = 0; i != count; ++i)
    {
        data[first_slot + i] = in[i];
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
HaloExchangeCK<ExecutionPolicy>::HaloExchangeCK(
    RealBody &real_body, const DomainPartition &partition,
    const MPIEnvironment &mpi_env, Real halo_reserve_factor)
    : real_body_(real_body), particles_(real_body.getBaseParticles()),
      partition_(partition), mpi_env_(mpi_env),
      dv_pos_(particles_.getVariableByName<Vecd>("Position")),
      sv_total_owned_particles_(
          particles_.registerSingleVariable<UnsignedInt>("TotalOwnedParticles")),
      neighbour_ranks_(partition.NeighbourRanks())
{
    total_owned_particles_ = particles_.TotalRealParticles();
    sv_total_owned_particles_->setValue(total_owned_particles_);

    // Halo particles are written into the buffer region directly above the owned particles,
    // so that [0, n_owned + n_halo) is contiguous and the cell linked list can simply treat
    // the whole span as real. The case must therefore reserve buffer particles, e.g.
    //   body.generateParticles<BaseParticles, Lattice>(ParticleBuffer<ReserveSizeFactor>(0.5));
    halo_capacity_ = UnsignedInt(halo_reserve_factor * Real(total_owned_particles_));
    if (particles_.ParticlesBound() < total_owned_particles_ + halo_capacity_)
    {
        std::cout << "\n ERROR: HaloExchangeCK needs " << halo_capacity_
                  << " buffer slots above " << total_owned_particles_
                  << " owned particles, but ParticlesBound() is "
                  << particles_.ParticlesBound() << "!" << std::endl;
        std::cout << "\n Increase the particle buffer reserve for this body." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    const size_t number_of_neighbours = neighbour_ranks_.size();
    send_indices_.resize(number_of_neighbours);
    send_counts_.resize(number_of_neighbours, 0);
    recv_counts_.resize(number_of_neighbours, 0);
    recv_first_slot_.resize(number_of_neighbours, 0);
    send_buffers_.resize(number_of_neighbours);
    recv_buffers_.resize(number_of_neighbours);
}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename DataType>
void HaloExchangeCK<ExecutionPolicy>::addToExchange(const std::string &variable_name)
{
    DiscreteVariable<DataType> *dv_variable =
        particles_.getVariableByName<DataType>(variable_name);
    if (dv_variable == nullptr)
    {
        std::cout << "\n ERROR: HaloExchangeCK cannot exchange unknown variable "
                  << variable_name << "!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    packers_.push_back(std::make_unique<TypedHaloVariablePacker<DataType>>(dv_variable));
    bytes_per_particle_ += sizeof(DataType);
}
//=================================================================================================//
template <class ExecutionPolicy>
void HaloExchangeCK<ExecutionPolicy>::buildSendLists()
{
    // Owned count is authoritative here: a preceding migration step may have changed it.
    total_owned_particles_ = sv_total_owned_particles_->getValue();

    for (size_t k = 0; k != send_indices_.size(); ++k)
    {
        send_indices_[k].clear();
    }

    // TODO(device+parallel): serial host scan over all owned particles. The device version
    // should mirror UpdateCellLinkedList: a per-neighbour atomic counter pass, an
    // exclusive_scan, then a scatter pass building the index lists on device. Restricting
    // the scan to boundary cells of the cell linked list would cut it further, which is
    // what PeriodicConditionUsingGhostParticles does via bound_cells_data_.
    Vecd *pos = dv_pos_->Data();
    StdVec<int> targets;
    for (UnsignedInt index_i = 0; index_i != total_owned_particles_; ++index_i)
    {
        targets.clear();
        partition_.haloTargetRanks(pos[index_i], targets);
        for (size_t t = 0; t != targets.size(); ++t)
        {
            for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
            {
                if (neighbour_ranks_[k] == targets[t])
                {
                    send_indices_[k].push_back(index_i);
                    break;
                }
            }
        }
    }

    for (size_t k = 0; k != send_indices_.size(); ++k)
    {
        send_counts_[k] = UnsignedInt(send_indices_[k].size());
    }

    negotiateCounts();
}
//=================================================================================================//
template <class ExecutionPolicy>
void HaloExchangeCK<ExecutionPolicy>::negotiateCounts()
{
#if SPHINXSYS_USE_MPI
    const int number_of_neighbours = int(neighbour_ranks_.size());
    std::vector<MPI_Request> requests(2 * number_of_neighbours);

    for (int k = 0; k != number_of_neighbours; ++k)
    {
        // UnsignedInt is uint32_t or size_t depending on the build, so ship raw bytes
        // rather than commit to an MPI integer type here.
        MPI_Irecv(&recv_counts_[k], sizeof(UnsignedInt), MPI_BYTE, neighbour_ranks_[k],
                  tag_halo_count, MPI_COMM_WORLD, &requests[2 * k]);
        MPI_Isend(&send_counts_[k], sizeof(UnsignedInt), MPI_BYTE, neighbour_ranks_[k],
                  tag_halo_count, MPI_COMM_WORLD, &requests[2 * k + 1]);
    }
    MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
#endif // SPHINXSYS_USE_MPI

    // Lay the incoming particles out contiguously above the owned particles.
    halo_first_slot_ = total_owned_particles_;
    total_halo_particles_ = 0;
    for (size_t k = 0; k != recv_counts_.size(); ++k)
    {
        recv_first_slot_[k] = halo_first_slot_ + total_halo_particles_;
        total_halo_particles_ += recv_counts_[k];
    }
    checkHaloCapacity();

    // Make the halo visible to UpdateCellLinkedList and the relation update, which both
    // loop over [0, TotalRealParticles()).
    particles_.svTotalRealParticles()->setValue(total_owned_particles_ + total_halo_particles_);
}
//=================================================================================================//
template <class ExecutionPolicy>
void HaloExchangeCK<ExecutionPolicy>::checkHaloCapacity() const
{
    if (total_halo_particles_ > halo_capacity_)
    {
        std::cout << "\n ERROR: rank " << mpi_env_.Rank() << " received "
                  << total_halo_particles_ << " halo particles but only "
                  << halo_capacity_ << " slots are reserved!" << std::endl;
        std::cout << "\n Increase halo_reserve_factor and the particle buffer reserve." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void HaloExchangeCK<ExecutionPolicy>::exchangeState()
{
    if (packers_.empty())
    {
        std::cout << "\n ERROR: HaloExchangeCK has no registered variables; "
                  << "call addToExchange() at setup!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    // Buffers are laid out variable-major: all of variable 0 for every particle, then
    // variable 1, and so on. Each packer therefore writes one contiguous block.
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

#if SPHINXSYS_USE_MPI
    const int number_of_neighbours = int(neighbour_ranks_.size());
    std::vector<MPI_Request> requests(2 * number_of_neighbours);
    for (int k = 0; k != number_of_neighbours; ++k)
    {
        MPI_Irecv(recv_buffers_[k].data(), int(recv_buffers_[k].size()), MPI_BYTE,
                  neighbour_ranks_[k], tag_halo_state, MPI_COMM_WORLD, &requests[2 * k]);
        MPI_Isend(send_buffers_[k].data(), int(send_buffers_[k].size()), MPI_BYTE,
                  neighbour_ranks_[k], tag_halo_state, MPI_COMM_WORLD, &requests[2 * k + 1]);
    }
    // TODO(overlap): this is where interior-particle kernels should be launched, with the
    // Waitall deferred until just before the boundary-particle kernels need the halo.
    // The SYCL path already has ExecutionEvent-based async dispatch to express that.
    MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
#endif // SPHINXSYS_USE_MPI

    for (size_t k = 0; k != neighbour_ranks_.size(); ++k)
    {
        size_t offset = 0;
        for (size_t v = 0; v != packers_.size(); ++v)
        {
            packers_[v]->unpack(recv_first_slot_[k], recv_counts_[k],
                                recv_buffers_[k].data() + offset);
            offset += size_t(recv_counts_[k]) * packers_[v]->bytesPerParticle();
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void HaloExchangeCK<ExecutionPolicy>::exec(Real dt)
{
    (void)dt;
    buildSendLists();
    exchangeState();
}
//=================================================================================================//
} // namespace SPH
#endif // HALO_EXCHANGE_CK_HPP
