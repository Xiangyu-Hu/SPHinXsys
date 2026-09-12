#ifndef SUBDOMAIN_EXCHANGE_HPP
#define SUBDOMAIN_EXCHANGE_HPP

#include "subdomain_exchange.h"

#include "algorithm_primitive.h"
#include "base_configuration_dynamics.h"
#include "base_particles.hpp"
#include "cell_linked_list.h" // declares ConcurrentVec, which particle_iterators.h relies on
#include "particle_iterators.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, typename DataType>
VariableExchangeBuffer<ExecutionPolicy, DataType>::VariableExchangeBuffer(
    DiscreteVariable<DataType> *variable, UnsignedInt capacity)
    : variable_(variable), capacity_(capacity)
{
    const UnsignedInt width = variable->getWidth();
    auto create = [&](const std::string &suffix)
    {
        return width == 1
                   ? buffer_keeper_.template createPtr<DiscreteVariable<DataType>>(
                         variable->Name() + suffix, capacity)
                   : buffer_keeper_.template createPtr<DiscreteVariable<DataType>>(
                         variable->Name() + suffix, capacity, MultiEntryTag{}, width);
    };

    for (int side = 0; side < NumberOfSides; ++side)
    {
        send_buffer_[side] = create("_Send" + std::to_string(side));
    }

    // Pre-touch every replica, so that a later DelegatedData() from a neighbor's
    // thread only reads an existing pointer and cannot race on the lazy allocation.
    for (int subdomain_id = 0; subdomain_id < execution::numberOfSubdomains(); ++subdomain_id)
    {
        execution::SubdomainScope scope(subdomain_id);
        for (int side = 0; side < NumberOfSides; ++side)
        {
            send_buffer_[side]->DelegatedData(ExecutionPolicy{});
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename DataType>
DataType *VariableExchangeBuffer<ExecutionPolicy, DataType>::SendBuffer(int subdomain_id, int side)
{
    execution::SubdomainScope scope(subdomain_id);
    return send_buffer_[side]->DelegatedData(ExecutionPolicy{});
}
//=================================================================================================//
template <class ExecutionPolicy, typename DataType>
void VariableExchangeBuffer<ExecutionPolicy, DataType>::reserve(UnsignedInt capacity)
{
    if (capacity <= capacity_)
    {
        return;
    }
    for (int subdomain_id = 0; subdomain_id < execution::numberOfSubdomains(); ++subdomain_id)
    {
        execution::SubdomainScope scope(subdomain_id);
        for (int side = 0; side < NumberOfSides; ++side)
        {
            send_buffer_[side]->reallocateData(ExecutionPolicy{}, capacity);
        }
    }
    capacity_ = capacity;
}
//=================================================================================================//
// Operations on the variable assembles. Each is invoked once per contained data type.
//=================================================================================================//
/** Whether a variable belongs to the selected subset of an exchange. */
template <typename DataType>
inline bool isSelectedVariable(DiscreteVariables &selected, DiscreteVariable<DataType> *variable)
{
    auto &selected_list = std::get<DataContainerAddressKeeper<DiscreteVariable<DataType>>>(selected);
    return std::find(selected_list.begin(), selected_list.end(), variable) != selected_list.end();
}

template <class ExecutionPolicy>
struct CreateExchangeBuffers
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    typename SubdomainExchange<ExecutionPolicy>::ExchangeBufferPtrs &buffer_ptrs,
                    typename SubdomainExchange<ExecutionPolicy>::ExchangeBuffers &buffers,
                    UnsignedInt capacity)
    {
        using BufferType = VariableExchangeBuffer<ExecutionPolicy, DataType>;
        auto &keeper = std::get<DataContainerUniquePtrKeeper<BufferType>>(buffer_ptrs);
        auto &buffer_list = std::get<DataContainerAddressKeeper<BufferType>>(buffers);
        for (std::size_t k = 0; k != variables.size(); ++k)
        {
            buffer_list.push_back(keeper.template createPtr<BufferType>(variables[k], capacity));
        }
    }
};

template <class ExecutionPolicy>
struct ReserveExchangeBuffers
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<VariableExchangeBuffer<ExecutionPolicy, DataType>> &buffers,
                    UnsignedInt capacity)
    {
        for (std::size_t k = 0; k != buffers.size(); ++k)
        {
            buffers[k]->reserve(capacity);
        }
    }
};

/** Gather the flagged particles of the current subdomain into its send buffer. */
template <class ExecutionPolicy>
struct PackVariablesToSendBuffer
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<VariableExchangeBuffer<ExecutionPolicy, DataType>> &buffers,
                    DiscreteVariables &selected, int subdomain_id, int side,
                    const UnsignedInt *send_index, UnsignedInt count)
    {
        if (count == 0)
        {
            return;
        }
        for (std::size_t k = 0; k != buffers.size(); ++k)
        {
            if (!isSelectedVariable(selected, buffers[k]->Variable()))
            {
                continue;
            }
            const UnsignedInt width = buffers[k]->Variable()->getWidth();
            DataType *source = buffers[k]->Variable()->DelegatedData(ExecutionPolicy{});
            DataType *destination = buffers[k]->SendBuffer(subdomain_id, side);
            particle_for(ExecutionPolicy{}, IndexRange(0, count),
                         [=](std::size_t i)
                         {
                             for (UnsignedInt entry = 0; entry < width; ++entry)
                             {
                                 destination[i * width + entry] =
                                     source[send_index[i] * width + entry];
                             }
                         });
        }
    }
};

/** Copy the neighbors' packed values into the destination slots of this subdomain. */
template <class ExecutionPolicy>
struct PullVariablesFromNeighbor
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<VariableExchangeBuffer<ExecutionPolicy, DataType>> &buffers,
                    DiscreteVariables &selected, int subdomain_id, int neighbor_subdomain,
                    int neighbor_side, UnsignedInt destination_offset, UnsignedInt count)
    {
        if (count == 0)
        {
            return;
        }
        for (std::size_t k = 0; k != buffers.size(); ++k)
        {
            if (!isSelectedVariable(selected, buffers[k]->Variable()))
            {
                continue;
            }
            const UnsignedInt width = buffers[k]->Variable()->getWidth();
            DataType *destination = buffers[k]->Variable()->DelegatedData(ExecutionPolicy{}) +
                                    destination_offset * width;
            DataType *source = buffers[k]->SendBuffer(neighbor_subdomain, neighbor_side);
            // The only backend specific call of the exchange.
            execution::copyBetweenSubdomains(ExecutionPolicy{}, subdomain_id,
                                             destination, source, count * width);
        }
    }
};

/** Remove the departing particles by moving staying particles from the tail into
 *  their slots. Holes lie below the new owned count and donors at or above it, so no
 *  slot is both read and written and the copies are independent. Every exchanged
 *  variable is moved, whether or not it is part of the current halo subset: the slot
 *  now belongs to a different particle. */
template <class ExecutionPolicy>
struct FillHolesFromTail
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<VariableExchangeBuffer<ExecutionPolicy, DataType>> &buffers,
                    const UnsignedInt *hole_index, const UnsignedInt *donor_index, UnsignedInt fill_count)
    {
        if (fill_count == 0)
        {
            return;
        }
        for (std::size_t k = 0; k != buffers.size(); ++k)
        {
            const UnsignedInt width = buffers[k]->Variable()->getWidth();
            DataType *data = buffers[k]->Variable()->DelegatedData(ExecutionPolicy{});
            particle_for(ExecutionPolicy{}, IndexRange(0, fill_count),
                         [=](std::size_t i)
                         {
                             for (UnsignedInt entry = 0; entry < width; ++entry)
                             {
                                 data[hole_index[i] * width + entry] =
                                     data[donor_index[i] * width + entry];
                             }
                         });
        }
    }
};

/** Reorder the host arrays so that each subdomain's particles are contiguous. */
struct PermuteHostVariables
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    const StdVec<UnsignedInt> &order)
    {
        StdVec<DataType> temporary;
        for (std::size_t k = 0; k != variables.size(); ++k)
        {
            const UnsignedInt width = variables[k]->getWidth();
            DataType *data = variables[k]->Data();
            temporary.assign(order.size() * width, ZeroData<DataType>::value);
            for (std::size_t i = 0; i != order.size(); ++i)
            {
                for (UnsignedInt entry = 0; entry < width; ++entry)
                {
                    temporary[i * width + entry] = data[order[i] * width + entry];
                }
            }
            std::copy(temporary.begin(), temporary.end(), data);
        }
    }
};

/** Copy one contiguous slice of the host arrays into a subdomain's replica, or back. */
template <class ExecutionPolicy>
struct StageHostSlice
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
                    int subdomain_id, UnsignedInt host_offset, UnsignedInt count, bool to_subdomain)
    {
        if (count == 0)
        {
            return;
        }
        execution::SubdomainScope scope(subdomain_id);
        for (std::size_t k = 0; k != variables.size(); ++k)
        {
            const UnsignedInt width = variables[k]->getWidth();
            DataType *host_data = variables[k]->Data() + host_offset * width;
            DataType *replica = variables[k]->DelegatedData(ExecutionPolicy{});
            if (to_subdomain)
            {
                execution::copyBetweenSubdomains(ExecutionPolicy{}, subdomain_id,
                                                 replica, host_data, count * width);
            }
            else
            {
                execution::copyBetweenSubdomains(ExecutionPolicy{}, subdomain_id,
                                                 host_data, replica, count * width);
            }
        }
    }
};
//=================================================================================================//
namespace
{
/** Initial reservation of the exchange buffers, as a fraction of the particle bound.
 *  A slab holds roughly bound/n particles, of which the halo band is the fraction
 *  halo_width / slab_thickness; a tenth is a generous start and the buffers grow on
 *  demand anyway. */
constexpr Real initial_buffer_fraction = Real(0.1);
} // namespace
//=================================================================================================//
template <class ExecutionPolicy>
SubdomainExchange<ExecutionPolicy>::SubdomainExchange(
    SlabDecomposition &decomposition, BaseParticles &particles,
    DiscreteVariables &variables_to_exchange)
    : decomposition_(decomposition), particles_(particles),
      variables_to_exchange_(variables_to_exchange),
      buffer_capacity_(std::max<UnsignedInt>(
          UnsignedInt(initial_buffer_fraction * Real(particles.ParticlesBound())), 1024))
{
    const int number_of_subdomains = execution::numberOfSubdomains();
    const UnsignedInt particles_bound = particles_.ParticlesBound();

    OperationOnDataAssemble<DiscreteVariables, CreateExchangeBuffers<ExecutionPolicy>> create_buffers;
    create_buffers(variables_to_exchange_, exchange_buffer_ptrs_, exchange_buffers_, buffer_capacity_);

    dv_send_flag_ = particles_.template addUniqueDiscreteVariable<UnsignedInt>(
        "SubdomainSendFlag", particles_bound + 1);
    dv_send_scan_ = particles_.template addUniqueDiscreteVariable<UnsignedInt>(
        "SubdomainSendScan", particles_bound + 1);
    for (int side = 0; side < NumberOfSides; ++side)
    {
        dv_send_index_[side] = particles_.template addUniqueDiscreteVariable<UnsignedInt>(
            "SubdomainSendIndex" + std::to_string(side), particles_bound);
    }
    dv_hole_index_ = particles_.template addUniqueDiscreteVariable<UnsignedInt>(
        "SubdomainHoleIndex", particles_bound);
    dv_donor_index_ = particles_.template addUniqueDiscreteVariable<UnsignedInt>(
        "SubdomainDonorIndex", particles_bound);

    send_count_.resize(number_of_subdomains, {0, 0});
    recv_count_.resize(number_of_subdomains, {0, 0});
    halo_offset_.resize(number_of_subdomains, {0, 0});
    owned_count_.resize(number_of_subdomains, 0);
    keep_count_.resize(number_of_subdomains, 0);
    fill_count_.resize(number_of_subdomains, 0);

    // The particle counters are written per subdomain through setValue(), which
    // addresses the replica of the subdomain bound to the calling thread. A replica is
    // only created on the first DelegatedData() call; before that every subdomain's
    // delegate still aliases the host value, and per-subdomain writes would overwrite
    // each other there. Touching the replicas here makes every later write land on
    // the right subdomain.
    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        execution::SubdomainScope scope(subdomain_id);
        particles_.svTotalRealParticles()->DelegatedData(ExecutionPolicy{});
        particles_.svTotalLocalParticles()->DelegatedData(ExecutionPolicy{});
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename DataType>
void SubdomainExchange<ExecutionPolicy>::addExchangeVariable(DiscreteVariable<DataType> *variable)
{
    if (particles_.template addDiscreteVariableToList<DataType>(variables_to_exchange_, variable) == nullptr)
    {
        return; // already in the set
    }
    using BufferType = VariableExchangeBuffer<ExecutionPolicy, DataType>;
    auto &keeper = std::get<DataContainerUniquePtrKeeper<BufferType>>(exchange_buffer_ptrs_);
    auto &buffer_list = std::get<DataContainerAddressKeeper<BufferType>>(exchange_buffers_);
    buffer_list.push_back(keeper.template createPtr<BufferType>(variable, buffer_capacity_));
}
//=================================================================================================//
template <class ExecutionPolicy>
UnsignedInt SubdomainExchange<ExecutionPolicy>::TotalOwnedParticles() const
{
    return std::accumulate(owned_count_.begin(), owned_count_.end(), UnsignedInt(0));
}
//=================================================================================================//
template <class ExecutionPolicy>
UnsignedInt SubdomainExchange<ExecutionPolicy>::largestSendCount() const
{
    UnsignedInt largest = 0;
    for (int subdomain_id = 0; subdomain_id < execution::numberOfSubdomains(); ++subdomain_id)
    {
        for (int side = 0; side < NumberOfSides; ++side)
        {
            largest = std::max(largest, send_count_[subdomain_id][side]);
        }
    }
    return largest;
}
//=================================================================================================//
template <class ExecutionPolicy>
Real SubdomainExchange<ExecutionPolicy>::HaloLoadFactor() const
{
    return buffer_capacity_ > 0 ? Real(largestSendCount()) / Real(buffer_capacity_) : Real(0);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::reserveBuffers(UnsignedInt capacity)
{
    if (capacity <= buffer_capacity_)
    {
        return;
    }
    buffer_capacity_ = capacity + capacity / 4;
    OperationOnDataAssemble<ExchangeBuffers, ReserveExchangeBuffers<ExecutionPolicy>> reserve_buffers;
    reserve_buffers(exchange_buffers_, buffer_capacity_);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::scatterFromHost()
{
    const SubdomainMap &map = decomposition_.getSubdomainMap();
    const int number_of_subdomains = execution::numberOfSubdomains();
    const UnsignedInt total_particles = particles_.TotalRealParticles();
    Vecd *host_position = particles_.dvParticlePosition()->Data();
    auto &scatter_variables = particles_.EvolvingVariables();

    // Group the host arrays by owner, so that each subdomain becomes one contiguous
    // range and can be staged into its replica with a single copy per variable.
    StdVec<UnsignedInt> order;
    order.reserve(total_particles);
    StdVec<UnsignedInt> host_offset(number_of_subdomains + 1, 0);
    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        for (UnsignedInt index = 0; index < total_particles; ++index)
        {
            if (map.subdomainOf(host_position[index]) == subdomain_id)
            {
                order.push_back(index);
            }
        }
        owned_count_[subdomain_id] = order.size() - host_offset[subdomain_id];
        host_offset[subdomain_id + 1] = order.size();
    }

    OperationOnDataAssemble<DiscreteVariables, PermuteHostVariables> permute;
    permute(scatter_variables, order);

    OperationOnDataAssemble<DiscreteVariables, StageHostSlice<ExecutionPolicy>> stage;
    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        stage(scatter_variables, subdomain_id, host_offset[subdomain_id],
              owned_count_[subdomain_id], true);
    }

    execution::fanOutOverSubdomains(
        ExecutionPolicy{},
        [&]()
        {
            const int subdomain_id = execution::currentSubdomainID();
            particles_.svTotalRealParticles()->setValue(owned_count_[subdomain_id]);
            particles_.svTotalLocalParticles()->setValue(owned_count_[subdomain_id]);
        });

    std::cout << "SubdomainExchange::scatterFromHost(): ";
    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        std::cout << owned_count_[subdomain_id]
                  << (subdomain_id + 1 < number_of_subdomains ? " / " : "\n");
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::gatherToHost(DiscreteVariables &variables)
{
    // Only the owned particles are gathered; the halo copies are duplicates of
    // particles owned elsewhere and would otherwise be written out twice.
    OperationOnDataAssemble<DiscreteVariables, StageHostSlice<ExecutionPolicy>> stage;
    UnsignedInt host_offset = 0;
    for (int subdomain_id = 0; subdomain_id < execution::numberOfSubdomains(); ++subdomain_id)
    {
        stage(variables, subdomain_id, host_offset,
              owned_count_[subdomain_id], false);
        host_offset += owned_count_[subdomain_id];
    }
    // Host side readers (output, restart) read the counter through the delegate of
    // subdomain 0, since the host thread is bound to subdomain 0 outside a fan-out.
    // Publish the global count there, and restore the owned count of subdomain 0
    // with finishHostAccess() before the next fan-out.
    particles_.svTotalRealParticles()->setValue(0, host_offset);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::finishHostAccess()
{
    particles_.svTotalRealParticles()->setValue(0, owned_count_[0]);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::buildSendListsOnCurrentSubdomain()
{
    const int subdomain_id = execution::currentSubdomainID();
    const SubdomainMap map = decomposition_.getSubdomainMap();
    const UnsignedInt owned_particles = owned_count_[subdomain_id];

    Vecd *position = particles_.dvParticlePosition()->DelegatedData(ExecutionPolicy{});
    UnsignedInt *send_flag = dv_send_flag_->DelegatedData(ExecutionPolicy{});
    UnsignedInt *send_scan = dv_send_scan_->DelegatedData(ExecutionPolicy{});

    for (int side = 0; side < NumberOfSides; ++side)
    {
        const int neighbor = map.neighborOf(subdomain_id, side);
        if (neighbor < 0)
        {
            send_count_[subdomain_id][side] = 0;
            continue;
        }
        UnsignedInt *send_index = dv_send_index_[side]->DelegatedData(ExecutionPolicy{});

        // The extra entry at owned_particles makes the exclusive scan return the total.
        particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles + 1),
                     [=](std::size_t i)
                     {
                         send_flag[i] = (i < owned_particles && map.inHaloBandOf(position[i], neighbor))
                                            ? UnsignedInt(1)
                                            : UnsignedInt(0);
                     });

        const UnsignedInt count = exclusive_scan(
            ExecutionPolicy{}, send_flag, send_scan, owned_particles + 1,
            typename PlusUnsignedInt<ExecutionPolicy>::type());

        particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles),
                     [=](std::size_t i)
                     {
                         if (send_flag[i] != UnsignedInt(0))
                         {
                             send_index[send_scan[i]] = i;
                         }
                     });

        send_count_[subdomain_id][side] = count;
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::publishLocalCountOnCurrentSubdomain()
{
    const int subdomain_id = execution::currentSubdomainID();
    const SubdomainMap &map = decomposition_.getSubdomainMap();

    UnsignedInt local_particles = owned_count_[subdomain_id];
    for (int side = 0; side < NumberOfSides; ++side)
    {
        const int neighbor = map.neighborOf(subdomain_id, side);
        recv_count_[subdomain_id][side] =
            neighbor < 0 ? 0 : send_count_[neighbor][SubdomainMap::oppositeSide(side)];
        halo_offset_[subdomain_id][side] = local_particles;
        local_particles += recv_count_[subdomain_id][side];
    }

    if (local_particles > particles_.ParticlesBound())
    {
        std::cout << "\n Error: subdomain " << subdomain_id << " needs " << local_particles
                  << " local particles but the particle bound is " << particles_.ParticlesBound()
                  << ". Increase the reserve of the body. \n";
        exit(1);
    }

    particles_.svTotalRealParticles()->setValue(owned_count_[subdomain_id]);
    particles_.svTotalLocalParticles()->setValue(local_particles);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::packOnCurrentSubdomain(DiscreteVariables &variables)
{
    const int subdomain_id = execution::currentSubdomainID();
    OperationOnDataAssemble<ExchangeBuffers, PackVariablesToSendBuffer<ExecutionPolicy>> pack;
    for (int side = 0; side < NumberOfSides; ++side)
    {
        const UnsignedInt count = send_count_[subdomain_id][side];
        if (count == 0)
        {
            continue;
        }
        const UnsignedInt *send_index = dv_send_index_[side]->DelegatedData(ExecutionPolicy{});
        pack(exchange_buffers_, variables, subdomain_id, side, send_index, count);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::pullFromNeighborsOnCurrentSubdomain(
    DiscreteVariables &variables)
{
    const int subdomain_id = execution::currentSubdomainID();
    const SubdomainMap &map = decomposition_.getSubdomainMap();
    OperationOnDataAssemble<ExchangeBuffers, PullVariablesFromNeighbor<ExecutionPolicy>> pull;
    for (int side = 0; side < NumberOfSides; ++side)
    {
        const int neighbor = map.neighborOf(subdomain_id, side);
        const UnsignedInt count = recv_count_[subdomain_id][side];
        if (neighbor < 0 || count == 0)
        {
            continue;
        }
        pull(exchange_buffers_, variables, subdomain_id, neighbor,
             SubdomainMap::oppositeSide(side), halo_offset_[subdomain_id][side], count);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::updateHaloPlan()
{
    execution::fanOutOverSubdomains(ExecutionPolicy{}, [&]()
                                    { buildSendListsOnCurrentSubdomain(); });
    // Barrier: every send count is now visible to the neighbors.

    reserveBuffers(largestSendCount());

    execution::fanOutOverSubdomains(ExecutionPolicy{}, [&]()
                                    { publishLocalCountOnCurrentSubdomain(); });

    refreshHalo(variables_to_exchange_);
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::refreshHalo(DiscreteVariables &variables)
{
    execution::fanOutOverSubdomains(ExecutionPolicy{}, [&]()
                                    { packOnCurrentSubdomain(variables); });
    // Barrier: all send buffers are complete before anyone reads a neighbor's.
    execution::fanOutOverSubdomains(ExecutionPolicy{}, [&]()
                                    { pullFromNeighborsOnCurrentSubdomain(variables); });
}
//=================================================================================================//
template <class ExecutionPolicy>
void SubdomainExchange<ExecutionPolicy>::migrateParticles()
{
    const SubdomainMap map = decomposition_.getSubdomainMap();

    // 1. Split the owned particles into those that stay and those that leave. A
    //    particle may only move to an adjacent subdomain within one step; a larger jump
    //    means the time step or the slab thickness is wrong, and is reported by
    //    checkConsistency() as a particle count mismatch.
    execution::fanOutOverSubdomains(
        ExecutionPolicy{},
        [&]()
        {
            const int subdomain_id = execution::currentSubdomainID();
            const UnsignedInt owned_particles = owned_count_[subdomain_id];
            Vecd *position = particles_.dvParticlePosition()->DelegatedData(ExecutionPolicy{});
            UnsignedInt *send_flag = dv_send_flag_->DelegatedData(ExecutionPolicy{});
            UnsignedInt *send_scan = dv_send_scan_->DelegatedData(ExecutionPolicy{});
            UnsignedInt *hole_index = dv_hole_index_->DelegatedData(ExecutionPolicy{});
            UnsignedInt *donor_index = dv_donor_index_->DelegatedData(ExecutionPolicy{});

            for (int side = 0; side < NumberOfSides; ++side)
            {
                const int neighbor = map.neighborOf(subdomain_id, side);
                if (neighbor < 0)
                {
                    send_count_[subdomain_id][side] = 0;
                    continue;
                }
                UnsignedInt *send_index = dv_send_index_[side]->DelegatedData(ExecutionPolicy{});
                particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles + 1),
                             [=](std::size_t i)
                             {
                                 send_flag[i] = (i < owned_particles &&
                                                 map.subdomainOf(position[i]) == neighbor)
                                                    ? UnsignedInt(1)
                                                    : UnsignedInt(0);
                             });
                const UnsignedInt count = exclusive_scan(
                    ExecutionPolicy{}, send_flag, send_scan, owned_particles + 1,
                    typename PlusUnsignedInt<ExecutionPolicy>::type());
                particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles),
                             [=](std::size_t i)
                             {
                                 if (send_flag[i] != UnsignedInt(0))
                                 {
                                     send_index[send_scan[i]] = i;
                                 }
                             });
                send_count_[subdomain_id][side] = count;
            }

            // Owned count after the departures: the departing particles are counted
            // once, whichever side they leave to.
            particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles + 1),
                         [=](std::size_t i)
                         {
                             send_flag[i] = (i < owned_particles &&
                                             map.subdomainOf(position[i]) != subdomain_id)
                                                ? UnsignedInt(1)
                                                : UnsignedInt(0);
                         });
            const UnsignedInt departing = exclusive_scan(
                ExecutionPolicy{}, send_flag, send_scan, owned_particles + 1,
                typename PlusUnsignedInt<ExecutionPolicy>::type());
            const UnsignedInt new_owned = owned_particles - departing;
            keep_count_[subdomain_id] = new_owned;

            // Holes: departing slots below the new owned count, in ascending order.
            particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles + 1),
                         [=](std::size_t i)
                         {
                             send_flag[i] = (i < new_owned &&
                                             map.subdomainOf(position[i]) != subdomain_id)
                                                ? UnsignedInt(1)
                                                : UnsignedInt(0);
                         });
            const UnsignedInt holes = exclusive_scan(
                ExecutionPolicy{}, send_flag, send_scan, owned_particles + 1,
                typename PlusUnsignedInt<ExecutionPolicy>::type());
            particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles),
                         [=](std::size_t i)
                         {
                             if (send_flag[i] != UnsignedInt(0))
                             {
                                 hole_index[send_scan[i]] = i;
                             }
                         });

            // Donors: staying particles at or above the new owned count, ascending.
            // There are exactly as many of them as there are holes.
            particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles + 1),
                         [=](std::size_t i)
                         {
                             send_flag[i] = (i >= new_owned && i < owned_particles &&
                                             map.subdomainOf(position[i]) == subdomain_id)
                                                ? UnsignedInt(1)
                                                : UnsignedInt(0);
                         });
            const UnsignedInt donors = exclusive_scan(
                ExecutionPolicy{}, send_flag, send_scan, owned_particles + 1,
                typename PlusUnsignedInt<ExecutionPolicy>::type());
            particle_for(ExecutionPolicy{}, IndexRange(0, owned_particles),
                         [=](std::size_t i)
                         {
                             if (send_flag[i] != UnsignedInt(0))
                             {
                                 donor_index[send_scan[i]] = i;
                             }
                         });
            if (holes != donors)
            {
                std::cout << "\n Error: subdomain " << subdomain_id << " has " << holes
                          << " holes but " << donors << " donors on migration. \n";
                exit(1);
            }
            fill_count_[subdomain_id] = holes;
        });

    reserveBuffers(largestSendCount());

    // 2. Pack the departing particles, then remove them by filling their slots from
    //    the tail. Packing reads the departing slots and the filling overwrites them,
    //    so the order matters; no barrier is needed between the two, since the filling
    //    touches only local memory.
    execution::fanOutOverSubdomains(
        ExecutionPolicy{},
        [&]()
        {
            const int subdomain_id = execution::currentSubdomainID();
            packOnCurrentSubdomain(variables_to_exchange_);

            OperationOnDataAssemble<ExchangeBuffers, FillHolesFromTail<ExecutionPolicy>> fill_holes;
            const UnsignedInt *hole_index = dv_hole_index_->DelegatedData(ExecutionPolicy{});
            const UnsignedInt *donor_index = dv_donor_index_->DelegatedData(ExecutionPolicy{});
            fill_holes(exchange_buffers_, hole_index, donor_index, fill_count_[subdomain_id]);
        });
    // Barrier: every send buffer holds the migrants before anyone pulls.

    // 3. Append the arrivals after the staying particles; they become owned.
    execution::fanOutOverSubdomains(
        ExecutionPolicy{},
        [&]()
        {
            const int subdomain_id = execution::currentSubdomainID();
            UnsignedInt owned_particles = keep_count_[subdomain_id];
            OperationOnDataAssemble<ExchangeBuffers, PullVariablesFromNeighbor<ExecutionPolicy>> pull;
            for (int side = 0; side < NumberOfSides; ++side)
            {
                const int neighbor = map.neighborOf(subdomain_id, side);
                const UnsignedInt count =
                    neighbor < 0 ? 0 : send_count_[neighbor][SubdomainMap::oppositeSide(side)];
                if (count == 0)
                {
                    continue;
                }
                pull(exchange_buffers_, variables_to_exchange_, subdomain_id, neighbor,
                     SubdomainMap::oppositeSide(side), owned_particles, count);
                owned_particles += count;
            }

            owned_count_[subdomain_id] = owned_particles;
            particles_.svTotalRealParticles()->setValue(owned_particles);
            particles_.svTotalLocalParticles()->setValue(owned_particles);
        });
}
//=================================================================================================//
template <class ExecutionPolicy>
std::string SubdomainExchange<ExecutionPolicy>::checkConsistency() const
{
    std::ostringstream stream;
    const int number_of_subdomains = execution::numberOfSubdomains();
    const SubdomainMap &map = decomposition_.getSubdomainMap();

    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        if (owned_count_[subdomain_id] > particles_.ParticlesBound())
        {
            stream << "subdomain " << subdomain_id << " owns " << owned_count_[subdomain_id]
                   << " particles, past the bound " << particles_.ParticlesBound();
            return stream.str();
        }
    }

    // Every owned particle must lie inside its own slab. This is the check that catches
    // a missed migration, which is otherwise silent until the physics drifts.
    for (int subdomain_id = 0; subdomain_id < number_of_subdomains; ++subdomain_id)
    {
        execution::SubdomainScope scope(subdomain_id);
        Vecd *position = const_cast<BaseParticles &>(particles_)
                             .dvParticlePosition()
                             ->DelegatedData(ExecutionPolicy{});
        if (!std::is_same<ExecutionPolicy, MultiHostPolicy>::value &&
            !std::is_same<ExecutionPolicy, SequencedMultiHostPolicy>::value)
        {
            continue; // device replicas are not host readable; use gatherToHost() first
        }
        for (UnsignedInt i = 0; i < owned_count_[subdomain_id]; ++i)
        {
            if (map.subdomainOf(position[i]) != subdomain_id)
            {
                stream << "particle " << i << " of subdomain " << subdomain_id
                       << " belongs to subdomain " << map.subdomainOf(position[i]);
                return stream.str();
            }
        }
    }
    return std::string();
}
//=================================================================================================//
} // namespace SPH
#endif // SUBDOMAIN_EXCHANGE_HPP
