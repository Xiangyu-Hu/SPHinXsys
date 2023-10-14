#include "particle_sorting.h"

#include "base_body.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
SwapSortableParticleData::SwapSortableParticleData(BaseParticles &base_particles)
    : sequence_(base_particles.sequence_),
      unsorted_id_(base_particles.unsorted_id_),
      sortable_data_(base_particles.sortable_data_) {}
//=================================================================================================//
void SwapSortableParticleData::operator()(size_t *a, size_t *b)
{
    std::swap(*a, *b);

    size_t index_a = a - sequence_.data();
    size_t index_b = b - sequence_.data();
    std::swap(unsorted_id_[index_a], unsorted_id_[index_b]);
    swap_particle_data_value_(sortable_data_, index_a, index_b);
}
//=================================================================================================//
MoveSortableParticleDeviceData::MoveSortableParticleDeviceData(BaseParticles &base_particles)
    : base_particles_(base_particles), initialized_swap_variables_(false) {}
//=================================================================================================//
void MoveSortableParticleDeviceData::operator()(size_t *index_permutation, size_t size)
{
    if (!initialized_swap_variables_)
    {
        swap_particle_device_data_value_.init(base_particles_.sortable_device_variables_);
        swap_unsorted_id_.init(base_particles_.unsorted_id_device_, base_particles_.total_real_particles_);
        initialized_swap_variables_ = true;
    }
    auto swap_events = swap_particle_device_data_value_(index_permutation, size);
    swap_events.emplace_back(swap_unsorted_id_(index_permutation, size));
    sycl::event::wait(swap_events);
}
//=================================================================================================//
ParticleSorting::ParticleSorting(BaseParticles &base_particles)
    : base_particles_(base_particles), index_sorting_device_variables_(nullptr),
      swap_sortable_particle_data_(base_particles), move_sortable_particle_device_data_(base_particles), compare_(),
      quick_sort_particle_range_(base_particles_.sequence_.data(), 0, compare_, swap_sortable_particle_data_),
      quick_sort_particle_body_() {}
//=================================================================================================//
void ParticleSorting::sortingParticleData(size_t *begin, size_t size)
{
    quick_sort_particle_range_.begin_ = begin;
    quick_sort_particle_range_.size_ = size;
    parallel_for(quick_sort_particle_range_, quick_sort_particle_body_, ap);
    updateSortedId();
}
//=================================================================================================//
void ParticleSorting::updateSortedId()
{
    const StdLargeVec<size_t> &unsorted_id = base_particles_.unsorted_id_;
    StdLargeVec<size_t> &sorted_id = base_particles_.sorted_id_;
    size_t total_real_particles = base_particles_.total_real_particles_;
    parallel_for(
        IndexRange(0, total_real_particles),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                sorted_id[unsorted_id[i]] = i;
            }
        },
        ap);
}
//=================================================================================================//
void ParticleSorting::updateSortedDeviceId() const
{
    size_t *unsorted_id_device = base_particles_.unsorted_id_device_;
    size_t *sorted_id_device = base_particles_.sorted_id_device_;
    size_t total_real_particles = base_particles_.total_real_particles_;
    execution::executionQueue.getQueue()
        .parallel_for(execution::executionQueue.getUniformNdRange(total_real_particles),
                      [=](sycl::nd_item<1> item)
                      {
                          size_t i = item.get_global_id();
                          if (i < total_real_particles)
                              sorted_id_device[unsorted_id_device[i]] = i;
                      })
        .wait();
}
//=================================================================================================//
template <>
void ParticleSorting::sortingParticleData(size_t *begin, size_t size, execution::ParallelSYCLDevicePolicy execution_policy)
{
    if (!index_sorting_device_variables_)
        index_sorting_device_variables_ = allocateDeviceData<size_t>(size);

    sort_by_key(begin, index_sorting_device_variables_, size, execution::executionQueue.getQueue(), 256, 4, [](size_t *data, size_t idx)
    { return idx; }).wait();

    move_sortable_particle_device_data_(index_sorting_device_variables_, size);

    updateSortedDeviceId();
}
//=================================================================================================//
size_t split_count(bool bit, sycl::nd_item<1> &item)
{
    const auto group_range = item.get_local_range().size();
    const size_t id = item.get_local_id();

    // Count 1s held by lower-numbered work-items and current work-item
    size_t true_before = sycl::inclusive_scan_over_group(item.get_group(), static_cast<size_t>(bit),
                                                         sycl::plus<size_t>{}); // prefix sum over group

    // Total number of 0s
    size_t false_totals = sycl::group_broadcast(item.get_group(), group_range - true_before,
                                                group_range - 1);

    // Return work-item's rank
    return bit ? true_before - 1 + false_totals : id - true_before;
}
//=================================================================================================//
size_t get_digit(size_t key, size_t d, size_t radix_bits)
{
    return (key >> d * radix_bits) & ((1ul << radix_bits) - 1);
}
//=================================================================================================//
size_t get_bit(size_t key, size_t b)
{
    return (key >> b) & 1;
}
//=================================================================================================//
} // namespace SPH
