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

    device_radix_sorting.sort_by_key(begin, index_sorting_device_variables_, size, execution::executionQueue.getQueue(), 512, 4).wait();

    move_sortable_particle_device_data_(index_sorting_device_variables_, size);

    updateSortedDeviceId();
}
//=================================================================================================//
template <class ValueType>
SYCL_EXTERNAL size_t DeviceRadixSort<ValueType>::split_count(bool bit, sycl::nd_item<1> &item)
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
template <class ValueType>
size_t DeviceRadixSort<ValueType>::get_digit(size_t key, size_t d, size_t radix_bits)
{
    return (key >> d * radix_bits) & ((1ul << radix_bits) - 1);
}
//=================================================================================================//
template <class ValueType>
size_t DeviceRadixSort<ValueType>::get_bit(size_t key, size_t b)
{
    return (key >> b) & 1;
}
//=================================================================================================//
template <class ValueType>
size_t DeviceRadixSort<ValueType>::find_max_element(const size_t *data, size_t size, size_t identity)
{
    size_t result = identity;
    auto &sycl_queue = execution::executionQueue.getQueue();
    {
        sycl::buffer<size_t> buffer_result(&result, 1);
        sycl_queue.submit([&](sycl::handler &cgh)
                          {
                              auto reduction_operator = sycl::reduction(buffer_result, cgh, sycl::maximum<>());
                              cgh.parallel_for(execution::executionQueue.getUniformNdRange(size), reduction_operator,
                                               [=](sycl::nd_item<1> item, auto& reduction) {
                                                   if(item.get_global_id() < size)
                                                       reduction.combine(data[item.get_global_linear_id()]);
                                               }); })
            .wait_and_throw();
    }
    return result;
}
//=================================================================================================//
template <class ValueType>
void DeviceRadixSort<ValueType>::resize(size_t data_size, size_t radix_bits, size_t workgroup_size)
{
    data_size_ = data_size;
    radix_bits_ = radix_bits;
    workgroup_size_ = workgroup_size;
    uniform_case_masking_ = data_size % workgroup_size;
    uniform_global_size_ = uniform_case_masking_ ? (data_size / workgroup_size + 1) * workgroup_size : data_size;
    kernel_range_ = {uniform_global_size_, workgroup_size};
    workgroups_ = kernel_range_.get_group_range().size();

    radix_ = 1ul << radix_bits;  // radix = 2^b

    sycl::range<2> buckets_column_major_range = {radix_, workgroups_}, buckets_row_major_range = {workgroups_, radix_};
    // Each entry contains global number of digits with the same value
    // Column-major, so buckets offsets can be computed by just applying a scan over it
    global_buckets_ = std::make_unique<sycl::buffer<size_t, 2>>(buckets_column_major_range);
    // Each entry contains global number of digits with the same and lower values
    global_buckets_offsets_ = std::make_unique<sycl::buffer<size_t, 2>>(buckets_column_major_range);
    local_buckets_offsets_buffer_ = std::make_unique<sycl::buffer<size_t, 2>>(buckets_row_major_range);  // save state of local accessor
    data_swap_buffer_ = std::make_unique<sycl::buffer<SortablePair>>(uniform_global_size_);  // temporary memory for swapping
    // Keep extra values to be swapped when kernel range has been made uniform
    uniform_extra_swap_buffer_ = std::make_unique<sycl::buffer<SortablePair>>(uniform_global_size_ - data_size);
}
//=================================================================================================//
template <class ValueType>
sycl::event DeviceRadixSort<ValueType>::sort_by_key(size_t *keys, ValueType *data, size_t data_size, sycl::queue &queue, size_t workgroup_size, size_t radix_bits)
{
    if(data_size_ != data_size || radix_bits_ != radix_bits || workgroup_size_ != workgroup_size)
        resize(data_size, radix_bits, workgroup_size);

    // Largest key, increased by 1 if the workgroup is not homogeneous with the data vector,
    // the new maximum will be used for those work-items out of data range, that will then be excluded once sorted
    const size_t max_key = find_max_element(keys, data_size, 0ul) + (uniform_case_masking_ ? 1 : 0);
    const size_t bits_max_key = std::floor(std::log2(max_key)) + 1.0;                                    // bits needed to represent max_key
    const size_t length = max_key ? bits_max_key / radix_bits + (bits_max_key % radix_bits ? 1 : 0) : 1; // max number of radix digits

    sycl::event sort_event{};
    for (int digit = 0; digit < length; ++digit)
    {

        auto buckets_event = queue.submit([&](sycl::handler &cgh)
                                          {
                                              cgh.depends_on(sort_event);
                                              auto data_swap_acc = data_swap_buffer_->get_access(cgh, sycl::write_only, sycl::no_init);
                                              auto local_buckets = sycl::local_accessor<size_t>(radix_, cgh);
                                              auto local_output = sycl::local_accessor<SortablePair>(kernel_range_.get_local_range(), cgh);
                                              auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_write, sycl::no_init);
                                              auto local_buckets_offsets_accessor = local_buckets_offsets_buffer_->get_access(cgh, sycl::write_only,
                                                                                                                              sycl::no_init);

                                              cgh.parallel_for(kernel_range_, [=, radix=radix_](sycl::nd_item<1> item) {
                                                                   const size_t workgroup = item.get_group_linear_id(),
                                                                                global_id = item.get_global_id();

                                                                   SortablePair number;
                                                                   // Initialize key-data pair, with masking in case of non-homogeneous data_size/workgroup_size
                                                                   if(global_id < data_size)
                                                                       number = {keys[global_id],
                                                                                 // Give possibility to initialize data here to avoid calling
                                                                                 // another kernel before sort_by_key in order to initialize it
                                                                                 digit ? data[global_id] : global_id};
                                                                   else  // masking extra indexes
                                                                       // Initialize exceeding values to the largest key considered
                                                                       number.first = (1 << bits_max_key) - 1;  // max key for given number of bits


                                                                   // Locally sort digit with split primitive
                                                                   auto radix_digit = get_digit(number.first, digit, radix_bits);
                                                                   auto rank = split_count(get_bit(radix_digit, 0), item);  // sorting first bit
                                                                   local_output[rank] = number;
                                                                   for (size_t b = 1; b < radix_bits; ++b) {  // sorting remaining bits
                                                                       item.barrier(sycl::access::fence_space::local_space);
                                                                       number = local_output[item.get_local_id()];
                                                                       radix_digit = get_digit(number.first, digit, radix_bits);

                                                                       rank = split_count(get_bit(radix_digit, b), item);
                                                                       local_output[rank] = number;
                                                                   }

                                                                   // Initialize local buckets to zero, since they are uninitialized by default
                                                                   for (size_t r = 0; r < radix; ++r)
                                                                       local_buckets[r] = 0;

                                                                   item.barrier(sycl::access::fence_space::local_space);
                                                                   {
                                                                       sycl::atomic_ref<size_t, sycl::memory_order_relaxed, sycl::memory_scope_work_group,
                                                                                        sycl::access::address_space::local_space> bucket_r{local_buckets[radix_digit]};
                                                                       ++bucket_r;
                                                                       item.barrier(sycl::access::fence_space::local_space);
                                                                   }

                                                                   // Save local buckets to global memory, with one row per work-group (in column-major order)
                                                                   for (size_t r = 0; r < radix; ++r)
                                                                       global_buckets_accessor[r][workgroup] = local_buckets[r];

                                                                   if(global_id < data_size)
                                                                       data_swap_acc[workgroup_size * workgroup + rank] = number;  // save local sorting back to data

                                                                   // Compute local buckets offsets
                                                                   size_t *begin = local_buckets.get_pointer(), *end = begin + radix,
                                                                          *outBegin = local_buckets_offsets_accessor.get_pointer().get() + workgroup * radix;
                                                                   sycl::joint_exclusive_scan(item.get_group(), begin, end, outBegin, sycl::plus<size_t>{});
                                                               }); });

        // Global synchronization to make sure that all locally computed buckets have been copied to global memory

        sycl::event scan_event = queue.submit([&](sycl::handler &cgh) {
                                                  cgh.depends_on(buckets_event);
                                                  auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_only);
                                                  auto global_buckets_offsets_accessor = global_buckets_offsets_->get_access(cgh, sycl::write_only);
                                                  cgh.parallel_for(kernel_range_, [=](sycl::nd_item<1> item) {
                                                                       // Compute global buckets offsets
                                                                       if(item.get_group_linear_id() == 0) {
                                                                           size_t *begin = global_buckets_accessor.get_pointer(), *end = begin + global_buckets_accessor.size();
                                                                           sycl::joint_exclusive_scan(item.get_group(), begin, end,
                                                                                                      global_buckets_offsets_accessor.get_pointer(), sycl::plus<size_t>{});
                                                                       }
                                                                   });
                                              });

        sort_event = queue.submit([&](sycl::handler &cgh)
                                  {
                                      cgh.depends_on(scan_event);
                                      auto data_swap_acc = data_swap_buffer_->get_access(cgh, sycl::read_only);
                                      auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_only);
                                      auto global_buckets_offsets_accessor = global_buckets_offsets_->get_access(cgh, sycl::read_write);
                                      auto local_buckets_offsets_accessor = local_buckets_offsets_buffer_->get_access(cgh, sycl::read_only);
                                      cgh.parallel_for(kernel_range_, [=](sycl::nd_item<1> item) {
                                                           // Compute global buckets offsets
                                                           size_t *begin = global_buckets_accessor.get_pointer(), *end = begin + global_buckets_accessor.size();
                                                           sycl::joint_exclusive_scan(item.get_group(), begin, end,
                                                                                      global_buckets_offsets_accessor.get_pointer(), sycl::plus<size_t>{});

                                                           // Mask only relevant indexes. All max_keys added to homogenize the computations
                                                           // should be owned by work-items with global_id >= data_size
                                                           if(item.get_global_id() < data_size) {
                                                               // Retrieve position and sorted data from swap memory
                                                               const size_t rank = item.get_local_id(), workgroup = item.get_group_linear_id();
                                                               const SortablePair number = data_swap_acc[workgroup_size * workgroup + rank];
                                                               const size_t radix_digit = get_digit(number.first, digit, radix_bits);

                                                               // Compute sorted position based on global and local buckets
                                                               const size_t data_offset = global_buckets_offsets_accessor[radix_digit][workgroup] + rank -
                                                                                          local_buckets_offsets_accessor[workgroup][radix_digit];

                                                               // Copy to original data pointers
                                                               keys[data_offset] = number.first;
                                                               data[data_offset] = number.second;
                                                           }
                                                       }); });
    }
    return sort_event;
}
//=================================================================================================//
} // namespace SPH
