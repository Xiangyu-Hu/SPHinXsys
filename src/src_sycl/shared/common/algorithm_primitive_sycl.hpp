#include "algorithm_primitive_sycl.h"

#ifndef ALGORITHM_PRIMITIVE_SYCL_HPP
#define ALGORITHM_PRIMITIVE_SYCL_HPP
namespace SPH
{
//=================================================================================================//
template <class DataType>
UnsignedInt DeviceRadixSort<DataType>::split_count(bool bit, sycl::nd_item<1> &item)
{
    const auto group_range = item.get_local_range().size();
    const UnsignedInt id = item.get_local_id();

    // Count 1s held by lower-numbered work-items and current work-item
    UnsignedInt true_before =
        sycl::inclusive_scan_over_group(item.get_group(), static_cast<UnsignedInt>(bit),
                                        sycl::plus<>{}); // prefix sum over group

    // Total number of 0s
    UnsignedInt false_totals =
        sycl::group_broadcast(item.get_group(), group_range - true_before,
                              group_range - 1);

    // Return work-item's rank
    return bit ? true_before - 1 + false_totals : id - true_before;
}
//=================================================================================================//
template <class DataType>
UnsignedInt DeviceRadixSort<DataType>::get_digit(
    UnsignedInt key, UnsignedInt d, UnsignedInt radix_bits)
{
    return (key >> d * radix_bits) & ((static_cast<UnsignedInt>(1) << radix_bits) - 1);
}
//=================================================================================================//
template <class DataType>
UnsignedInt DeviceRadixSort<DataType>::get_bit(UnsignedInt key, UnsignedInt b)
{
    return (key >> b) & 1;
}
//=================================================================================================//
template <class DataType>
UnsignedInt DeviceRadixSort<DataType>::find_max_element(
    const UnsignedInt *data, UnsignedInt size, UnsignedInt identity)
{
    UnsignedInt result = identity;
    auto &sycl_queue = execution_instance.getQueue();
    {
        sycl::buffer<UnsignedInt> buffer_result(&result, 1);
        sycl_queue.submit(
                      [&](sycl::handler &cgh)
                      {
                              auto reduction_operator = sycl::reduction(buffer_result, cgh, sycl::maximum<>());
                              cgh.parallel_for(execution_instance.getUniformNdRange(size), reduction_operator,
                                               [=](sycl::nd_item<1> item, auto &reduction)
                                               {
                                                   if(item.get_global_id() < size)
                                                       reduction.combine(data[item.get_global_linear_id()]);
                                               }); })
            .wait_and_throw();
    }
    return result;
}
//=================================================================================================//
template <class DataType>
void DeviceRadixSort<DataType>::resize(UnsignedInt data_size, UnsignedInt radix_bits, UnsignedInt workgroup_size)
{
    data_size_ = data_size;
    radix_bits_ = radix_bits;
    workgroup_size_ = workgroup_size;
    uniform_case_masking_ = data_size % workgroup_size;
    uniform_global_size_ = uniform_case_masking_ ? (data_size / workgroup_size + 1) * workgroup_size : data_size;
    kernel_range_ = {uniform_global_size_, workgroup_size};
    workgroups_ = kernel_range_.get_group_range().size();

    radix_ = static_cast<UnsignedInt>(1) << radix_bits; // radix = 2^b

    sycl::range<2> buckets_column_major_range = {static_cast<size_t>(radix_), static_cast<size_t>(workgroups_)},
                   buckets_row_major_range = {static_cast<size_t>(workgroups_), static_cast<size_t>(radix_)};
    // Each entry contains global number of digits with the same value
    // Column-major, so buckets offsets can be computed by just applying a scan over it
    global_buckets_ = std::make_unique<sycl::buffer<UnsignedInt, 2>>(buckets_column_major_range);
    // Each entry contains global number of digits with the same and lower values
    global_buckets_offsets_ = std::make_unique<sycl::buffer<UnsignedInt, 2>>(buckets_column_major_range);
    local_buckets_offsets_buffer_ = std::make_unique<sycl::buffer<UnsignedInt, 1>>(buckets_row_major_range.size()); // save state of local accessor
    data_swap_buffer_ = std::make_unique<sycl::buffer<SortablePair>>(uniform_global_size_);                         // temporary memory for swapping
    // Keep extra values to be swapped when kernel range has been made uniform
    uniform_extra_swap_buffer_ = std::make_unique<sycl::buffer<SortablePair>>(uniform_global_size_ - data_size);
}
//=================================================================================================//
template <class DataType>
sycl::event DeviceRadixSort<DataType>::sort_by_key(
    UnsignedInt *keys, DataType *data, UnsignedInt data_size, sycl::queue &queue,
    UnsignedInt workgroup_size, UnsignedInt radix_bits)
{
    if (data_size_ != data_size || radix_bits_ != radix_bits || workgroup_size_ != workgroup_size)
        resize(data_size, radix_bits, workgroup_size);

    // Largest key, increased by 1 if the workgroup is not homogeneous with the data vector,
    // the new maximum will be used for those work-items out of data range, that will then be excluded once sorted
    const UnsignedInt max_key = find_max_element(keys, data_size, 0ul) + (uniform_case_masking_ ? 1 : 0);
    const UnsignedInt bits_max_key = std::floor(std::log2(max_key)) + 1.0;                                    // bits needed to represent max_key
    const UnsignedInt length = max_key ? bits_max_key / radix_bits + (bits_max_key % radix_bits ? 1 : 0) : 1; // max number of radix digits

    sycl::event sort_event{};
    for (int digit = 0; digit < length; ++digit)
    {

        auto buckets_event = queue.submit(
            [&](sycl::handler &cgh)
            {
                cgh.depends_on(sort_event);
                auto data_swap_acc = data_swap_buffer_->get_access(cgh, sycl::write_only, sycl::no_init);
                auto local_buckets = sycl::local_accessor<UnsignedInt>(radix_, cgh);
                auto local_output = sycl::local_accessor<SortablePair>(kernel_range_.get_local_range(), cgh);
                auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_write, sycl::no_init);
                auto local_buckets_offsets_accessor = local_buckets_offsets_buffer_->get_access(cgh, sycl::write_only, sycl::no_init);

                cgh.parallel_for(
                    kernel_range_, 
                    [=, radix=radix_](sycl::nd_item<1> item) {
                        const UnsignedInt workgroup = item.get_group_linear_id(),
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
                        for (UnsignedInt b = 1; b < radix_bits; ++b) 
                        {   // sorting remaining bits
                            item.barrier(sycl::access::fence_space::local_space);
                            number = local_output[item.get_local_id()];
                            radix_digit = get_digit(number.first, digit, radix_bits);

                            rank = split_count(get_bit(radix_digit, b), item);
                            local_output[rank] = number;
                        }

                        // Initialize local buckets to zero, since they are uninitialized by default
                        for (UnsignedInt r = 0; r < radix; ++r)
                            local_buckets[r] = 0;

                        item.barrier(sycl::access::fence_space::local_space);
                        {
                            sycl::atomic_ref<UnsignedInt, sycl::memory_order_relaxed, sycl::memory_scope_work_group,
                                             sycl::access::address_space::local_space>
                                bucket_r{local_buckets[radix_digit]};
                            ++bucket_r;
                            item.barrier(sycl::access::fence_space::local_space);
                        }

                        // Save local buckets to global memory, with one row per work-group (in column-major order)
                        for (UnsignedInt r = 0; r < radix; ++r)
                            global_buckets_accessor[r][workgroup] = local_buckets[r];

                        if (global_id < data_size)
                            data_swap_acc[workgroup_size * workgroup + rank] = number; // save local sorting back to data

                        // Compute local buckets offsets
                        UnsignedInt *begin = local_buckets.get_multi_ptr<sycl::access::decorated::no>().get(), *end = begin + radix,
                                    *outBegin = local_buckets_offsets_accessor.get_multi_ptr<sycl::access::decorated::no>().get() + workgroup * radix;
                        sycl::joint_exclusive_scan(item.get_group(), begin, end, outBegin, sycl::plus<>{});
                                                               }); });

        // Global synchronization to make sure that all locally computed buckets have been copied to global memory

        sycl::event scan_event = queue.submit(
            [&](sycl::handler &cgh)
            {
                cgh.depends_on(buckets_event);
                auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_only);
                auto global_buckets_offsets_accessor = global_buckets_offsets_->get_access(cgh, sycl::write_only);
                cgh.parallel_for(
                    kernel_range_,
                    [=](sycl::nd_item<1> item)
                    {
                        // Compute global buckets offsets
                        if (item.get_group_linear_id() == 0)
                        {
                            const UnsignedInt *begin = global_buckets_accessor.get_multi_ptr<sycl::access::decorated::no>().get(),
                                              *end = begin + global_buckets_accessor.size();
                            sycl::joint_exclusive_scan(item.get_group(), begin, end,
                                                       global_buckets_offsets_accessor.get_multi_ptr<sycl::access::decorated::no>().get(),
                                                       sycl::plus<>{});
                        }
                    }); });

        sort_event = queue.submit(
            [&](sycl::handler &cgh)
            {
                cgh.depends_on(scan_event);
                auto data_swap_acc = data_swap_buffer_->get_access(cgh, sycl::read_only);
                auto global_buckets_accessor = global_buckets_->get_access(cgh, sycl::read_only);
                auto global_buckets_offsets_accessor = global_buckets_offsets_->get_access(cgh, sycl::read_write);
                auto local_buckets_offsets_accessor = local_buckets_offsets_buffer_->get_access(cgh, sycl::read_only);
                cgh.parallel_for(
                    kernel_range_,
                    [=, radix = radix_](sycl::nd_item<1> item)
                    {
                        // Compute global buckets offsets
                        const UnsignedInt *begin = global_buckets_accessor.get_multi_ptr<sycl::access::decorated::no>().get(),
                                          *end = begin + global_buckets_accessor.size();
                        sycl::joint_exclusive_scan(item.get_group(), begin, end,
                                                   global_buckets_offsets_accessor.get_multi_ptr<sycl::access::decorated::no>().get(),
                                                   sycl::plus<>{});

                        // Mask only relevant indexes. All max_keys added to homogenize the computations
                        // should be owned by work-items with global_id >= data_size
                        if (item.get_global_id() < data_size)
                        {
                            // Retrieve position and sorted data from swap memory
                            const UnsignedInt rank = item.get_local_id(), workgroup = item.get_group_linear_id();
                            const SortablePair number = data_swap_acc[workgroup_size * workgroup + rank];
                            const UnsignedInt radix_digit = get_digit(number.first, digit, radix_bits);

                            // Compute sorted position based on global and local buckets
                            const UnsignedInt data_offset = global_buckets_offsets_accessor[radix_digit][workgroup] + rank -
                                                            local_buckets_offsets_accessor[workgroup * radix + radix_digit];

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
#endif // ALGORITHM_PRIMITIVE_SYCL_HPP