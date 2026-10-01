#include "device_environment_sycl.h"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace SPH
{
namespace execution
{
//=================================================================================================//
void DeviceEnvironment::initialize(int requested_devices)
{
    if (initialized_)
    {
        std::cout << "\n Warning: DeviceEnvironment is already initialized with "
                  << number_of_devices_ << " device(s); the request is ignored."
                  << " Call initialize() before any particle data is generated. \n";
        return;
    }

    std::vector<sycl::device> gpus;
    try
    {
        gpus = sycl::device::get_devices(sycl::info::device_type::gpu);
    }
    catch (const sycl::exception &)
    {
        gpus.clear();
    }

    if (gpus.empty())
    { // no GPU exposed: keep the previous single default-device behavior
        gpus.push_back(sycl::device(sycl::default_selector_v));
    }

    int available = static_cast<int>(gpus.size());
    number_of_devices_ = requested_devices > 0 ? std::min(requested_devices, available) : available;
    if (number_of_devices_ > MaxSubdomains)
    {
        std::cout << "\n Warning: " << number_of_devices_ << " devices found but the build "
                  << "supports at most " << MaxSubdomains
                  << "; raise execution::MaxSubdomains to use them all. \n";
        number_of_devices_ = MaxSubdomains;
    }

    devices_.assign(gpus.begin(), gpus.begin() + number_of_devices_);
    // A single context over all devices keeps USM pointers mutually valid.
    context_ = makeUnique<sycl::context>(devices_);

    queues_.clear();
    work_group_sizes_.clear();
    for (int device_id = 0; device_id < number_of_devices_; ++device_id)
    {
        queues_.emplace_back(makeUnique<sycl::queue>(*context_, devices_[device_id]));
        const unsigned long max_workgroup_size =
            devices_[device_id].get_info<sycl::info::device::max_work_group_size>();
        work_group_sizes_.push_back(std::min<unsigned long>(max_workgroup_size, 64UL));
    }

    enablePeerAccess();

    numberOfSubdomainsRef() = number_of_devices_;
    initialized_ = true;

    if (number_of_devices_ > 1)
    {
        worker_pool_ = makeUnique<SubdomainWorkerPool>(number_of_devices_);
    }

    std::cout << describe() << std::endl;
}
//=================================================================================================//
void DeviceEnvironment::enablePeerAccess()
{
    peer_access_.assign(number_of_devices_ * number_of_devices_, 0);
    for (int i = 0; i < number_of_devices_; ++i)
    {
        peer_access_[i * number_of_devices_ + i] = 1;
    }
#ifdef SYCL_EXT_ONEAPI_PEER_ACCESS
    for (int dst = 0; dst < number_of_devices_; ++dst)
    {
        for (int src = 0; src < number_of_devices_; ++src)
        {
            if (src == dst)
                continue;
            try
            {
                const bool can_access = devices_[dst].ext_oneapi_can_access_peer(
                    devices_[src], sycl::ext::oneapi::peer_access::access_supported);
                if (can_access)
                {
                    devices_[dst].ext_oneapi_enable_peer_access(devices_[src]);
                    peer_access_[src * number_of_devices_ + dst] = 1;
                }
            }
            catch (const sycl::exception &e)
            { // not fatal: copies then stage through the host
                std::cout << "\n Note: peer access " << src << " -> " << dst
                          << " unavailable (" << e.what() << "). \n";
            }
        }
    }
#endif // SYCL_EXT_ONEAPI_PEER_ACCESS
}
//=================================================================================================//
bool DeviceEnvironment::peerAccessEnabled(int src_device, int dst_device)
{
    ensureInitialized();
    return peer_access_[src_device * number_of_devices_ + dst_device] != 0;
}
//=================================================================================================//
void DeviceEnvironment::synchronizeAllDevices()
{
    ensureInitialized();
    for (int device_id = 0; device_id < number_of_devices_; ++device_id)
    {
        queues_[device_id]->wait_and_throw();
    }
}
//=================================================================================================//
std::string DeviceEnvironment::describe()
{
    std::ostringstream stream;
    stream << "SPHinXsys device environment: " << number_of_devices_ << " device(s)\n";
    for (int device_id = 0; device_id < number_of_devices_; ++device_id)
    {
        stream << "  [" << device_id << "] "
               << devices_[device_id].get_info<sycl::info::device::name>()
               << ", work group size " << work_group_sizes_[device_id] << "\n";
    }
    if (number_of_devices_ > 1)
    {
        stream << "  peer access matrix (row: source, column: destination)\n";
        for (int src = 0; src < number_of_devices_; ++src)
        {
            stream << "    ";
            for (int dst = 0; dst < number_of_devices_; ++dst)
            {
                stream << (peer_access_[src * number_of_devices_ + dst] ? "1 " : "0 ");
            }
            stream << "\n";
        }
    }
    return stream.str();
}
//=================================================================================================//
} // namespace execution
} // namespace SPH
