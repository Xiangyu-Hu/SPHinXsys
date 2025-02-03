#ifndef MESH_WITH_DATA_PACKAGE_HPP
#define MESH_WITH_DATA_PACKAGE_HPP

#include "mesh_with_data_packages.h"
#include "execution_sycl.h"

template <int PKG_SIZE>
IndexHandler* MeshWithGridDataPackages<PKG_SIZE>::
    getIndexHandler(const ParallelDevicePolicy &par_device) const
{
    if(device_index_handler_ == nullptr)
    {
        device_index_handler_ = allocateDeviceOnly<IndexHandler>(1);
        copyToDevice(index_handler_, device_index_hadnler_, 1);
    }

    return device_index_handler_;
}

#endif