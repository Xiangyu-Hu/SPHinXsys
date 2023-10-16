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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	neighborhood.h
 * @brief 	There are the classes for particle neighborhood.
 * It saves the information for carrying out inter-particle (or pair) interaction,
 * and also considered as the local configuration of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "all_kernels.h"
#include "base_data_package.h"
#include "sph_data_containers.h"
#include "execution_selector.hpp"
#include "device_executable.hpp"

namespace SPH
{

class SPHBody;
class BodyPart;
class SPHAdaptation;
class NeighborhoodDevice;

/**
 * @class Neighborhood
 * @brief A neighborhood around particle i.
 */
class Neighborhood
{
  public:
    size_t current_size_;   /**< the current number of neighbors */
    size_t allocated_size_; /**< the limit of neighbors does not require memory allocation  */

    StdLargeVec<size_t> j_;      /**< index of the neighbor particle. */
    StdLargeVec<Real> W_ij_;     /**< kernel value or particle volume contribution */
    StdLargeVec<Real> dW_ijV_j_; /**< derivative of kernel function or inter-particle surface contribution */
    StdLargeVec<Real> r_ij_;     /**< distance between j and i. */
    StdLargeVec<Vecd> e_ij_;     /**< unit vector pointing from j to i or inter-particle surface direction */

    Neighborhood() : current_size_(0), allocated_size_(0){};
    ~Neighborhood(){};

    void removeANeighbor(size_t neighbor_n);

    inline size_t current_size() const { return current_size_; }

    Neighborhood& operator=(const NeighborhoodDevice& device);
};
using ParticleConfiguration = StdLargeVec<Neighborhood>;

class NeighborhoodDevice {
  public:
    size_t* current_size_;	/**< the current number of neighbors */
    size_t allocated_size_;  /**< the limit of neighbors does not require memory allocation  */

    size_t *j_;  /**< index of the neighbor particle. */
    DeviceReal *W_ij_;  /**< kernel value or particle volume contribution */
    DeviceReal *dW_ijV_j_;  /**< derivative of kernel function or inter-particle surface contribution */
    DeviceReal *r_ij_;  /**< distance between j and i. */
    DeviceVecd *e_ij_;  /**< unit vector pointing from j to i or inter-particle surface direction */

    inline size_t current_size() const { return *current_size_; }

    NeighborhoodDevice();
    ~NeighborhoodDevice();

    NeighborhoodDevice& operator=(const Neighborhood& host);
};

class NeighborBuilderKernel {
  public:
    explicit NeighborBuilderKernel(Kernel &kernel) : smoothing_kernel(kernel) {}

  protected:
    template<class NeighborhoodType, class RealType, class VecType>
    void initializeNeighbor(NeighborhoodType &neighborhood, const size_t current_size, const RealType &distance,
                            const VecType &displacement, const size_t index_j, const RealType &Vol_j) const
    {
        neighborhood.j_[current_size] = index_j;
        neighborhood.W_ij_[current_size] = smoothing_kernel.W(distance, displacement);
        neighborhood.dW_ijV_j_[current_size] = smoothing_kernel.dW(distance, displacement) * Vol_j;
        neighborhood.r_ij_[current_size] = distance;
        neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
    }

    DeviceKernelWendlandC2 smoothing_kernel;
};

/**
 * @class NeighborBuilder
 * @brief Base class for building a neighbor particle j around particles i.
 */
class NeighborBuilder
{
  protected:
    Kernel *kernel_;
    //----------------------------------------------------------------------
    //	Below are for constant smoothing length.
    //----------------------------------------------------------------------
    void createNeighbor(Neighborhood &neighborhood, const Real &distance,
                        const Vecd &displacement, size_t j_index, const Real &Vol_j);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                            const Vecd &displacement, size_t j_index, const Real &Vol_j);
    //----------------------------------------------------------------------
    //	Below are for variable smoothing length.
    //----------------------------------------------------------------------
    void createNeighbor(Neighborhood &neighborhood, const Real &distance,
                        const Vecd &displacement, size_t j_index, const Real &Vol_j, Real i_h_ratio, Real h_ratio_min);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                            const Vecd &displacement, size_t j_index, const Real &Vol_j, Real i_h_ratio, Real h_ratio_min);

  public:
    NeighborBuilder() : kernel_(nullptr){};
    virtual ~NeighborBuilder(){};
};


class NeighborBuilderInnerKernel : public NeighborBuilderKernel {
  public:
    NeighborBuilderInnerKernel(Kernel& kernel) : NeighborBuilderKernel(kernel) {}

    SYCL_EXTERNAL void operator()(NeighborhoodDevice &neighborhood, const DeviceVecd &pos_i, const size_t index_i,
                    const size_t index_j, const DeviceVecd pos_j, const DeviceReal Vol_j) const;
};

/**
 * @class NeighborBuilderInner
 * @brief A inner neighbor builder functor.
 */
class NeighborBuilderInner : public NeighborBuilder,
                             public execution::DeviceExecutable<NeighborBuilderInner, NeighborBuilderInnerKernel>
{
  public:
    explicit NeighborBuilderInner(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);
};

/**
 * @class NeighborBuilderInnerAdaptive
 * @brief A inner neighbor builder functor when the particles have different smoothing lengths.
 */
class NeighborBuilderInnerAdaptive : public NeighborBuilder
{
  public:
    explicit NeighborBuilderInnerAdaptive(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

  protected:
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class NeighborBuilderSelfContact
 * @brief A self-contact neighbor builder functor.
 */
class NeighborBuilderSelfContact : public NeighborBuilder
{
  public:
    explicit NeighborBuilderSelfContact(SPHBody &body);
    virtual ~NeighborBuilderSelfContact(){};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

  protected:
    StdLargeVec<Vecd> &pos0_;
};


class NeighborBuilderContactKernel : public NeighborBuilderKernel {
  public:
    NeighborBuilderContactKernel(Kernel& kernel) : NeighborBuilderKernel(kernel) {}

    void operator()(NeighborhoodDevice &neighborhood, const DeviceVecd &pos_i, const size_t index_i,
                    const size_t index_j, const DeviceVecd pos_j, const DeviceReal Vol_j) const
    {
        const DeviceVecd displacement = pos_i - pos_j;
        const DeviceReal distance = sycl::length(displacement);
        if (distance < smoothing_kernel.CutOffRadius())
        {
            auto current_size_atomic = sycl::atomic_ref<size_t, sycl::memory_order::relaxed,
                                                        sycl::memory_scope::device,
                                                        sycl::access::address_space::global_space>(*neighborhood.current_size_);
            initializeNeighbor(neighborhood, current_size_atomic++, distance, displacement, index_j, Vol_j);
        }
    }
};

/**
 * @class NeighborBuilderContact
 * @brief A contact neighbor builder functor for contact relation.
 */
class NeighborBuilderContact : public NeighborBuilder,
                               public execution::DeviceExecutable<NeighborBuilderContact, NeighborBuilderContactKernel>
{
  public:
    NeighborBuilderContact(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderContact(){};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);
};

/**
 * @class NeighborBuilderSurfaceContact
 * @brief A solid contact neighbor builder functor when bodies having surface contact.
 */
class NeighborBuilderSurfaceContact : public NeighborBuilderContact
{
  private:
    UniquePtrKeeper<Kernel> kernel_keeper_;

  public:
    NeighborBuilderSurfaceContact(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderSurfaceContact(){};
};

/**
 * @class NeighborBuilderContactBodyPart
 * @brief A contact neighbor builder functor when particles j belongs a body part.
 */
class NeighborBuilderContactBodyPart : public NeighborBuilder
{
  public:
    NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part);
    virtual ~NeighborBuilderContactBodyPart(){};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

  protected:
    StdLargeVec<int> part_indicator_; /**< indicator of the body part */
};

/**
 * @class NeighborBuilderContactAdaptive
 * @brief A contact neighbor builder functor when the particles have different smoothing lengths.
 */
class NeighborBuilderContactAdaptive : public NeighborBuilder
{
  public:
    explicit NeighborBuilderContactAdaptive(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderContactAdaptive(){};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

  protected:
    SPHAdaptation &adaptation_, &contact_adaptation_;
    Real relative_h_ref_;
};

} // namespace SPH
#endif // NEIGHBORHOOD_H
