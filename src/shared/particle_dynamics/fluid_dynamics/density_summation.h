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
 * @file density_summation.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef DENSITY_SUMMATION_INNER_H
#define DENSITY_SUMMATION_INNER_H

#include "base_fluid_dynamics.h"
#include "device_implementation.hpp"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class DensitySummation;

template <typename... InteractionTypes>
class DensitySummationKernel;

template <>
class DensitySummationKernel<Base>
{
  public:
    explicit DensitySummationKernel(BaseParticles* particles, DeviceReal rho0, DeviceReal invSigma0, DeviceReal W0)
        : rho_(particles->getDeviceVariableByName<DeviceReal>("Density")),
          rho_sum_(particles->registerDeviceVariable<DeviceReal>("DensitySummation", particles->total_real_particles_)),
          mass_(particles->getDeviceVariableByName<DeviceReal>("Mass")),
          Vol_(particles->getDeviceVariableByName<DeviceReal>("Volume")),
          rho0_(rho0), inv_sigma0_(invSigma0), W0_(W0) {}

  protected:
    DeviceReal *rho_, *rho_sum_, *mass_, *Vol_;
    DeviceReal rho0_, inv_sigma0_, W0_;
};

template <class DataDelegationType>
class DensitySummation<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit DensitySummation(BaseRelationType &base_relation);
    virtual ~DensitySummation(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_, &rho_sum_, &Vol_;
    Real rho0_, inv_sigma0_, W0_;
};

template <>
class DensitySummationKernel<Inner<Base>> : public DensitySummationKernel<Base>
{
  public:
    template<class ...Args>
    explicit DensitySummationKernel(BaseInnerRelation &inner_relation, Args &&...baseArgs)
        : DensitySummationKernel<Base>(std::forward<Args>(baseArgs)...),
          cell_linked_list_(inner_relation.getInnerCellLinkedListDevice()),
          inner_neighbor_builder_(inner_relation.getInnerNeighborBuilderDevice()) {}

  protected:
    CellLinkedListKernel* cell_linked_list_;
    NeighborBuilderInnerKernel *inner_neighbor_builder_;
};

template <>
class DensitySummation<Inner<Base>> : public DensitySummation<Base, FluidDataInner>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation)
        : DensitySummation<Base, FluidDataInner>(inner_relation){};
    virtual ~DensitySummation(){};
};

template <>
class DensitySummationKernel<Inner<>> : public DensitySummationKernel<Inner<Base>>
{
  public:
    template<class ...Args>
    explicit DensitySummationKernel(Args &&...baseArgs)
        : DensitySummationKernel<Inner<Base>>(std::forward<Args>(baseArgs)...) {}

    void interaction(DeviceInt index_i, DeviceReal dt = 0.0) {
        DeviceReal sigma = W0_;
        const auto &neighbor_builder = *inner_neighbor_builder_;
        cell_linked_list_->forEachInnerNeighbor(index_i, [&](const DeviceVecd &pos_i, size_t index_j, const DeviceVecd &pos_j)
                                                {
                                                    if(neighbor_builder.isWithinCutoff(pos_i, pos_j) && index_i != index_j)
                                                        sigma += neighbor_builder.W_ij(pos_i, pos_j);
                                                });
        rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
    }

    void update(DeviceInt index_i, DeviceReal dt = 0.0) {
        rho_[index_i] = rho_sum_[index_i];
        Vol_[index_i] = mass_[index_i] / rho_[index_i];
    }
};

template <>
class DensitySummation<Inner<>> : public DensitySummation<Inner<Base>>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation)
        : DensitySummation<Inner<Base>>(inner_relation),
          device_kernel(inner_relation, this->particles_, this->rho0_, this->inv_sigma0_, this->W0_){};
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    execution::DeviceImplementation<DensitySummationKernel<Inner<>>> device_kernel;
};
using DensitySummationInner = DensitySummation<Inner<>>;

template <>
class DensitySummation<Inner<Adaptive>> : public DensitySummation<Inner<Base>>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation);
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

template <>
class DensitySummationKernel<Contact<Base>> : public DensitySummationKernel<Base>
{
  public:
    template<class ...Args>
    explicit DensitySummationKernel(const BaseContactRelation& contact_relation,
                                    DeviceReal *contactInvRho0, DeviceReal **contactMass,
                                    Args &&... baseArgs)
        : DensitySummationKernel<Base>(std::forward<Args>(baseArgs)...),
          contact_inv_rho0_(contactInvRho0), contact_mass_(contactMass),
          particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")),
          contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
          contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
          contact_bodies_size_(contact_relation.contact_bodies_.size()) {}

  protected:
    DeviceReal ContactSummation(DeviceInt index_i)
    {
        DeviceReal sigma{0.0};
        for (auto k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceReal* contact_mass_k = contact_mass_[k];
            const DeviceReal& contact_inv_rho0_k = contact_inv_rho0_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, DeviceInt index_j, const DeviceVecd &pos_j)
                                                           {
                                                               if(neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                                   sigma += neighbor_builder.W_ij(pos_i, pos_j) *
                                                                            contact_inv_rho0_k * contact_mass_k[index_j];
                                                           });
        }
        return sigma;
    };

    DeviceReal *contact_inv_rho0_, **contact_mass_;
    DeviceVecd *particles_position_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    const DeviceInt contact_bodies_size_;
};

template <>
class DensitySummation<Contact<Base>> : public DensitySummation<Base, FluidContactData>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    StdSharedVec<DeviceReal> contact_inv_rho0_device_;
    StdSharedVec<DeviceReal*> contact_mass_device_;
    Real ContactSummation(size_t index_i);

    execution::DeviceImplementation<DensitySummationKernel<Contact<Base>>> device_kernel;
};

template <>
class DensitySummationKernel<Contact<>> : public DensitySummationKernel<Contact<Base>>
{
  public:
    template<class ...Args>
    explicit DensitySummationKernel(Args &&...baseArgs)
        : DensitySummationKernel<Contact<Base>>(std::forward<Args>(baseArgs)...) {}

    void interaction(size_t index_i, Real dt = 0.0) {
        DeviceReal sigma = DensitySummationKernel<Contact<Base>>::ContactSummation(index_i);
        rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
    }
};

template <>
class DensitySummation<Contact<>> : public DensitySummation<Contact<Base>>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation)
        : DensitySummation<Contact<Base>>(contact_relation),
          device_kernel(contact_relation, contact_inv_rho0_device_.data(), contact_mass_device_.data(),
                        this->particles_, this->rho0_, this->inv_sigma0_, this->W0_)  {}
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

    execution::DeviceImplementation<DensitySummationKernel<Contact<>>> device_kernel;
};

template <>
class DensitySummation<Contact<Adaptive>> : public DensitySummation<Contact<Base>>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Real> &h_ratio_;
};

template <typename... SummationType>
class DensitySummationKernel<Inner<FreeSurface, SummationType...>> : public DensitySummationKernel<Inner<SummationType...>>
{
  public:
    template<class ...Args>
    explicit DensitySummationKernel(Args &&...baseArgs)
        : DensitySummationKernel<Inner<SummationType...>>(std::forward<Args>(baseArgs)...) {}

    void update(size_t index_i, Real dt)
    {
        this->rho_[index_i] = sycl::fmax(this->rho_sum_[index_i], this->rho0_);
    }
};

template <typename... SummationType>
class DensitySummation<Inner<FreeSurface, SummationType...>> : public DensitySummation<Inner<SummationType...>>
{
  public:
    template <typename... Args>
    explicit DensitySummation(Args &&...args);
    virtual ~DensitySummation(){};
    void update(size_t index_i, Real dt = 0.0);

    execution::DeviceImplementation<DensitySummationKernel<Inner<FreeSurface, SummationType...>>> device_kernel;
};
using DensitySummationFreeSurfaceInner = DensitySummation<Inner<FreeSurface>>;

struct FreeStream
{
    Real operator()(Real rho_sum, Real rho0, Real rho)
    {
        if (rho_sum < rho)
        {
            return rho_sum + SMAX(Real(0), (rho - rho_sum)) * rho0 / rho;
        }
        return rho_sum;
    };
};

struct NotNearSurface
{
    Real operator()(Real rho_sum, Real rho0, Real rho)
    {
        return rho;
    };
};

template <typename NearSurfaceType, typename... SummationType>
class DensitySummation<Inner<NearSurfaceType, SummationType...>>
    : public DensitySummation<Inner<SummationType...>>
{
  public:
    template <typename... Args>
    explicit DensitySummation(Args &&...args);
    virtual ~DensitySummation(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    NearSurfaceType near_surface_rho_;
    StdLargeVec<int> &indicator_;
    bool isNearFreeSurface(size_t index_i);
};
using DensitySummationInnerNotNearSurface = DensitySummation<Inner<NotNearSurface>>;
using DensitySummationInnerFreeStream = DensitySummation<Inner<FreeStream>>;

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseDensitySummationComplex = ComplexInteraction<DensitySummation<InnerInteractionType, ContactInteractionTypes...>>;

using DensitySummationComplex = BaseDensitySummationComplex<Inner<>, Contact<>>;
using DensitySummationComplexAdaptive = BaseDensitySummationComplex<Inner<Adaptive>, Contact<Adaptive>>;
using DensitySummationComplexFreeSurface = BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>>;
using DensitySummationFreeSurfaceComplexAdaptive = BaseDensitySummationComplex<Inner<FreeSurface, Adaptive>, Contact<Adaptive>>;
using DensitySummationFreeStreamComplex = BaseDensitySummationComplex<Inner<FreeStream>, Contact<>>;
using DensitySummationFreeStreamComplexAdaptive = BaseDensitySummationComplex<Inner<FreeStream, Adaptive>, Contact<Adaptive>>;
using DensitySummationNotNearSurfaceComplex = BaseDensitySummationComplex<Inner<NotNearSurface>, Contact<>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H