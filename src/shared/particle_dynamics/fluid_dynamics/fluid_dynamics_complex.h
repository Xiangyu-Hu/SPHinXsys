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
 * @file 	fluid_dynamics_complex.h
 * @brief 	Here, we define the algorithm classes for complex fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_COMPLEX_H
#define FLUID_DYNAMICS_COMPLEX_H

#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"
#include "solid_body.h"
#include "solid_particles.h"
#include "device_implementation.hpp"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateContact<BaseParticles, SolidParticles, DataDelegateEmptyBase>
    FluidWallData;
typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
    FluidContactData;
typedef DataDelegateContact<BaseParticles, SolidParticles> FSIContactData;
/**
 * @class InteractionWithWall
 * @brief Base class adding interaction with wall to general relaxation process
 */
template <class BaseIntegrationType>
class InteractionWithWall : public BaseIntegrationType, public FluidWallData
{
  public:
    template <class BaseBodyRelationType, typename... Args>
    InteractionWithWall(BaseContactRelation &wall_contact_relation,
                        BaseBodyRelationType &base_body_relation, Args &&...args);
    template <typename... Args>
    InteractionWithWall(ComplexRelation &fluid_wall_relation, Args &&...args)
        : InteractionWithWall(fluid_wall_relation.getContactRelation(),
                              fluid_wall_relation.getInnerRelation(), std::forward<Args>(args)...) {}
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<Real> wall_inv_rho0_;
    StdSharedVec<DeviceReal> wall_inv_rho0_device_;
    StdVec<StdLargeVec<Real> *> wall_mass_;
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;

    StdSharedVec<DeviceReal*> wall_mass_device_;
    StdSharedVec<DeviceVecd*> wall_vel_ave_device_, wall_acc_ave_device_, wall_n_device_;
};

class BaseDensitySummationComplexKernel {
  public:
    BaseDensitySummationComplexKernel(DeviceReal *contactInvRho0, DeviceReal **contactMass,
                                      NeighborhoodDevice **contactConfiguration,
                                      size_t contactConfigurationSize) :
                                      contact_inv_rho0_(contactInvRho0),
                                      contact_mass_(contactMass),
                                      contact_configuration_(contactConfiguration),
                                      contact_configuration_size_(contactConfigurationSize) {}

    BaseDensitySummationComplexKernel(DeviceReal *contactInvRho0, DeviceReal **contactMass,
                                      const BaseContactRelation& contact_relation) :
                                      contact_inv_rho0_(contactInvRho0),
                                      contact_mass_(contactMass),
                                      particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")),
                                      contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
                                      contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
                                      contact_bodies_size_(contact_relation.contact_bodies_.size()) {}

    template<class RealType, class ContactMassFunc, class ContactConfigFunc>
    static RealType ContactSummation(size_t index_i, std::size_t contact_configuration_size,
                                  const RealType* contact_inv_rho0, ContactMassFunc&& contactMassFunc,
                                  ContactConfigFunc&& contactConfigFunc)
    {
        RealType sigma(0.0);
        for (size_t k = 0; k < contact_configuration_size; ++k)
        {
            const RealType* contact_mass_k = contactMassFunc(k);
            const auto& contact_inv_rho0_k = contact_inv_rho0[k];
            const auto& contact_neighborhood = contactConfigFunc(k, index_i);
            for (size_t n = 0; n != contact_neighborhood.current_size(); ++n)
                sigma += contact_neighborhood.W_ij_[n] * contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
            }
            return sigma;
    }

    DeviceReal ContactSummation(size_t index_i)
    {
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
        return BaseDensitySummationComplexKernel::ContactSummation(index_i, contact_configuration_size_,
                contact_inv_rho0_, [&](auto k){ return contact_mass_[k]; },
                [&](auto k, auto index_i) -> NeighborhoodDevice& { return contact_configuration_[k][index_i]; });
#else
        DeviceReal sigma(0.0);
        for (size_t k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceReal* contact_mass_k = contact_mass_[k];
            const DeviceReal& contact_inv_rho0_k = contact_inv_rho0_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, size_t index_j,
                                                               const DeviceVecd &pos_j, const DeviceReal& Vol_j)
                                                           {
                                                               if(neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                                   sigma += neighbor_builder.W_ij(pos_i, pos_j) *
                                                                            contact_inv_rho0_k * contact_mass_k[index_j];
                                                           });
        }
        return sigma;
#endif
    };

  private:
    DeviceReal *contact_inv_rho0_, **contact_mass_;
    NeighborhoodDevice** contact_configuration_;
    std::size_t contact_configuration_size_;

    DeviceVecd *particles_position_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    size_t contact_bodies_size_;
};

/**
 * @class DensitySummation
 * @brief computing density by summation considering contribution from contact bodies
 */
template <class DensitySummationInnerType>
class BaseDensitySummationComplex
    : public BaseInteractionComplex<DensitySummationInnerType, FluidContactData>
{
  public:
    template <typename... Args>
    explicit BaseDensitySummationComplex(Args &&...args);
    virtual ~BaseDensitySummationComplex(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdSharedVec<DeviceReal> contact_inv_rho0_device_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    StdSharedVec<DeviceReal*> contact_mass_device_;

    Real ContactSummation(size_t index_i);

  public:
    execution::DeviceImplementation<BaseDensitySummationComplexKernel> device_kernel;
};

class DensitySummationComplexKernel : public BaseDensitySummationComplexKernel, public DensitySummationInnerKernel {
  public:
    DensitySummationComplexKernel(const BaseDensitySummationComplexKernel& complexKernel,
                                  const DensitySummationInnerKernel& innerKernel)
        : BaseDensitySummationComplexKernel(complexKernel), DensitySummationInnerKernel(innerKernel) {}

    template<class RealType, class InnerInteractionFunc, class ContactSummationFunc>
    static void interaction(size_t index_i, Real dt, RealType* rho_sum, RealType rho0, RealType inv_sigma0, const RealType* mass,
                            InnerInteractionFunc&& innerInteraction, ContactSummationFunc&& ContactSummation)
    {
        innerInteraction(index_i, dt);
        RealType sigma = ContactSummation(index_i);
        rho_sum[index_i] += sigma * rho0 * rho0 * inv_sigma0 / mass[index_i];
    }

    void interaction(size_t index_i, Real dt)
    {
        interaction(index_i, dt, rho_sum_, rho0_, inv_sigma0_, mass_,
                    [&](auto idx, auto delta) { DensitySummationInnerKernel::interaction(idx, delta); },
                    [&](auto idx) { return BaseDensitySummationComplexKernel::ContactSummation(idx); });
    }

};

/**
* @class DensitySummationComplex
* @brief computing density by summation considering contribution from contact bodies
*/
class DensitySummationComplex
    : public BaseDensitySummationComplex<DensitySummationInner>
{
  public:
    template <typename... Args>
    explicit DensitySummationComplex(Args &&...args)
        : BaseDensitySummationComplex<DensitySummationInner>(std::forward<Args>(args)...),
          device_kernel(*BaseDensitySummationComplex<DensitySummationInner>::device_kernel.get_ptr(),
                       *DensitySummationInner::device_kernel.get_ptr()){};
    virtual ~DensitySummationComplex(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    execution::DeviceImplementation<DensitySummationComplexKernel> device_kernel;
};

/**
 * @class DensitySummationComplexAdaptive
 * @brief computing density by summation considering  contribution from contact bodies
 */
class DensitySummationComplexAdaptive
    : public BaseDensitySummationComplex<DensitySummationInnerAdaptive>
{
  public:
    template <typename... Args>
    explicit DensitySummationComplexAdaptive(Args &&...args)
        : BaseDensitySummationComplex<DensitySummationInnerAdaptive>(std::forward<Args>(args)...){};
    virtual ~DensitySummationComplexAdaptive(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ViscousWithWall
 * @brief  template class viscous acceleration with wall boundary
 */
template <class ViscousAccelerationInnerType>
class BaseViscousAccelerationWithWall : public InteractionWithWall<ViscousAccelerationInnerType>
{
  public:
    template <typename... Args>
    BaseViscousAccelerationWithWall(Args &&...args)
        : InteractionWithWall<ViscousAccelerationInnerType>(std::forward<Args>(args)...){};
    virtual ~BaseViscousAccelerationWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

using ViscousAccelerationWithWall = BaseViscousAccelerationWithWall<ViscousAccelerationInner>;

/**
 * @class TransportVelocityCorrectionComplex
 * @brief  transport velocity correction considering the contribution from contact bodies
 */
class TransportVelocityCorrectionComplex
    : public BaseInteractionComplex<TransportVelocityCorrectionInner, FluidContactData>
{
  public:
    template <typename... Args>
    TransportVelocityCorrectionComplex(Args &&...args)
        : BaseInteractionComplex<TransportVelocityCorrectionInner, FluidContactData>(
              std::forward<Args>(args)...){};
    virtual ~TransportVelocityCorrectionComplex(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class TransportVelocityCorrectionComplexAdaptive
 * @brief  transport velocity correction considering the contribution from contact bodies
 */
class TransportVelocityCorrectionComplexAdaptive
    : public BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactData>
{
  public:
    template <typename... Args>
    TransportVelocityCorrectionComplexAdaptive(Args &&...args)
        : BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactData>(
              std::forward<Args>(args)...){};
    virtual ~TransportVelocityCorrectionComplexAdaptive(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

template <class BaseIntegration1stHalfType>
class BaseIntegration1stHalfWithWallKernel : public BaseIntegration1stHalfType {
  public:
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
    template<class ...BaseArgs>
    BaseIntegration1stHalfWithWallKernel(StdSharedVec<NeighborhoodDevice*> &contact_configuration,
                                         DeviceVecd** wall_acc_ave, BaseArgs&& ...args) :
                                         BaseIntegration1stHalfType(std::forward<BaseArgs>(args)...),
                                         contact_configuration_(contact_configuration),
                                         wall_acc_ave_(wall_acc_ave) {}
#endif

    template<class ...BaseArgs>
    BaseIntegration1stHalfWithWallKernel(const BaseContactRelation& contact_relation,
                                         DeviceVecd** wall_acc_ave, BaseArgs&& ...args)
        : BaseIntegration1stHalfType(std::forward<BaseArgs>(args)...),
          wall_acc_ave_(wall_acc_ave),
          contact_bodies_size_(contact_relation.contact_bodies_.size()),
          contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
          contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
          particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")) {}

    template<class RealType, class VecType, class RiemannSolver, class WallNeighborhoodFunc,
             class NonConservativeAccFunc, class WallAccAveFunc>
    static void interaction(size_t index_i, Real dt, RealType *p, RealType *rho, RealType *drho_dt, VecType *acc,
                            RiemannSolver& riemann_solver, std::size_t contact_configuration_size,
                            NonConservativeAccFunc&& computeNonConservativeAcceleration,
                            WallAccAveFunc&& getWallAccAve, WallNeighborhoodFunc&& getWallNeighborhood) {
        const VecType acc_prior_i = computeNonConservativeAcceleration(index_i);
        VecType acceleration = VecdZero<VecType>();
        RealType rho_dissipation{0}, min_external_acc{0};
        for (size_t k = 0; k < contact_configuration_size; ++k)
        {
            const VecType* acc_ave_k = getWallAccAve(k);
            const auto &wall_neighborhood = getWallNeighborhood(k, index_i);
            for (size_t n = 0; n < wall_neighborhood.current_size(); ++n)
            {
                const auto& index_j = wall_neighborhood.j_[n];
                const auto& e_ij = wall_neighborhood.e_ij_[n];
                const auto& dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
                const auto& r_ij = wall_neighborhood.r_ij_[n];

                const RealType face_wall_external_acceleration = VecdDot(VecType(acc_prior_i - acc_ave_k[index_j]), VecType(-e_ij));
                const auto p_in_wall = p[index_i] + rho[index_i] * r_ij * SMAX(min_external_acc, face_wall_external_acceleration);
                acceleration -= (p[index_i] + p_in_wall) * dW_ijV_j * e_ij;
                rho_dissipation += riemann_solver.DissipativeUJump(p[index_i] - p_in_wall) * dW_ijV_j;
            }
        }
        acc[index_i] += acceleration / rho[index_i];
        drho_dt[index_i] += rho_dissipation * rho[index_i];
    }

    void interaction(size_t index_i, Real dt = 0.0) {
        BaseIntegration1stHalfType::interaction(index_i, dt);

#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
        interaction(index_i, dt, this->p_, this->rho_, this->drho_dt_, this->acc_, this->riemann_solver_,
                    contact_configuration_.size(), [&](auto index_i){ return this->acc_prior_[index_i]; },
                    [&](auto k){ return this->wall_acc_ave_[k]; },
                    [&](auto k, auto index_i) -> const NeighborhoodDevice&
                        { return this->contact_configuration_[k][index_i]; });
#else
        const DeviceVecd acc_prior_i = this->acc_prior_[index_i];
        DeviceVecd acceleration = VecdZero<DeviceVecd>();
        DeviceReal rho_dissipation{0}, min_external_acc{0};
        const DeviceReal pressure_i{this->p_[index_i]}, rho_i{this->rho_[index_i]};
        for (size_t k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceVecd* acc_ave_k = this->wall_acc_ave_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, size_t index_j,
                                                               const DeviceVecd &pos_j, const DeviceReal& Vol_j)
                                                           {
                                                               if (neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                               {
                                                                   const auto e_ij = neighbor_builder.e_ij(pos_i, pos_j);
                                                                   const auto dW_ijV_j = neighbor_builder.dW_ijV_j(pos_i, pos_j, Vol_j);
                                                                   const auto r_ij = neighbor_builder.r_ij(pos_i, pos_j);

                                                                   const DeviceReal face_wall_external_acceleration = VecdDot(DeviceVecd(acc_prior_i - acc_ave_k[index_j]), DeviceVecd(-e_ij));
                                                                   const auto p_in_wall = pressure_i + rho_i * r_ij * SMAX(min_external_acc, face_wall_external_acceleration);
                                                                   acceleration -= (pressure_i + p_in_wall) * dW_ijV_j * e_ij;
                                                                   rho_dissipation += this->riemann_solver_.DissipativeUJump(pressure_i - p_in_wall) * dW_ijV_j;
                                                               }
                                                           });
        }
        this->acc_[index_i] += acceleration / this->rho_[index_i];
        this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
#endif
    }
  private:
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
    StdSharedVec<NeighborhoodDevice*> &contact_configuration_;
#endif
    DeviceVecd** wall_acc_ave_;

    size_t contact_bodies_size_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    DeviceVecd *particles_position_;
};

/**
 * @class BaseIntegration1stHalfWithWall
 * @brief  template class pressure relaxation scheme together with wall boundary
 */
template <class BaseIntegration1stHalfType>
class BaseIntegration1stHalfWithWall : public InteractionWithWall<BaseIntegration1stHalfType>
{
    using BaseIntegration1stHalfTypeKernel = typename decltype(BaseIntegration1stHalfType::device_kernel)::KernelType;

  public:
    template <typename... Args>
    BaseIntegration1stHalfWithWall(BaseContactRelation& contact_relation, BaseInnerRelation& inner_relation, Args &&...args)
        : InteractionWithWall<BaseIntegration1stHalfType>(contact_relation, inner_relation, std::forward<Args>(args)...),
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
          device_kernel(*this->contact_configuration_device_, this->wall_acc_ave_device_.data(),
                        BaseIntegration1stHalfType::particles_,
                        BaseIntegration1stHalfType::inner_configuration_device_ ?
                            BaseIntegration1stHalfType::inner_configuration_device_->data() : nullptr,
                        this->riemann_solver_) {}
#else
          device_kernel(contact_relation, this->wall_acc_ave_device_.data(),
                        inner_relation, BaseIntegration1stHalfType::particles_,
                        this->riemann_solver_) {}
#endif

    template <typename... Args>
    BaseIntegration1stHalfWithWall(ComplexRelation& complex_relation, Args &&...args)
        : InteractionWithWall<BaseIntegration1stHalfType>(complex_relation, std::forward<Args>(args)...),
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
          device_kernel(*this->contact_configuration_device_, this->wall_acc_ave_device_.data(),
                        BaseIntegration1stHalfType::particles_,
                        BaseIntegration1stHalfType::inner_configuration_device_ ?
                            BaseIntegration1stHalfType::inner_configuration_device_->data() : nullptr,
                        this->riemann_solver_) {}
#else
          device_kernel(complex_relation.getContactRelation(),
                        this->wall_acc_ave_device_.data(),
                        complex_relation.getInnerRelation(),
                        BaseIntegration1stHalfType::particles_,
                        this->riemann_solver_) {}
#endif

    virtual ~BaseIntegration1stHalfWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    execution::DeviceImplementation<BaseIntegration1stHalfWithWallKernel<BaseIntegration1stHalfTypeKernel>> device_kernel;

  protected:
    virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
};

using Integration1stHalfWithWall = BaseIntegration1stHalfWithWall<Integration1stHalf>;
using Integration1stHalfRiemannWithWall = BaseIntegration1stHalfWithWall<Integration1stHalfRiemann>;

/**
 * @class BaseExtendIntegration1stHalfWithWall
 * @brief template class for pressure relaxation scheme with wall boundary
 * and considering non-conservative acceleration term and wall penalty to prevent
 * particle penetration.
 */
template <class BaseIntegration1stHalfType>
class BaseExtendIntegration1stHalfWithWall : public BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>
{
  public:
    template <class BaseBodyRelationType, typename... Args>
    BaseExtendIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation,
                                         BaseBodyRelationType &base_body_relation,
                                         Args &&...args, Real penalty_strength = 1.0)
        : BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>(
              wall_contact_relation, base_body_relation, std::forward<Args>(args)...),
          penalty_strength_(penalty_strength)
    {
        this->particles_->registerVariable(non_cnsrv_acc_, "NonConservativeAcceleration");
    };
    template <typename... Args>
    BaseExtendIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation,
                                         Args &&...args, Real penalty_strength = 1.0)
        : BaseExtendIntegration1stHalfWithWall(fluid_wall_relation.getContactRelation(),
                                               fluid_wall_relation.getInnerRelation(),
                                               std::forward<Args>(args)..., penalty_strength){};
    virtual ~BaseExtendIntegration1stHalfWithWall(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real penalty_strength_;
    StdLargeVec<Vecd> non_cnsrv_acc_;

    virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
};

using ExtendIntegration1stHalfRiemannWithWall = BaseExtendIntegration1stHalfWithWall<Integration1stHalfRiemann>;


template <class BaseIntegration2ndHalfType>
class BaseIntegration2ndHalfWithWallKernel : public BaseIntegration2ndHalfType {
  public:
    template<class ...BaseArgs>
    BaseIntegration2ndHalfWithWallKernel(const BaseContactRelation& contact_relation,
                                         DeviceVecd** wall_vel_ave, DeviceVecd** wall_n, BaseArgs&& ...args)
        : BaseIntegration2ndHalfType(std::forward<BaseArgs>(args)...),
          wall_vel_ave_(wall_vel_ave), wall_n_(wall_n),
          contact_bodies_size_(contact_relation.contact_bodies_.size()),
          contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
          contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
          particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")) {}

#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
    template<class ...BaseArgs>
    BaseIntegration2ndHalfWithWallKernel(StdSharedVec<NeighborhoodDevice*> &contact_configuration,
                                         DeviceVecd** wall_vel_ave, DeviceVecd** wall_n, BaseArgs&& ...args) :
            BaseIntegration2ndHalfType(std::forward<BaseArgs>(args)...),
            contact_configuration_(contact_configuration),
            wall_vel_ave_(wall_vel_ave), wall_n_(wall_n) {}
#endif

    template<class RealType, class VecType, class RiemannSolver, class WallNeighborhoodFunc,
             class WallVelAveFunc, class WallNormalFunc>
    static void interaction(size_t index_i, Real dt, RealType *rho, RealType *drho_dt, VecType* vel, VecType *acc,
                            RiemannSolver& riemann_solver, std::size_t contact_configuration_size,
                            WallVelAveFunc&& getWallVelAve, WallNormalFunc&& getWallNormal,
                            WallNeighborhoodFunc&& getWallNeighborhood) {
        RealType density_change_rate{0};
        auto p_dissipation = VecdZero<VecType>();
        for (size_t k = 0; k < contact_configuration_size; ++k)
        {
            VecType *vel_ave_k = getWallVelAve(k);
            VecType *n_k = getWallNormal(k);
            const auto &wall_neighborhood = getWallNeighborhood(k, index_i);
            for (size_t n = 0; n != wall_neighborhood.current_size(); ++n)
            {
                const auto &index_j = wall_neighborhood.j_[n];
                const auto &e_ij = wall_neighborhood.e_ij_[n];
                const auto &dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

                const VecType vel_in_wall = static_cast<RealType>(2.0) * vel_ave_k[index_j] - vel[index_i];
                density_change_rate += VecdDot(VecType(vel[index_i] - vel_in_wall), e_ij) * dW_ijV_j;
                const RealType u_jump = static_cast<RealType>(2.0) * VecdDot(VecType(vel[index_i] - vel_ave_k[index_j]), n_k[index_j]);
                p_dissipation += static_cast<RealType>(riemann_solver.DissipativePJump(u_jump)) * dW_ijV_j * n_k[index_j];
            }
        }
        drho_dt[index_i] += density_change_rate * rho[index_i];
        acc[index_i] += p_dissipation / rho[index_i];
    }

    void interaction(size_t index_i, Real dt = 0.0) {
        BaseIntegration2ndHalfType::interaction(index_i, dt);

#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
        interaction(index_i, dt, this->rho_, this->drho_dt_, this->vel_, this->acc_, this->riemann_solver_,
                    contact_configuration_.size(), [&](auto k){ return this->wall_vel_ave_[k]; },
                    [&](auto k){ return this->wall_n_[k]; },
                    [&](auto k, auto index_i) -> const NeighborhoodDevice&
                        { return this->contact_configuration_[k][index_i]; });
#else
        DeviceReal density_change_rate{0};
        auto p_dissipation = VecdZero<DeviceVecd>();
        const DeviceVecd vel_i = this->vel_[index_i];
        for (size_t k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceVecd *vel_ave_k = this->wall_vel_ave_[k];
            const DeviceVecd *n_k = this->wall_n_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, size_t index_j,
                                                               const DeviceVecd &pos_j, const DeviceReal &Vol_j)
                                                           {
                                                               if (neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                               {
                                                                   const auto e_ij = neighbor_builder.e_ij(pos_i, pos_j);
                                                                   const auto dW_ijV_j = neighbor_builder.dW_ijV_j(pos_i, pos_j, Vol_j);

                                                                   const DeviceVecd vel_in_wall = static_cast<DeviceReal>(2.0) * vel_ave_k[index_j] - vel_i;
                                                                   density_change_rate += VecdDot(DeviceVecd(vel_i - vel_in_wall), e_ij) * dW_ijV_j;
                                                                   const DeviceReal u_jump = static_cast<DeviceReal>(2.0) * VecdDot(DeviceVecd(vel_i - vel_ave_k[index_j]), n_k[index_j]);
                                                                   p_dissipation += static_cast<DeviceReal>(this->riemann_solver_.DissipativePJump(u_jump)) * dW_ijV_j * n_k[index_j];
                                                               }
                                                           });
        }
        this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
        this->acc_[index_i] += p_dissipation / this->rho_[index_i];
#endif
    }
  private:
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
    StdSharedVec<NeighborhoodDevice*> &contact_configuration_;
#endif
    DeviceVecd **wall_vel_ave_, **wall_n_;

    size_t contact_bodies_size_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    DeviceVecd *particles_position_;
};


/**
 * @class BaseIntegration2ndHalfWithWall
 * @brief template density relaxation scheme without using different Riemann solvers.
 * The difference from the free surface version is that no Riemann problem is applied
 */
template <class BaseIntegration2ndHalfType>
class BaseIntegration2ndHalfWithWall : public InteractionWithWall<BaseIntegration2ndHalfType>
{
    using BaseIntegration2ndHalfTypeKernel = typename decltype(BaseIntegration2ndHalfType::device_kernel)::KernelType;

  public:
    template <typename... Args>
    BaseIntegration2ndHalfWithWall(BaseContactRelation& contact_relation, BaseInnerRelation& inner_relation, Args &&...args)
        : InteractionWithWall<BaseIntegration2ndHalfType>(contact_relation, inner_relation, std::forward<Args>(args)...),
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
          device_kernel(*this->contact_configuration_device_, this->wall_vel_ave_device_.data(),
                        this->wall_n_device_.data(), BaseIntegration2ndHalfType::particles_,
                        BaseIntegration2ndHalfType::inner_configuration_device_ ?
                            BaseIntegration2ndHalfType::inner_configuration_device_->data() : nullptr,
                        this->riemann_solver_) {}
#else
          device_kernel(contact_relation, this->wall_vel_ave_device_.data(),
                        this->wall_n_device_.data(), inner_relation,
                        BaseIntegration2ndHalfType::particles_, this->riemann_solver_) {}
#endif

          template <typename... Args>
          BaseIntegration2ndHalfWithWall(ComplexRelation& complex_relation, Args &&...args)
              : InteractionWithWall<BaseIntegration2ndHalfType>(complex_relation, std::forward<Args>(args)...),
#ifdef SPHINXSYS_SYCL_COMPUTE_NEIGHBORHOOD
                device_kernel(*this->contact_configuration_device_, this->wall_vel_ave_device_.data(),
                              this->wall_n_device_.data(), BaseIntegration2ndHalfType::particles_,
                              BaseIntegration2ndHalfType::inner_configuration_device_ ?
                                                                                      BaseIntegration2ndHalfType::inner_configuration_device_->data() : nullptr,
                              this->riemann_solver_) {}
#else
                device_kernel(complex_relation.getContactRelation(), this->wall_vel_ave_device_.data(),
                              this->wall_n_device_.data(),
                              complex_relation.getInnerRelation(),
                              BaseIntegration2ndHalfType::particles_,
                              this->riemann_solver_) {}
#endif

    virtual ~BaseIntegration2ndHalfWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    execution::DeviceImplementation<BaseIntegration2ndHalfWithWallKernel<BaseIntegration2ndHalfTypeKernel>> device_kernel;
};

using Integration2ndHalfWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalf>;
using Integration2ndHalfRiemannWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalfRiemann>;

/**
 * @class Oldroyd_BIntegration1stHalfWithWall
 * @brief  first half of the pressure relaxation scheme using Riemann solver.
 */
class Oldroyd_BIntegration1stHalfWithWall : public BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>
{
  public:
    explicit Oldroyd_BIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
        : BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>(fluid_wall_relation){};

    virtual ~Oldroyd_BIntegration1stHalfWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class Oldroyd_BIntegration2ndHalfWithWall
 * @brief  second half of the pressure relaxation scheme using Riemann solver.
 */
class Oldroyd_BIntegration2ndHalfWithWall : public BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>
{
  public:
    explicit Oldroyd_BIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
        : BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>(fluid_wall_relation){};

    virtual ~Oldroyd_BIntegration2ndHalfWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_COMPLEX_H
