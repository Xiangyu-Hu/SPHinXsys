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
 * @file 	fluid_integration.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_INTEGRATION_H
#define FLUID_INTEGRATION_H

#include "base_fluid_dynamics.h"
#include "base_local_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
class WeaklyCompressibleFluid;
namespace fluid_dynamics
{
class FluidInitialCondition : public LocalDynamics, public FluidDataSimple
{
  public:
    explicit FluidInitialCondition(SPHBody &sph_body);
    virtual ~FluidInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
};

class ContinuumVolumeUpdate : public LocalDynamics, public FluidDataSimple
{
public:
    explicit ContinuumVolumeUpdate(SPHBody& sph_body);
    virtual ~ContinuumVolumeUpdate() {};

    void update(size_t index_i, Real dt)
    {
        Vol_[index_i] = mass_[index_i] / rho_[index_i];
    }

protected:
    StdLargeVec<Real> &Vol_, &mass_, &rho_;
};

template<class MaterialType>
class BaseIntegrationKernel {
    using MaterialTypeKernel = typename decltype(MaterialType::device_kernel)::KernelType;

  public:
    explicit BaseIntegrationKernel(BaseParticles *particles)
        : fluid_(*DynamicCast<MaterialType>(this, particles->getBaseMaterial()).device_kernel.get_ptr()),
          rho_(particles->getDeviceVariableByName<DeviceReal>("Density")),
          p_(particles->registerDeviceVariable<DeviceReal>("Pressure", particles->total_real_particles_)),
          drho_dt_(particles->registerDeviceVariable<DeviceReal>("DensityChangeRate", particles->total_real_particles_)),
          mass_(particles->getDeviceVariableByName<DeviceReal>("Mass")),
          Vol_(particles->getDeviceVariableByName<DeviceReal>("Volume")),
          pos_(particles->getDeviceVariableByName<DeviceVecd>("Position")),
          vel_(particles->getDeviceVariableByName<DeviceVecd>("Velocity")),
          force_(particles->getDeviceVariableByName<DeviceVecd>("Force")),
          force_prior_(particles->getDeviceVariableByName<DeviceVecd>("ForcePrior")) {}

  protected:
    const MaterialTypeKernel fluid_;
    DeviceReal *rho_, *p_, *drho_dt_, *mass_, *Vol_;
    DeviceVecd *pos_, *vel_, *force_, *force_prior_;
};

template <class DataDelegationType>
class BaseIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegration(BaseRelationType &base_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &mass_, &Vol_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &force_, &force_prior_;
};

template <typename... InteractionTypes>
class Integration1stHalf;

template <typename... InteractionTypes>
class Integration1stHalfKernel;

template<class RiemannSolverType, class MaterialType>
class Integration1stHalfKernel<Inner<>, RiemannSolverType, MaterialType> : public BaseIntegrationKernel<MaterialType>
{
  public:
    Integration1stHalfKernel(BaseInnerRelation& inner_relation, BaseParticles* particles,
                             const RiemannSolverType &riemann_solver)
        : BaseIntegrationKernel<MaterialType>(particles), correction_(particles),
          riemann_solver_(*riemann_solver.device_kernel.get_ptr()),
          cell_linked_list_(inner_relation.getInnerCellLinkedListDevice()),
          inner_neighbor_builder_(inner_relation.getInnerNeighborBuilderDevice()) {}

    void initialization(DeviceInt index_i, DeviceReal dt = 0.0) {
        this->rho_[index_i] += this->drho_dt_[index_i] * dt * DeviceReal(0.5);
        this->p_[index_i] = this->fluid_.getPressure(this->rho_[index_i]);
        this->pos_[index_i] += this->vel_[index_i] * dt * DeviceReal(0.5);
    }

    void update(DeviceInt index_i, DeviceReal dt = 0.0) {
        this->vel_[index_i] += (this->force_prior_[index_i] + this->force_[index_i]) / this->mass_[index_i] * dt;
    }

    void interaction(DeviceInt index_i, DeviceReal dt = 0.0) {
        auto force = VecdZero<DeviceVecd>();
        DeviceReal rho_dissipation(0);
        const auto &pressure_i = this->p_[index_i];
        const auto &mass_i = this->mass_[index_i];
        const auto correction_i = correction_(index_i);
        const auto& neighbor_builder = *inner_neighbor_builder_;
        cell_linked_list_->forEachInnerNeighbor(index_i, [&](const DeviceVecd &pos_i, size_t index_j, const DeviceVecd &pos_j)
                                                {
                                                    if(neighbor_builder.isWithinCutoff(pos_i, pos_j) && index_i != index_j)
                                                    {
                                                        const auto& dW_ijV_j = neighbor_builder.dW_ij(pos_i, pos_j) * this->Vol_[index_j];
                                                        const auto &e_ij = neighbor_builder.e_ij(pos_i, pos_j);

                                                        force -= mass_i * (pressure_i * correction_i + this->p_[index_j] * correction_(index_j)) * dW_ijV_j * e_ij;
                                                        rho_dissipation += this->riemann_solver_.DissipativeUJump(pressure_i - this->p_[index_j]) * dW_ijV_j;
                                                    }
                                                });
        this->force_[index_i] += force * this->Vol_[index_i] / this->mass_[index_i];
        this->drho_dt_[index_i] = rho_dissipation * this->mass_[index_i] / this->Vol_[index_i];
    }

  protected:
    using RiemannSolverTypeKernel = typename decltype(RiemannSolverType::device_kernel)::KernelType;

    const RiemannSolverTypeKernel riemann_solver_;
    const NoKernelCorrection correction_;
    const CellLinkedListKernel* cell_linked_list_;
    const NeighborBuilderInnerKernel* inner_neighbor_builder_;
};


template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegration<FluidDataInner>
{
  public:
    explicit Integration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;

  public:
    execution::DeviceImplementation<
        Integration1stHalfKernel<Inner<>, RiemannSolverType, WeaklyCompressibleFluid>> device_kernel;
};
using Integration1stHalfInnerNoRiemann = Integration1stHalf<Inner<>, NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerRiemann = Integration1stHalf<Inner<>, AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfCorrectionInnerRiemann = Integration1stHalf<Inner<>, AcousticRiemannSolver, LinearGradientCorrection>;

// The following is used to avoid the C3200 error triggered in Visual Studio.
// Please refer: https://developercommunity.visualstudio.com/t/c-invalid-template-argument-for-template-parameter/831128
using BaseIntegrationWithWall = InteractionWithWall<BaseIntegration>;

template<class MaterialType>
class BaseIntegrationWithWallKernel : public BaseIntegrationKernel<MaterialType>
{
  public:
    BaseIntegrationWithWallKernel(BaseParticles *particles, DeviceReal** wall_mass,
                                  DeviceReal** wall_Vol, DeviceVecd** wall_vel_ave,
                                  DeviceVecd** wall_force_ave, DeviceVecd** wall_n)
        : BaseIntegrationKernel<MaterialType>(particles),
          wall_mass_(wall_mass), wall_Vol_(wall_Vol), wall_vel_ave_(wall_vel_ave),
          wall_force_ave_(wall_force_ave), wall_n_(wall_n) {}

  protected:
    DeviceReal** wall_mass_, **wall_Vol_;
    DeviceVecd** wall_vel_ave_, **wall_force_ave_, **wall_n_;
};

template <class RiemannSolverType, class MaterialType>
class Integration1stHalfKernel<Contact<Wall>, RiemannSolverType, MaterialType>
    : public BaseIntegrationWithWallKernel<MaterialType>
{
  public:
    template<class ...BaseArgs>
    Integration1stHalfKernel(const BaseContactRelation& contact_relation, BaseParticles* particles,
                             const RiemannSolverType& riemann_solver, BaseArgs&& ...args)
        : BaseIntegrationWithWallKernel<MaterialType>(particles, std::forward<BaseArgs>(args)...),
          correction_(particles), riemann_solver_(*riemann_solver.device_kernel.get_ptr()),
          contact_bodies_size_(contact_relation.contact_bodies_.size()),
          contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
          contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
          particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")) {}

    void interaction(DeviceInt index_i, DeviceReal dt = 0.0) {
        DeviceVecd force = VecdZero<DeviceVecd>();
        DeviceReal rho_dissipation{0};
        const DeviceVecd &force_prior_i = this->force_prior_[index_i];
        const DeviceReal &mass_i = this->mass_[index_i],
                         &Vol_i = this->Vol_[index_i],
                         &pressure_i{this->p_[index_i]},
                         &rho_i{this->rho_[index_i]};
        const auto correction_i = correction_(index_i);
        for (auto k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceVecd* force_ave_k = this->wall_force_ave_[k];
            const DeviceReal* wall_mass_k = this->wall_mass_[k];
            const DeviceReal* wall_Vol_k = this->wall_Vol_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, DeviceInt index_j, const DeviceVecd &pos_j)
                                                           {
                                                               if (neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                               {
                                                                   const auto e_ij = neighbor_builder.e_ij(pos_i, pos_j);
                                                                   const auto dW_ijV_j = neighbor_builder.dW_ij(pos_i, pos_j) * wall_Vol_k[index_j];
                                                                   const auto r_ij = neighbor_builder.r_ij(pos_i, pos_j);

                                                                   const DeviceReal face_wall_external_acceleration = VecdDot(DeviceVecd(force_prior_i / mass_i - force_ave_k[index_j]/ wall_mass_k[index_j]), DeviceVecd(-e_ij));
                                                                   const auto p_in_wall = pressure_i + mass_i / Vol_i * r_ij * sycl::max(DeviceReal(0), face_wall_external_acceleration);
                                                                   force -= mass_i * (pressure_i + p_in_wall) * correction_i * dW_ijV_j * e_ij;
                                                                   rho_dissipation += this->riemann_solver_.DissipativeUJump(pressure_i - p_in_wall) * dW_ijV_j;
                                                               }
                                                           });
        }
        this->force_[index_i] += force * this->Vol_[index_i] / this->mass_[index_i];
        this->drho_dt_[index_i] += rho_dissipation * this->mass_[index_i] / this->Vol_[index_i];
    }

  protected:
    using RiemannSolverTypeKernel = typename decltype(RiemannSolverType::device_kernel)::KernelType;

    const NoKernelCorrection correction_;
    const RiemannSolverTypeKernel riemann_solver_;
    const DeviceInt contact_bodies_size_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    const DeviceVecd *particles_position_;
};

template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegrationWithWall
{
  public:
    explicit Integration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Integration1stHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;

  public:
    execution::DeviceImplementation<
        Integration1stHalfKernel<Contact<Wall>, RiemannSolverType, WeaklyCompressibleFluid>> device_kernel;
};

template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Contact<>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegration<FluidContactData>
{
  public:
    explicit Integration1stHalf(BaseContactRelation &contact_relation);
    virtual ~Integration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    StdVec<KernelCorrectionType> contact_corrections_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<StdLargeVec<Real> *> contact_p_;
    StdVec<StdLargeVec<Real>*> contact_Vol_;
};

template <class RiemannSolverType, class KernelCorrectionType>
using Integration1stHalfWithWall = ComplexInteraction<Integration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType, KernelCorrectionType>;

using Integration1stHalfWithWallNoRiemann = Integration1stHalfWithWall<NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfWithWallRiemann = Integration1stHalfWithWall<AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfCorrectionWithWallRiemann = Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrection>;

using MultiPhaseIntegration1stHalfWithWallRiemann =
    ComplexInteraction<Integration1stHalf<Inner<>, Contact<>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>;

template <typename... InteractionTypes>
class Integration2ndHalf;

template <typename... InteractionTypes>
class Integration2ndHalfKernel;

template<class RiemannSolverType, class MaterialType>
class Integration2ndHalfKernel<Inner<>, RiemannSolverType, MaterialType> : public BaseIntegrationKernel<MaterialType>
{
  public:
    Integration2ndHalfKernel(BaseInnerRelation& inner_relation, BaseParticles* particles,
                             const RiemannSolverType &riemann_solver)
        : BaseIntegrationKernel<MaterialType>(particles),
          riemann_solver_(*riemann_solver.device_kernel.get_ptr()),
          cell_linked_list_(inner_relation.getInnerCellLinkedListDevice()),
          inner_neighbor_builder_(inner_relation.getInnerNeighborBuilderDevice()) {}

    void initialization(DeviceInt index_i, DeviceReal dt = 0.0) {
        this->pos_[index_i] += this->vel_[index_i] * dt * DeviceReal(0.5);
    }

    void update(DeviceInt index_i, DeviceReal dt = 0.0) {
        this->rho_[index_i] += this->drho_dt_[index_i] * dt * DeviceReal(0.5);
    }

    void interaction(DeviceInt index_i, DeviceReal dt = 0.0)
    {
        DeviceReal density_change_rate(0);
        auto p_dissipation = VecdZero<DeviceVecd>();
        const DeviceVecd &vel_i = this->vel_[index_i];
        const DeviceReal &mass_i = this->mass_[index_i];
        const auto &neighbor_builder = *inner_neighbor_builder_;
        cell_linked_list_->forEachInnerNeighbor(index_i, [&](const DeviceVecd &pos_i, size_t index_j, const DeviceVecd &pos_j)
                                                {
                                                    if(neighbor_builder.isWithinCutoff(pos_i, pos_j) && index_i != index_j)
                                                    {
                                                        const auto &e_ij = neighbor_builder.e_ij(pos_i, pos_j);
                                                        const auto &dW_ijV_j = neighbor_builder.dW_ij(pos_i, pos_j) * this->Vol_[index_j];

                                                        const DeviceReal u_jump = VecdDot(DeviceVecd(vel_i - this->vel_[index_j]), e_ij);
                                                        density_change_rate += u_jump * dW_ijV_j;
                                                        p_dissipation += mass_i * static_cast<DeviceReal>(this->riemann_solver_.DissipativePJump(u_jump)) * dW_ijV_j * e_ij;
                                                    } });
        this->drho_dt_[index_i] += density_change_rate * this->mass_[index_i] / this->Vol_[index_i];
        this->force_[index_i] = p_dissipation * this->Vol_[index_i] / this->mass_[index_i];
    }

  protected:
    using RiemannSolverTypeKernel = typename decltype(RiemannSolverType::device_kernel)::KernelType;

    const RiemannSolverTypeKernel riemann_solver_;
    CellLinkedListKernel* cell_linked_list_;
    NeighborBuilderInnerKernel* inner_neighbor_builder_;
};

template <class RiemannSolverType>
class Integration2ndHalf<Inner<>, RiemannSolverType>
    : public BaseIntegration<FluidDataInner>
{
  public:
    typedef RiemannSolverType RiemannSolver;

    explicit Integration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Integration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Real> &mass_, & Vol_;

  public:
    execution::DeviceImplementation<Integration2ndHalfKernel<Inner<>, RiemannSolverType, WeaklyCompressibleFluid>> device_kernel;
};
using Integration2ndHalfInnerRiemann = Integration2ndHalf<Inner<>, AcousticRiemannSolver>;
using Integration2ndHalfInnerNoRiemann = Integration2ndHalf<Inner<>, NoRiemannSolver>;
using Integration2ndHalfInnerDissipativeRiemann = Integration2ndHalf<Inner<>, DissipativeRiemannSolver>;

template <class RiemannSolverType, class MaterialType>
class Integration2ndHalfKernel<Contact<Wall>, RiemannSolverType, MaterialType>
    : public BaseIntegrationWithWallKernel<MaterialType>
{
  public:
    template<class ...BaseArgs>
    Integration2ndHalfKernel(const BaseContactRelation& contact_relation, BaseParticles* particles,
                             const RiemannSolverType &riemann_solver, BaseArgs&& ...args)
        : BaseIntegrationWithWallKernel<MaterialType>(particles, std::forward<BaseArgs>(args)...),
          riemann_solver_(*riemann_solver.device_kernel.get_ptr()),
          contact_bodies_size_(contact_relation.contact_bodies_.size()),
          contact_cell_linked_lists_(contact_relation.getContactCellLinkedListsDevice()),
          contact_neighbor_builders_(contact_relation.getContactNeighborBuilderDevice()),
          particles_position_(contact_relation.base_particles_.getDeviceVariableByName<DeviceVecd>("Position")) {}

    void interaction(DeviceInt index_i, DeviceReal dt = 0.0) {
        DeviceReal density_change_rate{0};
        auto p_dissipation = VecdZero<DeviceVecd>();
        const DeviceVecd vel_i = this->vel_[index_i];
        const DeviceReal mass_i = this->mass_[index_i];
        for (auto k = 0; k < contact_bodies_size_; ++k)
        {
            const DeviceReal *wall_Vol_k = this->wall_Vol_[k];
            const DeviceVecd *vel_ave_k = this->wall_vel_ave_[k];
            const DeviceVecd *n_k = this->wall_n_[k];
            const auto& neighbor_builder = *contact_neighbor_builders_[k];
            contact_cell_linked_lists_[k]->forEachNeighbor(index_i, particles_position_,
                                                           [&](const DeviceVecd &pos_i, DeviceInt index_j, const DeviceVecd &pos_j)
                                                           {
                                                               if (neighbor_builder.isWithinCutoff(pos_i, pos_j))
                                                               {
                                                                   const auto e_ij = neighbor_builder.e_ij(pos_i, pos_j);
                                                                   const auto dW_ijV_j = neighbor_builder.dW_ij(pos_i, pos_j) * wall_Vol_k[index_j];

                                                                   const DeviceVecd vel_in_wall = static_cast<DeviceReal>(2.0) * vel_ave_k[index_j] - vel_i;
                                                                   density_change_rate += VecdDot(DeviceVecd(vel_i - vel_in_wall), e_ij) * dW_ijV_j;
                                                                   const DeviceReal u_jump = static_cast<DeviceReal>(2.0) * VecdDot(DeviceVecd(vel_i - vel_ave_k[index_j]), n_k[index_j]);
                                                                   p_dissipation += mass_i * static_cast<DeviceReal>(this->riemann_solver_.DissipativePJump(u_jump)) * dW_ijV_j * n_k[index_j];
                                                               }
                                                           });
        }
        this->drho_dt_[index_i] += density_change_rate * this->mass_[index_i] / this->Vol_[index_i];
        this->force_[index_i] += p_dissipation * this->Vol_[index_i] / this->mass_[index_i];
    }

  protected:
    using RiemannSolverTypeKernel = typename decltype(RiemannSolverType::device_kernel)::KernelType;

    const RiemannSolverTypeKernel riemann_solver_;
    const DeviceInt contact_bodies_size_;
    CellLinkedListKernel **contact_cell_linked_lists_;
    NeighborBuilderContactKernel **contact_neighbor_builders_;
    DeviceVecd *particles_position_;
};

template <class RiemannSolverType>
class Integration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit Integration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Integration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;

  public:
    execution::DeviceImplementation<Integration2ndHalfKernel<Contact<Wall>, RiemannSolverType, WeaklyCompressibleFluid>> device_kernel;
};

template <class RiemannSolverType>
class Integration2ndHalf<Contact<>, RiemannSolverType>
    : public BaseIntegration<FluidContactData>
{
  public:
    explicit Integration2ndHalf(BaseContactRelation &contact_relation);
    virtual ~Integration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<StdLargeVec<Real> *> contact_p_, contact_Vol_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

template <class RiemannSolverType>
using Integration2ndHalfWithWall = ComplexInteraction<Integration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;

using Integration2ndHalfWithWallNoRiemann = Integration2ndHalfWithWall<NoRiemannSolver>;
using Integration2ndHalfWithWallRiemann = Integration2ndHalfWithWall<AcousticRiemannSolver>;

using MultiPhaseIntegration2ndHalfWithWallRiemann =
    ComplexInteraction<Integration2ndHalf<Inner<>, Contact<>, Contact<Wall>>, AcousticRiemannSolver>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_INTEGRATION_H