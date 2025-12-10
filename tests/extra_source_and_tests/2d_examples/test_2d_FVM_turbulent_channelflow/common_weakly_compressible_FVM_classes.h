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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	common_weakly_compressible_FVM_classes.h
 * @brief 	Here, we define the common weakly compressible classes for fluid dynamics in FVM.
 * @author	Zhentong Wang and Xiangyu Hu
 */

#ifndef COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H
#define COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H

#include "sphinxsys.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class WCAcousticTimeStepSizeInFVM
 * @brief Computing the acoustic time step size
 */
class WCAcousticTimeStepSizeInFVM : public fluid_dynamics::AcousticTimeStep
{
  protected:
    Real *rho_, *p_;
    Vecd *vel_;
    Fluid &fluid_;
    Real min_distance_between_nodes_;

  public:
    explicit WCAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real min_distance_between_nodes, Real acousticCFL = 0.6);
    virtual ~WCAcousticTimeStepSizeInFVM() {};
    virtual Real outputResult(Real reduced_value) override;
    Real acousticCFL_;
};

/**
 * @class BaseForceFromFluidInFVM
 * @brief Base class for computing the forces from the fluid.
 * Note that In FVM , we need DataDelegateInner class to calculate force between solid and fluid.
 */
class BaseForceFromFluidInFVM : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit BaseForceFromFluidInFVM(BaseInnerRelation &inner_relation);
    virtual ~BaseForceFromFluidInFVM() {};
    Vecd *getForceFromFluid() { return force_from_fluid_; };

  protected:
    Real *Vol_;
    Vecd *force_from_fluid_;
};

/**
 * @class ViscousForceFromFluidInFVM
 * @brief Computing the viscous force from the fluid
 */
class ViscousForceFromFluidInFVM : public BaseForceFromFluidInFVM
{
  public:
    explicit ViscousForceFromFluidInFVM(BaseInnerRelation &inner_relation, StdVec<StdVec<size_t>> each_boundary_type_contact_real_index);
    virtual ~ViscousForceFromFluidInFVM() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Viscosity &viscosity_;
    Vecd *vel_;
    Real mu_;
    StdVec<StdVec<size_t>> each_boundary_type_contact_real_index_;
};

/**
 * @class PressureForceFromFluidInFVM
 * @brief Template class fro computing the pressure force from the fluid with different Riemann solvers in FVM.
 * The pressure force is added on the viscous force of the latter is computed.
 * time step size compared to the fluid dynamics
 */
template <class EulerianIntegration2ndHalfType>
class PressureForceFromFluidInFVM : public BaseForceFromFluidInFVM
{
    using RiemannSolverType = typename EulerianIntegration2ndHalfType::RiemannSolver;

  public:
    explicit PressureForceFromFluidInFVM(BaseInnerRelation &inner_relation, StdVec<StdVec<size_t>> each_boundary_type_contact_real_index)
        : BaseForceFromFluidInFVM(inner_relation),
          fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          p_(particles_->getVariableDataByName<Real>("Pressure")),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          riemann_solver_(fluid_, fluid_),
          each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index)
    {
        force_from_fluid_ = particles_->registerStateVariableData<Vecd>("PressureForceOnSolid");
    };
    Fluid &fluid_;
    Vecd *vel_;
    Real *p_, *rho_;
    RiemannSolverType riemann_solver_;
    StdVec<StdVec<size_t>> each_boundary_type_contact_real_index_;
    virtual ~PressureForceFromFluidInFVM() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        for (size_t real_particle_num = 0; real_particle_num != each_boundary_type_contact_real_index_[3].size(); ++real_particle_num)
        {
            Vecd force = Vecd::Zero();
            if (index_i == each_boundary_type_contact_real_index_[3][real_particle_num])
            {
                Real Vol_i = Vol_[index_i];
                FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
                const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                size_t index_j = inner_neighborhood.j_[2];
                Vecd e_ij = inner_neighborhood.e_ij_[2];
                FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, e_ij);
                force -= 2.0 * (-e_ij) * interface_state.p_ * Vol_i * inner_neighborhood.dW_ij_[2] * Vol_[index_j];
                force_from_fluid_[index_i] = force;
            }
        }
    };
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // COMMON_WEAKLY_COMPRESSIBLE_FVM_CLASSES_H