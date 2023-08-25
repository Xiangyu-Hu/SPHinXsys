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
 * @file 	fluid_structure_interaction.h
 * @brief 	Here, we define the algorithm classes for fluid structure interaction.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_STRUCTURE_INTERACTION_H
#define FLUID_STRUCTURE_INTERACTION_H

#include "all_particle_dynamics.h"
#include "base_material.h"
#include "elastic_dynamics.h"
#include "fluid_dynamics_complex.h"
#include "riemann_solver.h"

namespace SPH
{
namespace solid_dynamics
{
typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
typedef DataDelegateContact<SolidParticles, BaseParticles> FSIContactData;

/**
 * @class BaseForceFromFluid
 * @brief Base class for computing the forces from the fluid
 */
class BaseForceFromFluid : public LocalDynamics, public FSIContactData
{
  public:
    explicit BaseForceFromFluid(BaseContactRelation &contact_relation);
    virtual ~BaseForceFromFluid(){};
    StdLargeVec<Vecd> &getForceFromFluid() { return force_from_fluid_; };

  protected:
    StdLargeVec<Real> &Vol_;
    StdVec<Fluid *> contact_fluids_;
    StdLargeVec<Vecd> force_from_fluid_;
};

/**
 * @class ViscousForceFromFluid
 * @brief Computing the viscous force from the fluid
 */
class ViscousForceFromFluid : public BaseForceFromFluid
{
  public:
    explicit ViscousForceFromFluid(BaseContactRelation &contact_relation);
    virtual ~ViscousForceFromFluid(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real Vol_i = Vol_[index_i];
        const Vecd &vel_ave_i = vel_ave_[index_i];

        Vecd force = Vecd::Zero();
        /** Contact interaction. */
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real mu_k = mu_[k];
            Real smoothing_length_k = smoothing_length_[k];
            StdLargeVec<Vecd> &vel_n_k = *(contact_vel_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];

                Vecd vel_derivative = 2.0 * (vel_ave_i - vel_n_k[index_j]) /
                                      (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);

                force += 2.0 * mu_k * vel_derivative * Vol_i * contact_neighborhood.dW_ijV_j_[n];
            }
        }

        force_from_fluid_[index_i] = force;
    };

  protected:
    StdLargeVec<Vecd> &vel_ave_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
    StdVec<Real> mu_;
    StdVec<Real> smoothing_length_;
};

/**
 * @class BasePressureForceAccelerationFromFluid
 * @brief Template class fro computing the pressure force from the fluid with different Riemann solvers.
 * The pressure force is added on the viscous force of the latter is computed.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
template <class RiemannSolverType>
class BasePressureForceAccelerationFromFluid : public BaseForceFromFluid
{
  public:
    explicit BasePressureForceAccelerationFromFluid(BaseContactRelation &contact_relation)
        : BasePressureForceAccelerationFromFluid(true, contact_relation)
    {
        particles_->registerVariable(force_from_fluid_, "PressureForceFromFluid");
    };
    virtual ~BasePressureForceAccelerationFromFluid(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
            StdLargeVec<Real> &p_k = *(contact_p_[k]);
            StdLargeVec<Vecd> &vel_k = *(contact_vel_[k]);
            StdLargeVec<Vecd> &acc_prior_k = *(contact_acc_prior_[k]);
            RiemannSolverType &riemann_solvers_k = riemann_solvers_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Real r_ij = contact_neighborhood.r_ij_[n];
                Real face_wall_external_acceleration = (acc_prior_k[index_j] - acc_ave_[index_i]).dot(e_ij);
                Real p_in_wall = p_k[index_j] + rho_n_k[index_j] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
                Real u_jump = 2.0 * (vel_k[index_j] - vel_ave_[index_i]).dot(n_[index_i]);
                force += (riemann_solvers_k.DissipativePJump(u_jump) * n_[index_i] - (p_in_wall + p_k[index_j]) * e_ij) * Vol_[index_i] * contact_neighborhood.dW_ijV_j_[n];
            }
        }
        force_from_fluid_[index_i] = force;
        acc_prior_[index_i] = force / particles_->ParticleMass(index_i); // TODO: to add gravity contribution
    };

  protected:
    StdLargeVec<Vecd> &vel_ave_, &acc_prior_, &acc_ave_, &n_;
    StdVec<StdLargeVec<Real> *> contact_rho_n_, contact_p_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_, contact_acc_prior_;
    StdVec<RiemannSolverType> riemann_solvers_;

    BasePressureForceAccelerationFromFluid(bool mostDerived, BaseContactRelation &contact_relation)
        : BaseForceFromFluid(contact_relation),
          vel_ave_(*particles_->AverageVelocity()),
          acc_prior_(particles_->acc_prior_),
          acc_ave_(*particles_->AverageAcceleration()), n_(particles_->n_)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_rho_n_.push_back(&(contact_particles_[k]->rho_));
            contact_vel_.push_back(&(contact_particles_[k]->vel_));
            contact_p_.push_back(contact_particles_[k]->template getVariableByName<Real>("Pressure"));
            contact_acc_prior_.push_back(&(contact_particles_[k]->acc_prior_));
            riemann_solvers_.push_back(RiemannSolverType(*contact_fluids_[k], *contact_fluids_[k]));
        }
    };
};
using PressureForceAccelerationFromFluid = BasePressureForceAccelerationFromFluid<NoRiemannSolver>;
using PressureForceAccelerationFromFluidRiemann = BasePressureForceAccelerationFromFluid<AcousticRiemannSolver>;

/**
 * @class BaseAllForceAccelerationFromFluid
 * @brief template class for computing force from fluid with updated viscous force
 */
template <class PressureForceType>
class BaseAllForceAccelerationFromFluid : public PressureForceType
{
  public:
    template <class ViscousForceFromFluidType>
    BaseAllForceAccelerationFromFluid(BaseContactRelation &contact_relation,
                                      ViscousForceFromFluidType &viscous_force_from_fluid)
        : PressureForceType(false, contact_relation),
          viscous_force_from_fluid_(viscous_force_from_fluid.getForceFromFluid())
    {
        this->particles_->registerVariable(this->force_from_fluid_, "AllForceFromFluid");
    };
    virtual ~BaseAllForceAccelerationFromFluid(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        PressureForceType::interaction(index_i, dt);
        this->force_from_fluid_[index_i] += viscous_force_from_fluid_[index_i];
        this->acc_prior_[index_i] += viscous_force_from_fluid_[index_i] / this->particles_->ParticleMass(index_i);
    };

  protected:
    StdLargeVec<Vecd> &viscous_force_from_fluid_;
};
using AllForceAccelerationFromFluid =
    BaseAllForceAccelerationFromFluid<PressureForceAccelerationFromFluid>;
using AllForceAccelerationFromFluidRiemann =
    BaseAllForceAccelerationFromFluid<PressureForceAccelerationFromFluidRiemann>;

/**
 * @class TotalForceFromFluid
 * @brief Computing the total force from fluid
 */
class TotalForceFromFluid : public LocalDynamicsReduce<Vecd, ReduceSum<Vecd>>
{
  protected:
    BaseDynamics<void> &force_from_fluid_dynamics_;
    StdLargeVec<Vecd> &force_from_fluid_;

  public:
    template <class ForceFromFluidDynamicsType>
    explicit TotalForceFromFluid(ForceFromFluidDynamicsType &force_from_fluid_dynamics, const std::string &force_name)
        : LocalDynamicsReduce<Vecd, ReduceSum<Vecd>>(force_from_fluid_dynamics.getSPHBody(), Vecd::Zero()),
          force_from_fluid_dynamics_(force_from_fluid_dynamics),
          force_from_fluid_(force_from_fluid_dynamics.getForceFromFluid())
    {
        quantity_name_ = force_name;
    };

    virtual ~TotalForceFromFluid(){};
    virtual void setupDynamics(Real dt = 0.0) override;
    Vecd reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class InitializeDisplacement
 * @brief initialize the displacement for computing average velocity.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
class InitializeDisplacement : public LocalDynamics, public ElasticSolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_temp_, &pos_;

  public:
    explicit InitializeDisplacement(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp);
    virtual ~InitializeDisplacement(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class UpdateAverageVelocityAndAcceleration
 * @brief Computing average velocity.
 * This class is for FSI applications to achieve smaller solid dynamics
 * time step size compared to the fluid dynamics
 */
class UpdateAverageVelocityAndAcceleration : public LocalDynamics, public ElasticSolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_temp_, &pos_, &vel_ave_, &acc_ave_;

  public:
    explicit UpdateAverageVelocityAndAcceleration(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp);
    virtual ~UpdateAverageVelocityAndAcceleration(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class AverageVelocityAndAcceleration
 * @brief Impose force matching between fluid and solid dynamics.
 * Note that the fluid time step should be larger than that of solid time step.
 * Otherwise numerical instability may occur.
 */
class AverageVelocityAndAcceleration
{
  protected:
    StdLargeVec<Vecd> pos_temp_;

  public:
    SimpleDynamics<InitializeDisplacement> initialize_displacement_;
    SimpleDynamics<UpdateAverageVelocityAndAcceleration> update_averages_;

    explicit AverageVelocityAndAcceleration(SolidBody &solid_body);
    ~AverageVelocityAndAcceleration(){};
};
} // namespace solid_dynamics
} // namespace SPH
#endif // FLUID_STRUCTURE_INTERACTION_H