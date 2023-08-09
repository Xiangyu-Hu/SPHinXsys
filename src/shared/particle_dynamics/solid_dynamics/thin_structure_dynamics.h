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
 * @file 	thin_structure_dynamics.h
 * @brief 	Here, we define the algorithm classes for thin structure dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Dong Wu, Chi Zhang and Xiangyu Hu
 */

#ifndef THIN_STRUCTURE_DYNAMICS_H
#define THIN_STRUCTURE_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "thin_structure_math.h"

namespace SPH
{
namespace thin_structure_dynamics
{
typedef DataDelegateSimple<ShellParticles> ShellDataSimple;
typedef DataDelegateInner<ShellParticles> ShellDataInner;

/**
 * @class ShellDynamicsInitialCondition
 * @brief  set initial condition for shell particles
 * This is a abstract class to be override for case specific initial conditions.
 */
class ShellDynamicsInitialCondition : public LocalDynamics, public ShellDataSimple
{
  public:
    explicit ShellDynamicsInitialCondition(SPHBody &sph_body);
    virtual ~ShellDynamicsInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &n0_, &n_, &pseudo_n_, &pos0_;
    StdLargeVec<Matd> &transformation_matrix_;
};

/**
 * @class ShellAcousticTimeStepSize
 * @brief Computing the acoustic time step size for shell
 */
class ShellAcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
                                  public ShellDataSimple
{
  protected:
    Real CFL_;
    StdLargeVec<Vecd> &vel_, &acc_, &angular_vel_, &dangular_vel_dt_, &acc_prior_;
    StdLargeVec<Real> &thickness_;
    Real rho0_, E0_, nu_, c0_;
    Real smoothing_length_;

  public:
    explicit ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
    virtual ~ShellAcousticTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShellCorrectConfiguration
 * @brief obtain the corrected initial configuration in strong form
 */
class ShellCorrectConfiguration : public LocalDynamics, public ShellDataInner
{
  public:
    explicit ShellCorrectConfiguration(BaseInnerRelation &inner_relation);
    virtual ~ShellCorrectConfiguration(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        /** A small number is added to diagonal to avoid dividing by zero. */
        Matd global_configuration = Eps * Matd::Identity();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            global_configuration += r_ji * gradW_ijV_j.transpose();
        }
        Matd local_configuration =
            transformation_matrix_[index_i] * global_configuration * transformation_matrix_[index_i].transpose();
        /** correction matrix is obtained from local configuration. */
        B_[index_i] = getCorrectionMatrix(local_configuration);
    };

  protected:
    StdLargeVec<Matd> &B_;
    StdLargeVec<Vecd> &n0_;
    StdLargeVec<Matd> &transformation_matrix_;
};

/**
 * @class ShellDeformationGradientTensor
 * @brief computing deformation gradient tensor for shell
 * TODO: need a test case for this.
 */
class ShellDeformationGradientTensor : public LocalDynamics, public ShellDataInner
{
  public:
    explicit ShellDeformationGradientTensor(BaseInnerRelation &inner_relation);
    virtual ~ShellDeformationGradientTensor(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pseudo_n_i = pseudo_n_[index_i];
        const Vecd &pos_n_i = pos_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix_[index_i];

        Matd deformation_part_one = Matd::Zero();
        Matd deformation_part_two = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            deformation_part_one -= (pos_n_i - pos_[index_j]) * gradW_ijV_j.transpose();
            deformation_part_two -= ((pseudo_n_i - n0_[index_i]) - (pseudo_n_[index_j] - n0_[index_j])) * gradW_ijV_j.transpose();
        }
        F_[index_i] = transformation_matrix_i * deformation_part_one * transformation_matrix_i.transpose() * B_[index_i];
        F_[index_i].col(Dimensions - 1) = transformation_matrix_i * pseudo_n_[index_i];
        F_bending_[index_i] = transformation_matrix_i * deformation_part_two * transformation_matrix_i.transpose() * B_[index_i];
    };

  protected:
    StdLargeVec<Vecd> &pos_, &pseudo_n_, &n0_;
    StdLargeVec<Matd> &B_, &F_, &F_bending_;
    StdLargeVec<Matd> &transformation_matrix_;
};

/**
 * @class BaseShellRelaxation
 * @brief abstract class for preparing shell relaxation
 */
class BaseShellRelaxation : public LocalDynamics, public ShellDataInner
{
  public:
    explicit BaseShellRelaxation(BaseInnerRelation &inner_relation);
    virtual ~BaseShellRelaxation(){};

  protected:
    StdLargeVec<Real> &rho_, &thickness_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
    StdLargeVec<Vecd> &n0_, &pseudo_n_, &dpseudo_n_dt_, &dpseudo_n_d2t_, &rotation_,
        &angular_vel_, &dangular_vel_dt_;
    StdLargeVec<Matd> &B_, &F_, &dF_dt_, &F_bending_, &dF_bending_dt_;
    StdLargeVec<Matd> &transformation_matrix_;
};

/**
 * @class ShellStressRelaxationFirstHalf
 * @brief computing stress relaxation process by Verlet time stepping
 * This is the first step
 */
class ShellStressRelaxationFirstHalf : public BaseShellRelaxation
{
  public:
    explicit ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                            int number_of_gaussian_points = 3, bool hourglass_control = false);
    virtual ~ShellStressRelaxationFirstHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &global_shear_stress_i = global_shear_stress_[index_i];
        const Matd &global_stress_i = global_stress_[index_i];
        const Matd &global_moment_i = global_moment_[index_i];

        Vecd acceleration = Vecd::Zero();
        Vecd pseudo_normal_acceleration = global_shear_stress_i;
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            if (hourglass_control_)
            {
                Vecd e_ij = inner_neighborhood.e_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
                Vecd pos_jump = getLinearVariableJump(e_ij, r_ij, pos_[index_i],
                                                      transformation_matrix_[index_i].transpose() * F_[index_i] * transformation_matrix_[index_i],
                                                      pos_[index_j],
                                                      transformation_matrix_[index_i].transpose() * F_[index_j] * transformation_matrix_[index_i]);
                Real limiter_pos = SMIN(2.0 * pos_jump.norm() / r_ij, 1.0);
                acceleration += hourglass_control_factor_ * weight * G0_ * pos_jump * Dimensions *
                                inner_neighborhood.dW_ijV_j_[n] * limiter_pos;

                Vecd pseudo_n_variation_i = pseudo_n_[index_i] - n0_[index_i];
                Vecd pseudo_n_variation_j = pseudo_n_[index_j] - n0_[index_j];
                Vecd pseudo_n_jump = getLinearVariableJump(e_ij, r_ij, pseudo_n_variation_i,
                                                           transformation_matrix_[index_i].transpose() * F_bending_[index_i] * transformation_matrix_[index_i],
                                                           pseudo_n_variation_j,
                                                           transformation_matrix_[index_j].transpose() * F_bending_[index_j] * transformation_matrix_[index_j]);
                Real limiter_pseudo_n = SMIN(2.0 * pseudo_n_jump.norm() / ((pseudo_n_variation_i- pseudo_n_variation_j).norm() + Eps), 1.0);
                pseudo_normal_acceleration += hourglass_control_factor_ * weight * G0_ * pseudo_n_jump * Dimensions *
                                              inner_neighborhood.dW_ijV_j_[n] * pow(thickness_[index_i], 2) * limiter_pseudo_n;
            }

            acceleration += (global_stress_i + global_stress_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            pseudo_normal_acceleration += (global_moment_i + global_moment_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }

        acc_[index_i] = acceleration * inv_rho0_ / thickness_[index_i];
        dpseudo_n_d2t_[index_i] = pseudo_normal_acceleration * inv_rho0_ * 12.0 / pow(thickness_[index_i], 3);

        /** the relation between pseudo-normal and rotations */
        Vecd local_dpseudo_n_d2t = transformation_matrix_[index_i] * dpseudo_n_d2t_[index_i];
        dangular_vel_dt_[index_i] = getRotationFromPseudoNormalForFiniteDeformation(local_dpseudo_n_d2t, rotation_[index_i], angular_vel_[index_i], dt);
    };

    void update(size_t index_i, Real dt = 0.0);

  protected:
    ElasticSolid &elastic_solid_;
    StdLargeVec<Matd> &global_stress_, &global_moment_, &mid_surface_cauchy_stress_, &numerical_damping_scaling_;
    StdLargeVec<Vecd> &global_shear_stress_, &n_;
    Real rho0_, inv_rho0_;
    Real smoothing_length_, E0_, G0_, nu_, hourglass_control_factor_;
    bool hourglass_control_;
    const Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);
    const Real shear_correction_factor_ = 5.0 / 6.0;

    const StdVec<Real> one_gaussian_point_ = {0.0};
    const StdVec<Real> one_gaussian_weight_ = {2.0};
    const StdVec<Real> three_gaussian_points_ = {0.0, 0.7745966692414834, -0.7745966692414834};
    const StdVec<Real> three_gaussian_weights_ = {0.8888888888888889, 0.5555555555555556, 0.5555555555555556};
    const StdVec<Real> five_gaussian_points_ = {0.0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640};
    const StdVec<Real> five_gaussian_weights_ = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891};
    int number_of_gaussian_points_;
    StdVec<Real> gaussian_point_;
    StdVec<Real> gaussian_weight_;
};

/**
 * @class ShellStressRelaxationSecondHalf
 * @brief computing stress relaxation process by verlet time stepping
 * This is the second step
 */
class ShellStressRelaxationSecondHalf : public BaseShellRelaxation
{
  public:
    explicit ShellStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseShellRelaxation(inner_relation){};
    virtual ~ShellStressRelaxationSecondHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &vel_n_i = vel_[index_i];
        const Vecd &dpseudo_n_dt_i = dpseudo_n_dt_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix_[index_i];

        Matd deformation_gradient_change_rate_part_one = Matd::Zero();
        Matd deformation_gradient_change_rate_part_two = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            deformation_gradient_change_rate_part_one -= (vel_n_i - vel_[index_j]) * gradW_ijV_j.transpose();
            deformation_gradient_change_rate_part_two -= (dpseudo_n_dt_i - dpseudo_n_dt_[index_j]) * gradW_ijV_j.transpose();
        }
        dF_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_one * transformation_matrix_i.transpose() * B_[index_i];
        dF_dt_[index_i].col(Dimensions - 1) = transformation_matrix_i * dpseudo_n_dt_[index_i];
        dF_bending_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_two * transformation_matrix_i.transpose() * B_[index_i];
    };

    void update(size_t index_i, Real dt = 0.0);
};

/**@class ConstrainShellBodyRegion
 * @brief Fix the position and angle of a shell body part.
 */
class ConstrainShellBodyRegion : public BaseLocalDynamics<BodyPartByParticle>, public ShellDataSimple
{
  public:
    ConstrainShellBodyRegion(BodyPartByParticle &body_part);
    virtual ~ConstrainShellBodyRegion(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_, &angular_vel_;
};

/**@class ConstrainShellBodyRegionAlongAxis
 * @brief The boundary conditions are denoted by SS1 according to the references.
 * The axis must be 0 or 1.
 * Note that the average values for FSI are prescribed also.
 */
class ConstrainShellBodyRegionAlongAxis : public BaseLocalDynamics<BodyPartByParticle>, public ShellDataSimple
{
  public:
    ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis);
    virtual ~ConstrainShellBodyRegionAlongAxis(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    const int axis_; /**< the axis direction for bounding*/
    StdLargeVec<Vecd> &pos_, &pos0_;
    StdLargeVec<Vecd> &vel_, &acc_;
    StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
};

/**
 * @class DistributingPointForcesToShell
 * @brief Distribute a series of point forces to its contact shell bodies.
 */
class DistributingPointForcesToShell : public LocalDynamics, public ShellDataSimple
{
  protected:
    std::vector<Vecd> point_forces_, reference_positions_, time_dependent_point_forces_;
    Real time_to_full_external_force_;
    Real particle_spacing_ref_, h_spacing_ratio_;
    StdLargeVec<Vecd> &pos0_, &acc_prior_;
    StdLargeVec<Real> &thickness_;
    std::vector<StdLargeVec<Real>> weight_;
    std::vector<Real> sum_of_weight_;

    void getWeight();

  public:
    DistributingPointForcesToShell(SPHBody &sph_body, std::vector<Vecd> point_forces,
                                   std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                                   Real particle_spacing_ref, Real h_spacing_ratio = 1.6);
    virtual ~DistributingPointForcesToShell(){};

    virtual void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace thin_structure_dynamics
} // namespace SPH
#endif // THIN_STRUCTURE_DYNAMICS_H
