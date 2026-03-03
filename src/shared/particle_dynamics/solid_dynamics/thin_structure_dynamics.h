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
#include "thin_structure_math.h"

namespace SPH
{
namespace thin_structure_dynamics
{
/**
 * @class UpdateShellNormalDirection
 * @brief update particle normal directions for shell
 */
class UpdateShellNormalDirection : public LocalDynamics
{
  protected:
    Vecd *n_;
    Matd *F_;
    Matd *transformation_matrix0_;

  public:
    explicit UpdateShellNormalDirection(SPHBody &sph_body);
    virtual ~UpdateShellNormalDirection() {};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShellAcousticTimeStepSize
 * @brief Computing the acoustic time step size for shell
 */
class ShellAcousticTimeStepSize : public LocalDynamicsReduce<ReduceMin>
{
  protected:
    Real CFL_;
    ElasticSolid &elastic_solid_;
    Vecd *vel_, *force_, *angular_vel_, *dangular_vel_dt_, *force_prior_;
    Real *thickness_, *mass_;
    Real rho0_, E0_, nu_, c0_;
    Real smoothing_length_;

  public:
    explicit ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
    virtual ~ShellAcousticTimeStepSize() {};

    Real reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShellCorrectConfiguration
 * @brief obtain the corrected initial configuration in strong form
 */
class ShellCorrectConfiguration : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ShellCorrectConfiguration(BaseInnerRelation &inner_relation);
    virtual ~ShellCorrectConfiguration() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        /** A small number is added to diagonal to avoid dividing by zero. */
        Matd global_configuration = Eps * Matd::Identity();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            global_configuration += r_ji * gradW_ijV_j.transpose();
        }
        Matd local_configuration =
            transformation_matrix0_[index_i] * global_configuration * transformation_matrix0_[index_i].transpose();
        /** correction matrix is obtained from local configuration. */
        B_[index_i] = getCorrectionMatrix(local_configuration);
    };

  protected:
    Real *Vol_;
    Matd *B_;
    Vecd *n0_;
    Matd *transformation_matrix0_;
};

/**
 * @class ShellDeformationGradientTensor
 * @brief computing deformation gradient tensor for shell
 * TODO: need a test case for this.
 */
class ShellDeformationGradientTensor : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ShellDeformationGradientTensor(BaseInnerRelation &inner_relation);
    virtual ~ShellDeformationGradientTensor() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pseudo_n_i = pseudo_n_[index_i];
        const Vecd &pos_n_i = pos_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix0_[index_i];

        Matd deformation_part_one = Matd::Zero();
        Matd deformation_part_two = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            deformation_part_one -= (pos_n_i - pos_[index_j]) * gradW_ijV_j.transpose();
            deformation_part_two -= ((pseudo_n_i - n0_[index_i]) - (pseudo_n_[index_j] - n0_[index_j])) * gradW_ijV_j.transpose();
        }
        F_[index_i] = transformation_matrix_i * deformation_part_one * transformation_matrix_i.transpose() * B_[index_i];
        F_[index_i].col(Dimensions - 1) = transformation_matrix_i * pseudo_n_[index_i];
        F_bending_[index_i] = transformation_matrix_i * deformation_part_two * transformation_matrix_i.transpose() * B_[index_i];
    };

  protected:
    Real *Vol_;
    Vecd *pos_, *pseudo_n_, *n0_;
    Matd *B_, *F_, *F_bending_;
    Matd *transformation_matrix0_;
};

/**
 * @class BaseShellRelaxation
 * @brief abstract class for preparing shell relaxation
 */
class BaseShellRelaxation : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit BaseShellRelaxation(BaseInnerRelation &inner_relation);
    virtual ~BaseShellRelaxation() {};

  protected:
    Real *thickness_, *Vol_;
    Vecd *pos_, *vel_, *force_, *force_prior_;
    Vecd *n0_, *pseudo_n_, *dpseudo_n_dt_, *dpseudo_n_d2t_, *rotation_,
        *angular_vel_, *dangular_vel_dt_;
    Matd *transformation_matrix0_; // Transformation matrix from global to local coordinates
    Matd *B_, *F_, *dF_dt_, *F_bending_, *dF_bending_dt_;
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
                                            int number_of_gaussian_points = 3, bool hourglass_control = false, Real hourglass_control_factor = 0.002);
    virtual ~ShellStressRelaxationFirstHalf() {};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &global_shear_stress_i = global_shear_stress_[index_i];
        const Matd &global_stress_i = global_stress_[index_i];
        const Matd &global_moment_i = global_moment_[index_i];

        Vecd force = Vecd::Zero();
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
                Vecd pos_jump = getLinearVariableJump(e_ij, r_ij, pos_[index_i], global_F_[index_i], pos_[index_j], global_F_[index_j]);
                Real limiter_pos = SMIN(2.0 * pos_jump.norm() / r_ij, 1.0);
                force += mass_[index_i] * hourglass_control_factor_ * weight * G0_ * pos_jump * Dimensions *
                         inner_neighborhood.dW_ij_[n] * Vol_[index_j] * limiter_pos;

                Vecd pseudo_n_variation_i = pseudo_n_[index_i] - n0_[index_i];
                Vecd pseudo_n_variation_j = pseudo_n_[index_j] - n0_[index_j];
                Vecd pseudo_n_jump = getLinearVariableJump(e_ij, r_ij, pseudo_n_variation_i, global_F_bending_[index_i],
                                                           pseudo_n_variation_j, global_F_bending_[index_j]);
                Real limiter_pseudo_n = SMIN(2.0 * pseudo_n_jump.norm() / ((pseudo_n_variation_i - pseudo_n_variation_j).norm() + Eps), 1.0);
                pseudo_normal_acceleration += hourglass_control_factor_ * weight * G0_ * pseudo_n_jump * Dimensions *
                                              inner_neighborhood.dW_ij_[n] * Vol_[index_j] * pow(thickness_[index_i], 2) * limiter_pseudo_n;
            }

            force += mass_[index_i] * (global_stress_i + global_stress_[index_j]) * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            pseudo_normal_acceleration += (global_moment_i + global_moment_[index_j]) * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        }

        force_[index_i] = force * inv_rho0_ / thickness_[index_i];
        dpseudo_n_d2t_[index_i] = pseudo_normal_acceleration * inv_rho0_ * 12.0 / pow(thickness_[index_i], 3);

        /** the relation between pseudo-normal and rotations */
        Vecd local_dpseudo_n_d2t = transformation_matrix0_[index_i] * dpseudo_n_d2t_[index_i];
        dangular_vel_dt_[index_i] = getRotationFromPseudoNormal(local_dpseudo_n_d2t, rotation_[index_i], angular_vel_[index_i], dt);
    };

    void update(size_t index_i, Real dt = 0.0);

  protected:
    ElasticSolid &elastic_solid_;
    Real rho0_, inv_rho0_;
    Real smoothing_length_;
    Matd numerical_damping_scaling_matrix_;
    Real *rho_, *mass_;
    Matd *global_stress_, *global_moment_, *mid_surface_cauchy_stress_;
    Vecd *global_shear_stress_;
    Matd *global_F_, *global_F_bending_;
    Real E0_, G0_, nu_, hourglass_control_factor_;
    bool hourglass_control_;
    const Real inv_W0_ = 1.0 / getSPHAdaptation().getKernel()->W0(ZeroVecd);
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
 * @brief computing stress relaxation process by Verlet time stepping
 * This is the second step
 */
class ShellStressRelaxationSecondHalf : public BaseShellRelaxation
{
  public:
    explicit ShellStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseShellRelaxation(inner_relation) {};
    virtual ~ShellStressRelaxationSecondHalf() {};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &vel_n_i = vel_[index_i];
        const Vecd &dpseudo_n_dt_i = dpseudo_n_dt_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix0_[index_i];

        Matd deformation_gradient_change_rate_part_one = Matd::Zero();
        Matd deformation_gradient_change_rate_part_two = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
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
class ConstrainShellBodyRegion : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    ConstrainShellBodyRegion(BodyPartByParticle &body_part);
    virtual ~ConstrainShellBodyRegion() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vel_, *angular_vel_;
};

/**@class ConstrainShellBodyRegionAlongAxis
 * @brief The boundary conditions are denoted by SS1 according to the references.
 * The axis must be 0 or 1.
 * Note that the average values for FSI are prescribed also.
 */
class ConstrainShellBodyRegionAlongAxis : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis);
    virtual ~ConstrainShellBodyRegionAlongAxis() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    const int axis_; /**< the axis direction for bounding*/
    Vecd *pos_, *pos0_;
    Vecd *vel_, *force_;
    Vecd *rotation_, *angular_vel_, *dangular_vel_dt_;
    Real *mass_;
};

/**
 * @class ShellInitialCurvature
 * @brief  Compute shell initial curvature
 */
class InitialShellCurvature : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit InitialShellCurvature(BaseInnerRelation &inner_relation);
    void update(size_t index_i, Real);

  private:
    Real *Vol_;
    Vecd *n0_;
    Matd *B_;
    Matd *transformation_matrix0_;
    Vecd *n_;
    Matd *F_;
    Matd *F_bending_;

    Real *k1_; // first principle curvature
    Real *k2_; // second principle curvature

    Matd *dn_0_;
};

/**
 * @class ShellCurvature
 * @brief  Update shell curvature during deformation
 */
class ShellCurvatureUpdate : public LocalDynamics
{
  public:
    explicit ShellCurvatureUpdate(SPHBody &sph_body);
    void update(size_t index_i, Real);

  private:
    Matd *transformation_matrix0_;
    Matd *F_;
    Matd *F_bending_;

    Real *k1_; // first principle curvature
    Real *k2_; // second principle curvature

    Matd *dn_0_;
};

/**
 * @class AverageShellCurvature
 * @brief  Calculate shell curvature using the cut-off radius of contact fluid body
 */
class AverageShellCurvature : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit AverageShellCurvature(BaseInnerRelation &inner_relation);
    void update(size_t index_i, Real);

  private:
    Real *Vol_;
    Vecd *n_;
    Real *k1_ave_; // first principle curvature
    Real *k2_ave_; // second principle curvature
};
} // namespace thin_structure_dynamics
} // namespace SPH
#endif // THIN_STRUCTURE_DYNAMICS_H
