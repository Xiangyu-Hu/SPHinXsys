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
 * @file 	slender_structure_dynamics.h
 * @brief 	Here, we define the algorithm classes for slender structure dynamics based on thin structure dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Xipeng Lyu
 */

#ifndef SLENDER_STRUCTURE_DYNAMICS_H
#define SLENDER_STRUCTURE_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "beam_particles.h"
#include "elastic_solid.h"
#include "slender_structure_math.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
namespace slender_structure_dynamics
{
typedef DataDelegateSimple<BarParticles> BarDataSimple;
typedef DataDelegateInner<BarParticles> BarDataInner;

/**
 * @class BarDynamicsInitialCondition
 * @brief  set initial condition for bar particles
 * This is a abstract class to be override for case specific initial conditions.
 */
class BarDynamicsInitialCondition : public LocalDynamics, public BarDataSimple
{
  public:
    explicit BarDynamicsInitialCondition(SPHBody &sph_body);
    virtual ~BarDynamicsInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &n0_, &n_, &pseudo_n_, &pos0_;
    StdLargeVec<Matd> &transformation_matrix_;
    StdLargeVec<Vecd> &b_n0_, &b_n_, &pseudo_b_n_;
};

/**
 * @class BarAcousticTimeStepSize
 * @brief Computing the acoustic time step size for bar
 */
class BarAcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
                                public BarDataSimple
{
  protected:
    Real CFL_;
    StdLargeVec<Vecd> &vel_, &acc_, &angular_vel_, &dangular_vel_dt_, &acc_prior_;
    StdLargeVec<Real> &thickness_;
    Real rho0_, E0_, nu_, c0_;
    Real smoothing_length_;

    StdLargeVec<Vecd> &angular_b_vel_, &dangular_b_vel_dt_;
    StdLargeVec<Real> &width_;

  public:
    explicit BarAcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
    virtual ~BarAcousticTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class BarCorrectConfiguration
 * @brief obtain the corrected initial configuration in strong form
 */
class BarCorrectConfiguration : public LocalDynamics, public BarDataInner
{
  public:
    explicit BarCorrectConfiguration(BaseInnerRelation &inner_relation);
    virtual ~BarCorrectConfiguration(){};

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
        B_[index_i] = getCorrectionMatrix_beam(local_configuration);
    };

  protected:
    StdLargeVec<Matd> &B_;
    StdLargeVec<Vecd> &n0_;
    StdLargeVec<Vecd> &b_n0_;
    StdLargeVec<Matd> &transformation_matrix_;
};

/**
 * @class BarDeformationGradientTensor
 * @brief computing deformation gradient tensor for bar
 * TODO: need a test case for this.
 */
class BarDeformationGradientTensor : public LocalDynamics, public BarDataInner
{
  public:
    explicit BarDeformationGradientTensor(BaseInnerRelation &inner_relation);
    virtual ~BarDeformationGradientTensor(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pseudo_n_i = pseudo_n_[index_i];
        const Vecd &pseudo_b_n_i = pseudo_b_n_[index_i];
        const Vecd &pos_n_i = pos_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix_[index_i];

        Matd deformation_part_one = Matd::Zero();
        Matd deformation_part_two = Matd::Zero();
        Matd deformation_part_three = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            deformation_part_one -= (pos_n_i - pos_[index_j]) * gradW_ijV_j.transpose();
            deformation_part_two -= ((pseudo_n_i - n0_[index_i]) - (pseudo_n_[index_j] - n0_[index_j])) * gradW_ijV_j.transpose();
            deformation_part_three -= ((pseudo_b_n_i - b_n0_[index_i]) - (pseudo_b_n_[index_j] - b_n0_[index_j])) * gradW_ijV_j.transpose();
        }
        F_[index_i] = transformation_matrix_i * deformation_part_one * transformation_matrix_i.transpose() * B_[index_i];
        F_[index_i].col(Dimensions - 1) = transformation_matrix_i * pseudo_n_[index_i];
        F_[index_i].col(Dimensions - 2) = transformation_matrix_i * pseudo_b_n_[index_i];
        F_bending_[index_i] = transformation_matrix_i * deformation_part_two * transformation_matrix_i.transpose() * B_[index_i];
        F_b_bending_[index_i] = transformation_matrix_i * deformation_part_three * transformation_matrix_i.transpose() * B_[index_i];
    };

  protected:
    StdLargeVec<Vecd> &pos_, &pseudo_n_, &n0_;
    StdLargeVec<Matd> &B_, &F_, &F_bending_;
    StdLargeVec<Matd> &transformation_matrix_;

    StdLargeVec<Vecd> &pseudo_b_n_, &b_n0_;
    StdLargeVec<Matd> &F_b_bending_;
};

/**
 * @class BaseBarRelaxation
 * @brief abstract class for preparing bar relaxation
 */
class BaseBarRelaxation : public LocalDynamics, public BarDataInner
{
  public:
    explicit BaseBarRelaxation(BaseInnerRelation &inner_relation);
    virtual ~BaseBarRelaxation(){};

  protected:
    StdLargeVec<Real> &rho_, &thickness_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
    StdLargeVec<Vecd> &n0_, &pseudo_n_, &dpseudo_n_dt_, &dpseudo_n_d2t_, &rotation_,
        &angular_vel_, &dangular_vel_dt_;
    StdLargeVec<Matd> &B_, &F_, &dF_dt_, &F_bending_, &dF_bending_dt_;
    StdLargeVec<Matd> &transformation_matrix_;

    StdLargeVec<Real> &width_;
    StdLargeVec<Vecd> &b_n0_, &pseudo_b_n_, &dpseudo_b_n_dt_, &dpseudo_b_n_d2t_, &rotation_b_,
        &angular_b_vel_, dangular_b_vel_dt_;
    StdLargeVec<Matd> &F_b_bending_, &dF_b_bending_dt_;
};

/**
 * @class BarStressRelaxationFirstHalf
 * @brief computing stress relaxation process by Verlet time stepping
 * This is the first step
 */
class BarStressRelaxationFirstHalf : public BaseBarRelaxation
{
  public:
    explicit BarStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                          int number_of_gaussian_points = 4, bool hourglass_control = false);
    virtual ~BarStressRelaxationFirstHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &global_shear_stress_i = global_shear_stress_[index_i];
        const Matd &global_stress_i = global_stress_[index_i];
        const Matd &global_moment_i = global_moment_[index_i];
        const Matd &global_b_moment_i = global_b_moment_[index_i];
        const Vecd &global_b_shear_stress_i = global_b_shear_stress_[index_i];

        Vecd acceleration = Vecd::Zero();
        Vecd pseudo_normal_acceleration = global_shear_stress_i;
        Vecd pseudo_b_normal_acceleration = global_b_shear_stress_i;

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            acceleration += (global_stress_i + global_stress_[index_j]) *
                            inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            pseudo_normal_acceleration += (global_moment_i + global_moment_[index_j]) *
                                          inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            pseudo_b_normal_acceleration += (global_b_moment_i + global_b_moment_[index_j]) *
                                            inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }

        acc_[index_i] = acceleration * inv_rho0_ / (thickness_[index_i] * width_[index_i]);
        dpseudo_n_d2t_[index_i] = pseudo_normal_acceleration * inv_rho0_ * 12.0 / pow(thickness_[index_i], 4);
        dpseudo_b_n_d2t_[index_i] = -pseudo_b_normal_acceleration * inv_rho0_ * 12.0 / pow(thickness_[index_i], 4);

        Vecd local_dpseudo_n_d2t = transformation_matrix_[index_i] * dpseudo_n_d2t_[index_i];
        Vecd local_dpseudo_b_n_d2t = transformation_matrix_[index_i] * dpseudo_b_n_d2t_[index_i];
        dangular_b_vel_dt_[index_i] = getRotationFromPseudoNormalForSmallDeformation_b(
            Vec3d(local_dpseudo_b_n_d2t), Vec3d(local_dpseudo_n_d2t), Vec3d(rotation_b_[index_i]), Vec3d(angular_b_vel_[index_i]), dt);
        dangular_vel_dt_[index_i] = getRotationFromPseudoNormalForSmallDeformation(
            Vec3d(local_dpseudo_b_n_d2t), Vec3d(local_dpseudo_n_d2t), Vec3d(rotation_[index_i]), Vec3d(angular_vel_[index_i]), dt);
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

    Real gpt = sqrt(3.0 / 5.0);
    const StdVec<Real> nine_gaussian_points_2d_vector_x{-gpt, gpt, gpt, -gpt, 0, gpt, 0, -gpt, 0};
    const StdVec<Real> nine_gaussian_points_2d_vector_y{-gpt, -gpt, gpt, gpt, -gpt, 0, gpt, 0, 0};

    const StdVec<Real> nine_gaussian_weights_2d_ =
        {25. / 81., 25. / 81., 25. / 81., 25. / 81., 40. / 81., 40. / 81., 40. / 81., 40. / 81., 64. / 81.};

    Real gpt_4 = sqrt(1. / 3.);
    const StdVec<Real> four_gaussian_points_2d_vector_x{-gpt_4, gpt_4, gpt_4, -gpt_4};
    const StdVec<Real> four_gaussian_points_2d_vector_y{-gpt_4, -gpt_4, gpt_4, gpt_4};

    const StdVec<Real> four_gaussian_weights_2d_ = {1, 1, 1, 1};

    int number_of_gaussian_points_;

    StdLargeVec<Matd> &global_b_stress_, &global_b_moment_;
    StdLargeVec<Vecd> &global_b_shear_stress_;
    StdLargeVec<Vecd> &b_n_;
    StdVec<Real> gaussian_point_x;
    StdVec<Real> gaussian_point_y;
    StdVec<Real> gaussian_weight_;
};

/**
 * @class BarStressRelaxationSecondHalf
 * @brief computing stress relaxation process by verlet time stepping
 * This is the second step
 */
class BarStressRelaxationSecondHalf : public BaseBarRelaxation
{
  public:
    explicit BarStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
        : BaseBarRelaxation(inner_relation){};
    virtual ~BarStressRelaxationSecondHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &vel_n_i = vel_[index_i];
        const Vecd &dpseudo_n_dt_i = dpseudo_n_dt_[index_i];
        const Vecd &dpseudo_b_n_dt_i = dpseudo_b_n_dt_[index_i];
        const Matd &transformation_matrix_i = transformation_matrix_[index_i];

        Matd deformation_gradient_change_rate_part_one = Matd::Zero();
        Matd deformation_gradient_change_rate_part_three = Matd::Zero();
        Matd deformation_gradient_change_rate_part_two = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            deformation_gradient_change_rate_part_one -= (vel_n_i - vel_[index_j]) * gradW_ijV_j.transpose();
            deformation_gradient_change_rate_part_two -= (dpseudo_n_dt_i - dpseudo_n_dt_[index_j]) * gradW_ijV_j.transpose();
            deformation_gradient_change_rate_part_three -= (dpseudo_b_n_dt_i - dpseudo_b_n_dt_[index_j]) * gradW_ijV_j.transpose();
        }
        dF_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_one *
                          transformation_matrix_i.transpose() * B_[index_i];
        dF_dt_[index_i].col(Dimensions - 1) = transformation_matrix_i * dpseudo_n_dt_[index_i];
        dF_dt_[index_i].col(Dimensions - 2) = transformation_matrix_i * dpseudo_b_n_dt_[index_i];
        dF_bending_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_two *
                                  transformation_matrix_i.transpose() * B_[index_i];
        dF_b_bending_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_three *
                                    transformation_matrix_i.transpose() * B_[index_i];
    };

    void update(size_t index_i, Real dt = 0.0);
};

/**@class ConstrainBarBodyRegion
 * @brief Fix the position and angle of a bar body part.
 */
class ConstrainBarBodyRegion : public BaseLocalDynamics<BodyPartByParticle>, public BarDataSimple
{
  public:
    ConstrainBarBodyRegion(BodyPartByParticle &body_part);
    virtual ~ConstrainBarBodyRegion(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_, &angular_vel_, &angular_b_vel_;
};

/**@class ConstrainBarBodyRegionAlongAxis
 * @brief The boundary conditions are denoted by SS1 according to the references.
 * The axis must be 0 or 1.
 * Note that the average values for FSI are prescribed also.
 */
class ConstrainBarBodyRegionAlongAxis : public BaseLocalDynamics<BodyPartByParticle>, public BarDataSimple
{
  public:
    ConstrainBarBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis);
    virtual ~ConstrainBarBodyRegionAlongAxis(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    const int axis_; /**< the axis direction for bounding*/
    StdLargeVec<Vecd> &pos_, &pos0_;
    StdLargeVec<Vecd> &vel_, &acc_;
    StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
    StdLargeVec<Vecd> &rotation_b_, &angular_b_vel_, &dangular_b_vel_dt_;
};

/**
 * @class DistributingPointForcesToBar
 * @brief Distribute a series of point forces to its contact Bar bodies.
 */
class DistributingPointForcesToBar : public LocalDynamics, public BarDataSimple
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
    DistributingPointForcesToBar(SPHBody &sph_body, std::vector<Vecd> point_forces,
                                 std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                                 Real particle_spacing_ref, Real h_spacing_ratio = 1.6);
    virtual ~DistributingPointForcesToBar(){};

    virtual void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace slender_structure_dynamics
} // namespace SPH
#endif // THIN_STRUCTURE_DYNAMICS_H
