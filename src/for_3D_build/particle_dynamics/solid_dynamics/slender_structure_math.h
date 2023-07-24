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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	slender_structure_dynamics.h
 * @brief 	Here, we define the math operation for slender structure dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Dong Wu, Chi Zhang and Xiangyu Hu
 */

#ifndef SLENDER_STRUCTURE_MATH_H
#define SLENDER_STRUCTURE_MATH_H

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "beam_particles.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
namespace SPH
{
namespace slender_structure_dynamics
{
/**
 * Each of these basic vector rotations appears counterclockwise
 * when the axis about which they occur points toward the observer,
 * and the coordinate system is right-handed.
 */
Vec3d getVectorAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles);

/** Vector change rate after rotation. */
Vec3d getVectorChangeRateAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles, const Vec3d &angular_vel);

/** get the rotation from pseudo-normal for finite deformation. */
Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt);

/** get the rotation from pseudo-normal for small deformation. */
Vec3d getRotationFromPseudoNormalForSmallDeformation(const Vec3d &dpseudo_b_n_d2t, const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt);
Vec3d getRotationFromPseudoNormalForSmallDeformation_b(const Vec3d &dpseudo_b_n_d2t, const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt);

/** get the current normal direction from deformation gradient tensor. */
Vec3d getNormalFromDeformationGradientTensor(const Mat3d &F);
Vec3d getBinormalFromDeformationGradientTensor(const Mat3d &F);

/** get the corrected Eulerian Almansi strain tensor according to plane stress problem. */
Mat3d getCorrectedAlmansiStrain(const Mat3d &current_local_almansi_strain, const Real &nu_);

/** get the correction matrix. */
Mat3d getCorrectionMatrix(const Mat3d &local_deformation_part_one);
Mat3d getCorrectionMatrix_beam(const Mat3d &local_deformation_part_one);
} // namespace slender_structure_dynamics
} // namespace SPH
#endif // THIN_STRUCTURE_MATH_H