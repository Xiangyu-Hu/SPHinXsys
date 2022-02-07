/**
 * @file 	thin_structure_dynamics.cpp
 * @author	Dong Wu and Xiangyu Hu
 */

#include "thin_structure_dynamics.h"
#include "thin_structure_math.h"

using namespace SimTK;

namespace SPH
{
	namespace thin_structure_dynamics
	{
		//=================================================================================================//
		ShellDynamicsInitialCondition::
			ShellDynamicsInitialCondition(SolidBody &solid_body)
			: ParticleDynamicsSimple(solid_body),
			  ShellDataSimple(solid_body),
			  n_0_(particles_->n_0_), n_(particles_->n_), pseudo_n_(particles_->pseudo_n_),
			  pos_0_(particles_->pos_0_),
			  transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		ShellAcousticTimeStepSize::ShellAcousticTimeStepSize(SolidBody &solid_body)
			: ParticleDynamicsReduce<Real, ReduceMin>(solid_body),
			  ShellDataSimple(solid_body),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  angular_vel_(particles_->angular_vel_), dangular_vel_dt_(particles_->dangular_vel_dt_),
			  shell_thickness_(particles_->shell_thickness_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()),
			  rho0_(material_->ReferenceDensity()),
			  E0_(material_->YoungsModulus()),
			  nu_(material_->PoissonRatio())

		{
			initial_reference_ = DBL_MAX;
		}
		//=================================================================================================//
		Real ShellAcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			// Since the particle does not change its configuration in pressure relaxation step,
			// I chose a time-step size according to Eulerian method.
			Real sound_speed = material_->ReferenceSoundSpeed();
			Real time_setp_0 = 0.6 * SMIN(sqrt(smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal)),
										  smoothing_length_ / (sound_speed + vel_n_[index_i].norm()));
			Real time_setp_1 = 0.6 * SMIN(sqrt(1.0 / (dangular_vel_dt_[index_i].norm() + TinyReal)),
										  1.0 / (angular_vel_[index_i].norm() + TinyReal));
			Real time_setp_2 = smoothing_length_ * sqrt(rho0_ * (1.0 - nu_ * nu_) / E0_ /
														(2.0 + (Pi * Pi / 12.0) * (1.0 - nu_) *
																   (1.0 + 1.5 * powerN(smoothing_length_ / shell_thickness_[index_i], 2))));
			return SMIN(time_setp_0, time_setp_1, time_setp_2);
		}
		//=================================================================================================//
		ShellCorrectConfiguration::
			ShellCorrectConfiguration(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  ShellDataInner(inner_relation),
			  Vol_(particles_->Vol_), B_(particles_->B_),
			  n_0_(particles_->n_0_), transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		void ShellCorrectConfiguration::Interaction(size_t index_i, Real dt)
		{
			/** A small number is added to diagonal to avoid dividing by zero. */
			Matd global_configuration(Eps);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				global_configuration += Vol_[index_j] * SimTK::outer(r_ji, gradw_ij);
			}
			Matd local_configuration =
				transformation_matrix_[index_i] * global_configuration * (~transformation_matrix_[index_i]);
			/** correction matrix is obtained from local configuration. */
			B_[index_i] = SimTK::inverse(local_configuration) * reduced_unit_matrix;
		}
		//=================================================================================================//
		ShellDeformationGradientTensor::
			ShellDeformationGradientTensor(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  ShellDataInner(inner_relation),
			  Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			  pseudo_n_(particles_->pseudo_n_), n_0_(particles_->n_0_),
			  B_(particles_->B_), F_(particles_->F_), F_bending_(particles_->F_bending_),
			  transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		void ShellDeformationGradientTensor::Interaction(size_t index_i, Real dt)
		{
			const Vecd &pseudo_n_i = pseudo_n_[index_i];
			const Vecd &pos_n_i = pos_n_[index_i];
			const Matd &transformation_matrix_i = transformation_matrix_[index_i];

			Matd deformation_part_one(0.0);
			Matd deformation_part_two(0.0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation_part_one -= Vol_[index_j] * SimTK::outer((pos_n_i - pos_n_[index_j]), gradw_ij);
				deformation_part_two -= Vol_[index_j] * SimTK::outer(
															((pseudo_n_i - n_0_[index_i]) - (pseudo_n_[index_j] - n_0_[index_j])), gradw_ij);
			}
			F_[index_i] = transformation_matrix_i * deformation_part_one * (~transformation_matrix_i) * B_[index_i];
			F_[index_i].col(Dimensions - 1) = transformation_matrix_i * pseudo_n_[index_i];
			F_bending_[index_i] = transformation_matrix_i * deformation_part_two * (~transformation_matrix_i) * B_[index_i];
		}
		//=================================================================================================//
		BaseShellRelaxation::BaseShellRelaxation(BaseBodyRelationInner &inner_relation)
			: ParticleDynamics1Level(*inner_relation.sph_body_),
			  ShellDataInner(inner_relation), Vol_(particles_->Vol_),
			  rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			  shell_thickness_(particles_->shell_thickness_),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			  dvel_dt_(particles_->dvel_dt_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  force_from_fluid_(particles_->force_from_fluid_),
			  n_0_(particles_->n_0_), pseudo_n_(particles_->pseudo_n_),
			  dpseudo_n_dt_(particles_->dpseudo_n_dt_), dpseudo_n_d2t_(particles_->dpseudo_n_d2t_),
			  rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
			  dangular_vel_dt_(particles_->dangular_vel_dt_),
			  B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
			  F_bending_(particles_->F_bending_), dF_bending_dt_(particles_->dF_bending_dt_),
			  transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		ShellStressRelaxationFirstHalf::
			ShellStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation,
				int number_of_gaussian_points, bool hourglass_control)
			: BaseShellRelaxation(inner_relation),
			stress_PK1_(particles_->stress_PK1_),
			global_stress_(particles_->global_stress_),
			global_moment_(particles_->global_moment_),
			global_shear_stress_(particles_->global_shear_stress_),
			n_(particles_->n_),
			number_of_gaussian_points_(number_of_gaussian_points),
			hourglass_control_(hourglass_control),
			rho0_(material_->ReferenceDensity()),
			inv_rho0_(1.0 / rho0_),
			smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()),
			E0_(material_->YoungsModulus()),
			G0_(material_->ShearModulus()),
			nu_(material_->PoissonRatio())
		{
			/** Note that, only three-point and five-point Gaussian quadrature rules are defined. */
			switch (number_of_gaussian_points)
			{
			case 5:
				gaussian_point_ = five_gaussian_points_;
				gaussian_weight_ = five_gaussian_weights_;
				break;
			default:
				gaussian_point_ = three_gaussian_points_;
				gaussian_weight_ = three_gaussian_weights_;
			}
			/** Define the factor of hourglass control algorithm according to the dimension. */
			if (Dimensions == 2)
			{
				hourglass_control_factor_ = 0.05;
			}
			else
			{
				hourglass_control_factor_ = 1.0e-4;
			}
		}
		//=================================================================================================//
		void ShellStressRelaxationFirstHalf::Initialization(size_t index_i, Real dt)
		{
			// Note that F_[index_i], F_bending_[index_i], dF_dt_[index_i], dF_bending_dt_[index_i]
			// and rotation_[index_i], angular_vel_[index_i], dangular_vel_dt_[index_i], B_[index_i]
			// are defined in local coordinates, while others in global coordinates.
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
			pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;

			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
			rho_n_[index_i] = rho0_ / det(F_[index_i]);

			/** Calculate the current normal direction of mid-surface. */
			n_[index_i] = (~transformation_matrix_[index_i]) * getNormalFromDeformationGradientTensor(F_[index_i]);
			/** Get transformation matrix from global coordinates to current local coordinates. */
			Matd current_transformation_matrix = getTransformationMatrix(n_[index_i]);

			/** Initialize the local stress to 0. */
			Matd resultant_stress(0);
			Matd resultant_moment(0);
			Vecd resultant_shear_stress(0);
			for (int i = 0; i != number_of_gaussian_points_; ++i)
			{
				Matd F_gaussian_point = F_[index_i] + gaussian_point_[i] * F_bending_[index_i] * shell_thickness_[index_i] * 0.5;
				Matd dF_gaussian_point_dt = dF_dt_[index_i] + gaussian_point_[i] * dF_bending_dt_[index_i] * shell_thickness_[index_i] * 0.5;
				Matd inverse_F_gaussion_point = SimTK::inverse(F_gaussian_point);
				Matd current_local_almansi_strain = current_transformation_matrix * (~transformation_matrix_[index_i])
											* 0.5 * (Matd(1.0) - ~inverse_F_gaussion_point * inverse_F_gaussion_point)
											* transformation_matrix_[index_i] * (~current_transformation_matrix);
				/** correct Almansi strain tensor according to plane stress problem. */
				current_local_almansi_strain = getCorrectedAlmansiStrain(current_local_almansi_strain, nu_);
				Matd cauchy_stress = material_->EulerianConstitutiveRelation(current_local_almansi_strain, F_gaussian_point, index_i)
					+ current_transformation_matrix * (~transformation_matrix_[index_i]) * F_gaussian_point
					* material_->NumericalDampingRightCauchy(F_gaussian_point, dF_gaussian_point_dt, smoothing_length_, index_i)
					* (~F_gaussian_point) * transformation_matrix_[index_i] * (~current_transformation_matrix) / det(F_gaussian_point);

				/** Impose modeling assumptions. */
				cauchy_stress.col(Dimensions - 1) *= shear_correction_factor_;
				cauchy_stress.row(Dimensions - 1) *= shear_correction_factor_;
				cauchy_stress[Dimensions - 1][Dimensions - 1] = 0.0;

				Matd stress_PK2_gaussian_point = det(F_gaussian_point) * SimTK::inverse(F_gaussian_point)
					* transformation_matrix_[index_i] * (~current_transformation_matrix) * cauchy_stress
					* current_transformation_matrix * (~transformation_matrix_[index_i])
					* (~SimTK::inverse(F_gaussian_point));

				/** Get the mid-surface stress to output the von-Mises equivalent stress. */
				if (i == 0) stress_PK1_[index_i] = F_gaussian_point * stress_PK2_gaussian_point;

				Vecd shear_stress_PK2_gaussian_point = -stress_PK2_gaussian_point.col(Dimensions - 1);
				Matd moment_PK2_gaussian_point = stress_PK2_gaussian_point * gaussian_point_[i] * shell_thickness_[index_i] * 0.5;

				resultant_stress +=
					0.5 * shell_thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * stress_PK2_gaussian_point;
				resultant_moment +=
					0.5 * shell_thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * moment_PK2_gaussian_point;
				resultant_shear_stress +=
					0.5 * shell_thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * shear_stress_PK2_gaussian_point;
			}
			/** Only one (for 2D) or two (for 3D) angular momentum equations left. */
			resultant_moment.col(Dimensions - 1) = Vecd(0);
			resultant_moment.row(Dimensions - 1) = ~Vecd(0);
			resultant_shear_stress[Dimensions - 1] = 0.0;

			/** stress and moment in global coordinates for pair interaction */
			global_stress_[index_i] =
				(~transformation_matrix_[index_i]) * resultant_stress * transformation_matrix_[index_i];
			global_moment_[index_i] =
				(~transformation_matrix_[index_i]) * resultant_moment * transformation_matrix_[index_i];
			global_shear_stress_[index_i] = (~transformation_matrix_[index_i]) * resultant_shear_stress;
		}
		//=================================================================================================//
		void ShellStressRelaxationFirstHalf::Interaction(size_t index_i, Real dt)
		{
			const Vecd &global_shear_stress_i = global_shear_stress_[index_i];
			const Matd &global_stress_i = global_stress_[index_i];
			const Matd &global_moment_i = global_moment_[index_i];

			Vecd acceleration(0.0);
			Vecd pseudo_normal_acceleration = global_shear_stress_i;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				/** Whether or not to use hourglass control algorithm. */
				if (hourglass_control_)
				{
					Vecd e_ij = inner_neighborhood.e_ij_[n];
					Real r_ij = inner_neighborhood.r_ij_[n];
					Real dim_inv_r_ij = Dimensions / r_ij;
					Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
					Vecd pos_jump = getLinearVariableJump(e_ij, r_ij, pos_n_[index_i], F_[index_i], pos_n_[index_j], F_[index_j]);
					acceleration += hourglass_control_factor_ * weight * E0_ * pos_jump * dim_inv_r_ij
									* inner_neighborhood.dW_ij_[n] * Vol_[index_j] * shell_thickness_[index_i];

					Vecd pseudo_n_jump = getLinearVariableJump(e_ij, r_ij, pseudo_n_[index_i] - n_0_[index_i],
											F_bending_[index_i], pseudo_n_[index_j] - n_0_[index_j], F_bending_[index_j]);
					Vecd rotation_jump = getRotationJump(pseudo_n_jump, transformation_matrix_[index_i]);
					pseudo_normal_acceleration += hourglass_control_factor_ / 3.0 * weight * Dimensions * r_ij * G0_
								* rotation_jump * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * shell_thickness_[index_i];
				}

				acceleration += (global_stress_i + global_stress_[index_j])
							    * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
				pseudo_normal_acceleration += (global_moment_i + global_moment_[index_j])
											  * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
			}
			/** including external force (body force) and force from fluid */
			dvel_dt_[index_i] = acceleration * inv_rho0_ / shell_thickness_[index_i]
				+ dvel_dt_prior_[index_i] + force_from_fluid_[index_i] / mass_[index_i] / shell_thickness_[index_i];
			dpseudo_n_d2t_[index_i] = pseudo_normal_acceleration * inv_rho0_
				* 12.0 / powerN(shell_thickness_[index_i], 3);

			/** the relation between pseudo-normal and rotations */
			Vecd local_dpseudo_n_d2t = transformation_matrix_[index_i] * dpseudo_n_d2t_[index_i];
			dangular_vel_dt_[index_i] = getRotationFromPseudoNormalForFiniteDeformation
										(local_dpseudo_n_d2t, rotation_[index_i], angular_vel_[index_i], dt);
		}
		//=================================================================================================//
		void ShellStressRelaxationFirstHalf::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
			angular_vel_[index_i] += dangular_vel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
			rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
			dpseudo_n_dt_[index_i] = (~transformation_matrix_[index_i]) *
									 getVectorChangeRateAfterThinStructureRotation(local_pseudo_n_0, rotation_[index_i], angular_vel_[index_i]);
			pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::Interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_n_i = vel_n_[index_i];
			const Vecd &dpseudo_n_dt_i = dpseudo_n_dt_[index_i];
			const Matd &transformation_matrix_i = transformation_matrix_[index_i];

			Matd deformation_gradient_change_rate_part_one(0.0);
			Matd deformation_gradient_change_rate_part_two(0.0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradw_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				deformation_gradient_change_rate_part_one -= Vol_[index_j] * SimTK::outer(
																				 (vel_n_i - vel_n_[index_j]), gradw_ij);
				deformation_gradient_change_rate_part_two -= Vol_[index_j] * SimTK::outer(
																				 (dpseudo_n_dt_i - dpseudo_n_dt_[index_j]), gradw_ij);
			}
			dF_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_one * (~transformation_matrix_i) * B_[index_i];
			dF_dt_[index_i].col(Dimensions - 1) = transformation_matrix_i * dpseudo_n_dt_[index_i];
			dF_bending_dt_[index_i] = transformation_matrix_i * deformation_gradient_change_rate_part_two * (~transformation_matrix_i) * B_[index_i];
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::Update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		ConstrainShellBodyRegion::
			ConstrainShellBodyRegion(SolidBody &solid_body, BodyPartByParticle &body_part)
			: PartSimpleDynamicsByParticle(solid_body, body_part), ShellDataSimple(solid_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_), n_(particles_->n_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_),
			  rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
			  dangular_vel_dt_(particles_->dangular_vel_dt_),
			  pseudo_n_(particles_->pseudo_n_), dpseudo_n_dt_(particles_->dpseudo_n_dt_)
		{
		}
		//=================================================================================================//
		void ConstrainShellBodyRegion::Update(size_t index_i, Real dt)
		{
			Vecd pos_0 = pos_0_[index_i];
			Vecd pos_n = pos_n_[index_i];
			Vecd vel_n = vel_n_[index_i];
			Vecd dvel_dt = dvel_dt_[index_i];
			Vecd rotation_0(0.0);
			Vecd angular_vel(0.0);
			Vecd dangular_vel_dt(0.0);
			Vecd dpseudo_normal_dt(0.0);

			pos_n_[index_i] = getDisplacement(pos_0, pos_n);
			vel_n_[index_i] = getVelocity(pos_0, pos_n, vel_n);
			dvel_dt_[index_i] = GetAcceleration(pos_0, pos_n, dvel_dt);
			rotation_[index_i] = GetRotationAngle(pos_0, pos_n, rotation_0);
			angular_vel_[index_i] = GetAngularVelocity(pos_0, pos_n, angular_vel);
			dangular_vel_dt_[index_i] = GetAngularAcceleration(pos_0, pos_n, dangular_vel_dt);
			pseudo_n_[index_i] = GetPseudoNormal(pos_0, pos_n, local_pseudo_n_0);
			dpseudo_n_dt_[index_i] = GetPseudoNormalChangeRate(pos_0, pos_n, dpseudo_normal_dt);

			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		FixedFreeRotateShellBoundary::
			FixedFreeRotateShellBoundary(BaseBodyRelationInner &inner_relation,
										 BodyPartByParticle &body_part, Vecd constrained_direction)
			: PartInteractionDynamicsByParticle1Level(*inner_relation.sph_body_, body_part),
			  ShellDataInner(inner_relation),
			  W0_(sph_adaptation_->getKernel()->W0(Vecd(0))),
			  constrain_matrix_(Matd(0)), recover_matrix_(Matd(1.0)),
			  Vol_(particles_->Vol_), vel_n_(particles_->vel_n_),
			  angular_vel_(particles_->angular_vel_)
		{
			particles_->registerAVariable<Vecd>(vel_n_temp_, "TemporaryVelocity");
			particles_->registerAVariable<Vecd>(angular_vel_temp_, "TemporaryAngularVelocity");
			for (int k = 0; k != Dimensions; ++k)
			{
				constrain_matrix_[k][k] = constrained_direction[k];
				recover_matrix_[k][k] = 1.0 - constrain_matrix_[k][k];
			}
		}
		//=================================================================================================//
		void FixedFreeRotateShellBoundary::Initialization(size_t index_i, Real dt)
		{
			vel_n_[index_i] = constrain_matrix_ * vel_n_[index_i];
			angular_vel_[index_i] = Vecd(0);
		}
		//=================================================================================================//
		void FixedFreeRotateShellBoundary::Interaction(size_t index_i, Real dt)
		{
			Real ttl_weight = W0_ * Vol_[index_i];
			Vecd vel_i = recover_matrix_ * vel_n_[index_i] * ttl_weight;
			Real ttl_weight_angular = W0_ * Vol_[index_i];
			Vecd angular_vel_i = angular_vel_[index_i] * ttl_weight_angular;

			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real weight_j = inner_neighborhood.W_ij_[n] * Vol_[index_j];

				ttl_weight += weight_j;
				vel_i += recover_matrix_ * vel_n_[index_j] * weight_j;
				//exclude boundary particles to achieve extrapolation
				if (angular_vel_[index_j].norm() >= Eps)
				{
					angular_vel_i += angular_vel_[index_j] * weight_j;
					ttl_weight_angular += weight_j;
				}
			}

			vel_n_temp_[index_i] = vel_i / ttl_weight;
			angular_vel_temp_[index_i] = angular_vel_i / ttl_weight_angular;
		}
		//=================================================================================================//
		void FixedFreeRotateShellBoundary::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += vel_n_temp_[index_i];
			angular_vel_[index_i] = angular_vel_temp_[index_i];
		}
		//=================================================================================================//
		ClampConstrainShellBodyRegion::
			ClampConstrainShellBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part)
			: PartInteractionDynamicsByParticle1Level(*inner_relation.sph_body_, body_part),
			  ShellDataInner(inner_relation),
			  Vol_(particles_->Vol_), vel_n_(particles_->vel_n_),
			  angular_vel_(particles_->angular_vel_)
		{
			particles_->registerAVariable<Vecd>(vel_n_temp_, "TemporaryVelocity");
			particles_->registerAVariable<Vecd>(angular_vel_temp_, "TemporaryAngularVelocity");
		}
		//=================================================================================================//
		void ClampConstrainShellBodyRegion::Initialization(size_t index_i, Real dt)
		{
			vel_n_[index_i] = Vecd(0);
			angular_vel_[index_i] = Vecd(0);
		}
		//=================================================================================================//
		void ClampConstrainShellBodyRegion::Interaction(size_t index_i, Real dt)
		{
			Real ttl_weight(Eps);
			Vecd vel_i = vel_n_[index_i];
			Vecd angular_vel_i = angular_vel_[index_i];

			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real weight_j = inner_neighborhood.W_ij_[n] * Vol_[index_j];

				ttl_weight += weight_j;
				vel_i += vel_n_[index_j] * weight_j;
				angular_vel_i += angular_vel_[index_j] * weight_j;
			}

			vel_n_temp_[index_i] = vel_i / ttl_weight;
			angular_vel_temp_[index_i] = angular_vel_i / ttl_weight;
		}
		//=================================================================================================//
		void ClampConstrainShellBodyRegion::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] = vel_n_temp_[index_i];
			angular_vel_[index_i] = angular_vel_temp_[index_i];
		}
		//=================================================================================================//
		ConstrainShellBodyRegionInAxisDirection::
			ConstrainShellBodyRegionInAxisDirection(SolidBody &solid_body, BodyPartByParticle &body_part, int axis_direction)
			: PartSimpleDynamicsByParticle(solid_body, body_part), ShellDataSimple(solid_body),
			  axis_(axis_direction), pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_),
			  rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
			  dangular_vel_dt_(particles_->dangular_vel_dt_) {}
		//=================================================================================================//
		void ConstrainShellBodyRegionInAxisDirection::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i][axis_] = 0.0;
			vel_n_[index_i][2] = 0.0;
			dvel_dt_[index_i][axis_] = 0.0;
			dvel_dt_[index_i][2] = 0.0;

			angular_vel_[index_i][1 - axis_] = 0.0;
			dangular_vel_dt_[index_i][1 - axis_] = 0.0;

			/** the average values are prescirbed also. */
			vel_ave_[index_i] = vel_n_[index_i];
			dvel_dt_ave_[index_i] = dvel_dt_[index_i];
		}
		//=================================================================================================//
		DistributingPointForcesToShell::
			DistributingPointForcesToShell(SolidBody &solid_body, std::vector<Vecd> point_forces,
				std::vector<Vecd> reference_positions, Real time_to_full_external_force,
				Real particle_spacing_ref, Real h_spacing_ratio)
			: ParticleDynamicsSimple(solid_body), ShellDataSimple(solid_body),
			point_forces_(point_forces), reference_positions_(reference_positions),
			time_to_full_external_force_(time_to_full_external_force),
			particle_spacing_ref_(particle_spacing_ref), h_spacing_ratio_(h_spacing_ratio),
			pos_0_(particles_->pos_0_), dvel_dt_prior_(particles_->dvel_dt_prior_),
			Vol_(particles_->Vol_), mass_(particles_->mass_), shell_thickness_(particles_->shell_thickness_)
		{
			for (int i = 0; i < point_forces_.size(); i++)
			{
				weight_.push_back(StdLargeVec<Real>(0.0));
				time_dependent_point_forces_.push_back(Vecd(0.0));
				sum_of_weight_.push_back(0.0);
				particles_->registerAVariable<Real>(weight_[i], "Weight_" + std::to_string(i));
			}
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::getWeight()
		{
			Kernel *kernel_ = sph_body_->sph_adaptation_->getKernel();
			Real reference_smoothing_length = sph_adaptation_->ReferenceSmoothingLength();
			Real smoothing_length = h_spacing_ratio_ * particle_spacing_ref_;
			Real h_ratio = reference_smoothing_length / smoothing_length;
			Real cutoff_radius_sqr = powerN(2.0 * smoothing_length, 2);
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				sum_of_weight_[i] = 0.0;
				for (size_t index = 0; index < particles_->total_real_particles_; ++index)
				{
					weight_[i][index] = 0.0;
					Vecd displacement = reference_positions_[i] - pos_0_[index];
					if (displacement.normSqr() <= cutoff_radius_sqr)
					{
						weight_[i][index] = kernel_->W(h_ratio, displacement.norm(), displacement);
						sum_of_weight_[i] += weight_[i][index];
					}
				}

			}
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::getForce()
		{
			Real current_time = GlobalStaticVariables::physical_time_;
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				time_dependent_point_forces_[i] = current_time < time_to_full_external_force_ ?
					current_time * point_forces_[i] / time_to_full_external_force_ : point_forces_[i];
			}
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::Update(size_t index_i, Real dt)
		{
			dvel_dt_prior_[index_i] = 0.0;
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				Vecd force = weight_[i][index_i] / (sum_of_weight_[i] + TinyReal) * time_dependent_point_forces_[i];
				dvel_dt_prior_[index_i] += force / mass_[index_i] / shell_thickness_[index_i];
			}
		}
		//=================================================================================================//
	}
}
