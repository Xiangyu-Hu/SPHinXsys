#include "thin_structure_dynamics.h"
#include "thin_structure_math.h"

namespace SPH
{
	//=====================================================================================================//
	namespace thin_structure_dynamics
	{
		//=================================================================================================//
		ShellDynamicsInitialCondition::ShellDynamicsInitialCondition(SPHBody &sph_body)
			: LocalDynamics(sph_body), ShellDataSimple(sph_body),
			  n0_(particles_->n0_), n_(particles_->n_), pseudo_n_(particles_->pseudo_n_),
			  pos0_(particles_->pos0_), transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		ShellAcousticTimeStepSize::ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL)
			: LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)),
			  ShellDataSimple(sph_body), CFL_(CFL), vel_(particles_->vel_), acc_(particles_->acc_),
			  angular_vel_(particles_->angular_vel_), dangular_vel_dt_(particles_->dangular_vel_dt_),
			  thickness_(particles_->thickness_),
			  smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
			  rho0_(particles_->elastic_solid_.ReferenceDensity()),
			  E0_(particles_->elastic_solid_.YoungsModulus()),
			  nu_(particles_->elastic_solid_.PoissonRatio()),
			  c0_(particles_->elastic_solid_.ReferenceSoundSpeed()) {}
		//=================================================================================================//
		Real ShellAcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			// Since the particle does not change its configuration in pressure relaxation step,
			// I chose a time-step size according to Eulerian method.
			Real time_setp_0 = CFL_ * SMIN(sqrt(smoothing_length_ / (acc_[index_i].norm() + TinyReal)),
										   smoothing_length_ / (c0_ + vel_[index_i].norm()));
			Real time_setp_1 = CFL_ * SMIN(sqrt(1.0 / (dangular_vel_dt_[index_i].norm() + TinyReal)),
										   1.0 / (angular_vel_[index_i].norm() + TinyReal));
			Real time_setp_2 = smoothing_length_ * sqrt(rho0_ * (1.0 - nu_ * nu_) / E0_ /
														(2.0 + (Pi * Pi / 12.0) * (1.0 - nu_) *
																   (1.0 + 1.5 * powerN(smoothing_length_ / thickness_[index_i], 2))));
			return SMIN(time_setp_0, time_setp_1, time_setp_2);
		}
		//=================================================================================================//
		ShellCorrectConfiguration::
			ShellCorrectConfiguration(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), ShellDataInner(inner_relation),
			  B_(particles_->B_),
			  n0_(particles_->n0_), transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		void ShellCorrectConfiguration::interaction(size_t index_i, Real dt)
		{
			/** A small number is added to diagonal to avoid dividing by zero. */
			Matd global_configuration = Eps * Matd::Identity();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
				global_configuration += r_ji * gradW_ijV_j.transpose();
			}
			Matd local_configuration =
				transformation_matrix_[index_i] * global_configuration * transformation_matrix_[index_i].transpose();
			/** correction matrix is obtained from local configuration. */
			B_[index_i] = getCorrectionMatrix(local_configuration);
		}
		//=================================================================================================//
		ShellDeformationGradientTensor::
			ShellDeformationGradientTensor(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), ShellDataInner(inner_relation),
			  pos_(particles_->pos_), pseudo_n_(particles_->pseudo_n_), n0_(particles_->n0_),
			  B_(particles_->B_), F_(particles_->F_), F_bending_(particles_->F_bending_),
			  transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		void ShellDeformationGradientTensor::interaction(size_t index_i, Real dt)
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
		}
		//=================================================================================================//
		BaseShellRelaxation::BaseShellRelaxation(BaseInnerRelation &inner_relation)
			: LocalDynamics(inner_relation.sph_body_), ShellDataInner(inner_relation),
			  rho_(particles_->rho_),
			  thickness_(particles_->thickness_),
			  pos_(particles_->pos_), vel_(particles_->vel_),
			  acc_(particles_->acc_),
			  acc_prior_(particles_->acc_prior_),
			  n0_(particles_->n0_), pseudo_n_(particles_->pseudo_n_),
			  dpseudo_n_dt_(particles_->dpseudo_n_dt_), dpseudo_n_d2t_(particles_->dpseudo_n_d2t_),
			  rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
			  dangular_vel_dt_(particles_->dangular_vel_dt_),
			  B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
			  F_bending_(particles_->F_bending_), dF_bending_dt_(particles_->dF_bending_dt_),
			  transformation_matrix_(particles_->transformation_matrix_) {}
		//=================================================================================================//
		ShellStressRelaxationFirstHalf::
			ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
										   int number_of_gaussian_points, bool hourglass_control)
			: BaseShellRelaxation(inner_relation),
			  elastic_solid_(particles_->elastic_solid_),
			  global_stress_(particles_->global_stress_),
			  global_moment_(particles_->global_moment_),
			  global_shear_stress_(particles_->global_shear_stress_),
			  n_(particles_->n_),
			  number_of_gaussian_points_(number_of_gaussian_points),
			  hourglass_control_(hourglass_control),
			  rho0_(elastic_solid_.ReferenceDensity()),
			  inv_rho0_(1.0 / rho0_),
			  smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
			  E0_(elastic_solid_.YoungsModulus()),
			  G0_(elastic_solid_.ShearModulus()),
			  nu_(elastic_solid_.PoissonRatio())
		{
			/** Note that, only three-point and five-point Gaussian quadrature rules are defined. */
			switch (number_of_gaussian_points)
			{
			case 1:
				gaussian_point_ = one_gaussian_point_;
				gaussian_weight_ = one_gaussian_weight_;
				break;
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
		void ShellStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
		{
			// Note that F_[index_i], F_bending_[index_i], dF_dt_[index_i], dF_bending_dt_[index_i]
			// and rotation_[index_i], angular_vel_[index_i], dangular_vel_dt_[index_i], B_[index_i]
			// are defined in local coordinates, while others in global coordinates.
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
			pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;

			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
			rho_[index_i] = rho0_ / F_[index_i].determinant();

			/** Calculate the current normal direction of mid-surface. */
			n_[index_i] = transformation_matrix_[index_i].transpose() * getNormalFromDeformationGradientTensor(F_[index_i]);
			/** Get transformation matrix from global coordinates to current local coordinates. */
			Matd current_transformation_matrix = getTransformationMatrix(n_[index_i]);

			Matd resultant_stress = Matd::Zero();
			Matd resultant_moment = Matd::Zero();
			Vecd resultant_shear_stress = Vecd::Zero();

			for (int i = 0; i != number_of_gaussian_points_; ++i)
			{
				Matd F_gaussian_point = F_[index_i] + gaussian_point_[i] * F_bending_[index_i] * thickness_[index_i] * 0.5;
				Matd dF_gaussian_point_dt = dF_dt_[index_i] + gaussian_point_[i] * dF_bending_dt_[index_i] * thickness_[index_i] * 0.5;
				Matd inverse_F_gaussian_point = F_gaussian_point.inverse();
				Matd current_local_almansi_strain = current_transformation_matrix * transformation_matrix_[index_i].transpose() * 0.5 * 
					(Matd::Identity() - inverse_F_gaussian_point.transpose() * inverse_F_gaussian_point) * 
					transformation_matrix_[index_i] * current_transformation_matrix.transpose();

				/** correct Almansi strain tensor according to plane stress problem. */
				current_local_almansi_strain = getCorrectedAlmansiStrain(current_local_almansi_strain, nu_);
				
				Matd cauchy_stress = elastic_solid_.StressCauchy(current_local_almansi_strain, F_gaussian_point, index_i) + 
					current_transformation_matrix * transformation_matrix_[index_i].transpose() * F_gaussian_point * 
					elastic_solid_.NumericalDampingRightCauchy(F_gaussian_point, dF_gaussian_point_dt, smoothing_length_, index_i) 
					* F_gaussian_point.transpose() * transformation_matrix_[index_i] * current_transformation_matrix.transpose() / F_gaussian_point.determinant();

				/** Impose modeling assumptions. */
				cauchy_stress.col(Dimensions - 1) *= shear_correction_factor_;
				cauchy_stress.row(Dimensions - 1) *= shear_correction_factor_;
				cauchy_stress(Dimensions - 1, Dimensions - 1) = 0.0;

				Matd stress_PK2_gaussian_point = F_gaussian_point.determinant() * F_gaussian_point.inverse() * transformation_matrix_[index_i] * 
					current_transformation_matrix.transpose() * cauchy_stress * current_transformation_matrix * transformation_matrix_[index_i].transpose() 
					* F_gaussian_point.inverse().transpose();

				Vecd shear_stress_PK2_gaussian_point = -stress_PK2_gaussian_point.col(Dimensions - 1);
				Matd moment_PK2_gaussian_point = stress_PK2_gaussian_point * gaussian_point_[i] * thickness_[index_i] * 0.5;

				resultant_stress +=
					0.5 * thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * stress_PK2_gaussian_point;
				resultant_moment +=
					0.5 * thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * moment_PK2_gaussian_point;
				resultant_shear_stress +=
					0.5 * thickness_[index_i] * gaussian_weight_[i] * F_gaussian_point * shear_stress_PK2_gaussian_point;
			}
			/** Only one (for 2D) or two (for 3D) angular momentum equations left. */
			resultant_moment.col(Dimensions - 1) = Vecd::Zero();
			resultant_moment.row(Dimensions - 1) = Vecd::Zero().transpose();
			resultant_shear_stress[Dimensions - 1] = 0.0;

			/** stress and moment in global coordinates for pair interaction */
			global_stress_[index_i] = transformation_matrix_[index_i].transpose() * resultant_stress * transformation_matrix_[index_i];
			global_moment_[index_i] = transformation_matrix_[index_i].transpose() * resultant_moment * transformation_matrix_[index_i];
			global_shear_stress_[index_i] = transformation_matrix_[index_i].transpose() * resultant_shear_stress;
		}
		//=================================================================================================//
		void ShellStressRelaxationFirstHalf::interaction(size_t index_i, Real dt)
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
					Real dim_inv_r_ij = Dimensions / r_ij;
					Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
					Vecd pos_jump = getLinearVariableJump(e_ij, r_ij, pos_[index_i], F_[index_i], pos_[index_j], F_[index_j]);
					acceleration += hourglass_control_factor_ * weight * E0_ * pos_jump * dim_inv_r_ij * 
									inner_neighborhood.dW_ijV_j_[n] * thickness_[index_i];

					Vecd pseudo_n_jump = getLinearVariableJump(e_ij, r_ij, pseudo_n_[index_i] - n0_[index_i],
						transformation_matrix_[index_i].transpose() * F_bending_[index_i] * transformation_matrix_[index_i],
						pseudo_n_[index_j] - n0_[index_j],
						transformation_matrix_[index_j].transpose() * F_bending_[index_j] * transformation_matrix_[index_j]);
					Vecd rotation_jump = getRotationJump(pseudo_n_jump, transformation_matrix_[index_i]);
					pseudo_normal_acceleration += hourglass_control_factor_ / 3.0 * weight * Dimensions * r_ij * G0_ * rotation_jump * 
												  inner_neighborhood.dW_ijV_j_[n] * thickness_[index_i];
				}

				acceleration += (global_stress_i + global_stress_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				pseudo_normal_acceleration += (global_moment_i + global_moment_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			}
	
			acc_[index_i] = acceleration * inv_rho0_ / thickness_[index_i];
			dpseudo_n_d2t_[index_i] = pseudo_normal_acceleration * inv_rho0_ * 12.0 / powerN(thickness_[index_i], 3);

			/** the relation between pseudo-normal and rotations */
			Vecd local_dpseudo_n_d2t = transformation_matrix_[index_i] * dpseudo_n_d2t_[index_i];
			dangular_vel_dt_[index_i] = getRotationFromPseudoNormalForFiniteDeformation(local_dpseudo_n_d2t, rotation_[index_i], angular_vel_[index_i], dt);
		}
		//=================================================================================================//
		void ShellStressRelaxationFirstHalf::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
			angular_vel_[index_i] += dangular_vel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
			rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
			dpseudo_n_dt_[index_i] = transformation_matrix_[index_i].transpose() *
									 getVectorChangeRateAfterThinStructureRotation(local_pseudo_n_0, rotation_[index_i], angular_vel_[index_i]);
			pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::interaction(size_t index_i, Real dt)
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
		}
		//=================================================================================================//
		void ShellStressRelaxationSecondHalf::update(size_t index_i, Real dt)
		{
			F_[index_i] += dF_dt_[index_i] * dt * 0.5;
			F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		ConstrainShellBodyRegion::
			ConstrainShellBodyRegion(BodyPartByParticle &body_part)
			: LocalDynamics(body_part.getSPHBody()), ShellDataSimple(sph_body_),
			  vel_(particles_->vel_), angular_vel_(particles_->angular_vel_)
		{}
		//=================================================================================================//
		void ConstrainShellBodyRegion::update(size_t index_i, Real dt)
		{
			vel_[index_i] = Vecd::Zero();
			angular_vel_[index_i] = Vecd::Zero();
		}
		//=================================================================================================//
		ConstrainShellBodyRegionAlongAxis::ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis)
			: LocalDynamics(body_part.getSPHBody())
			, ShellDataSimple(sph_body_)
			, axis_(axis), pos_(particles_->pos_)
			, pos0_(particles_->pos0_)
			, vel_(particles_->vel_)
			, acc_(particles_->acc_)
			, rotation_(particles_->rotation_)
			, angular_vel_(particles_->angular_vel_)
			, dangular_vel_dt_(particles_->dangular_vel_dt_) 
		{}
		//=================================================================================================//
		void ConstrainShellBodyRegionAlongAxis::update(size_t index_i, Real dt)
		{
			vel_[index_i][axis_] = 0.0;
			vel_[index_i][2] = 0.0;
			acc_[index_i][axis_] = 0.0;
			acc_[index_i][2] = 0.0;

			angular_vel_[index_i][1 - axis_] = 0.0;
			dangular_vel_dt_[index_i][1 - axis_] = 0.0;
		}
		//=================================================================================================//
		DistributingPointForcesToShell::
			DistributingPointForcesToShell(SPHBody &sph_body, std::vector<Vecd> point_forces,
										   std::vector<Vecd> reference_positions, Real time_to_full_external_force,
										   Real particle_spacing_ref, Real h_spacing_ratio)
			: LocalDynamics(sph_body), ShellDataSimple(sph_body),
			  point_forces_(point_forces), reference_positions_(reference_positions),
			  time_to_full_external_force_(time_to_full_external_force),
			  particle_spacing_ref_(particle_spacing_ref), h_spacing_ratio_(h_spacing_ratio),
			  pos0_(particles_->pos0_), acc_prior_(particles_->acc_prior_),
			  thickness_(particles_->thickness_)
		{
			for (int i = 0; i < point_forces_.size(); i++)
			{
				weight_.push_back(StdLargeVec<Real>(0.0));
				time_dependent_point_forces_.push_back(Vecd::Zero());
				sum_of_weight_.push_back(0.0);
				particles_->registerVariable(weight_[i], "Weight_" + std::to_string(i));
			}

			getWeight(); // TODO: should be revised and parallelized, using SimpleDynamics
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::getWeight()
		{
			Kernel *kernel_ = sph_body_.sph_adaptation_->getKernel();
			Real reference_smoothing_length = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
			Real smoothing_length = h_spacing_ratio_ * particle_spacing_ref_;
			Real h_ratio = reference_smoothing_length / smoothing_length;
			Real cutoff_radius_sqr = powerN(2.0 * smoothing_length, 2);
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				sum_of_weight_[i] = 0.0;
				for (size_t index = 0; index < particles_->total_real_particles_; ++index)
				{
					weight_[i][index] = 0.0;
					Vecd displacement = reference_positions_[i] - pos0_[index];
					if (displacement.squaredNorm() <= cutoff_radius_sqr)
					{
						weight_[i][index] = kernel_->W(h_ratio, displacement.norm(), displacement);
						sum_of_weight_[i] += weight_[i][index];
					}
				}
			}
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::setupDynamics(Real dt)
		{
			Real current_time = GlobalStaticVariables::physical_time_;
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				time_dependent_point_forces_[i] = current_time < time_to_full_external_force_
													  ? current_time * point_forces_[i] / time_to_full_external_force_
													  : point_forces_[i];
			}
		}
		//=================================================================================================//
		void DistributingPointForcesToShell::update(size_t index_i, Real dt)
		{
			acc_prior_[index_i] = Vecd::Zero();
			for (int i = 0; i < point_forces_.size(); ++i)
			{
				Vecd force = weight_[i][index_i] / (sum_of_weight_[i] + TinyReal) * time_dependent_point_forces_[i];
				acc_prior_[index_i] += force / particles_->ParticleMass(index_i);
			}
		}
		//=================================================================================================//
	}
}
