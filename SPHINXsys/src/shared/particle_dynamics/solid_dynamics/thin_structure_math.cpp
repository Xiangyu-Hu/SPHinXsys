/**
 * @file 	thin_structure_math.cpp
 * @author	Dong Wu, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "thin_structure_math.h"

using namespace SimTK;

namespace SPH
{
	namespace thin_structure_dynamics
	{
		//=================================================================================================//
		Vec2d getVectorAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles)
		{
			/**The rotation matrix. */
			Mat2d rotation_matrix(0.0);
			rotation_matrix[0][0] = cos(rotation_angles[0]);
			rotation_matrix[0][1] = sin(rotation_angles[0]);
			rotation_matrix[1][0] = -rotation_matrix[0][1];
			rotation_matrix[1][1] = rotation_matrix[0][0];

			return rotation_matrix * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles)
		{
			/**The rotation matrix about the X-axis. */
			Mat3d rotation_matrix_x(0.0);
			rotation_matrix_x[0][0] = 1.0;
			rotation_matrix_x[1][1] = cos(rotation_angles[0]);
			rotation_matrix_x[1][2] = -sin(rotation_angles[0]);
			rotation_matrix_x[2][1] = -rotation_matrix_x[1][2];
			rotation_matrix_x[2][2] = rotation_matrix_x[1][1];
			/**The rotation matrix about the Y-axis. */
			Mat3d rotation_matrix_y(0.0);
			rotation_matrix_y[0][0] = cos(rotation_angles[1]);
			rotation_matrix_y[0][2] = sin(rotation_angles[1]);
			rotation_matrix_y[1][1] = 1.0;
			rotation_matrix_y[2][0] = -rotation_matrix_y[0][2];
			rotation_matrix_y[2][2] = rotation_matrix_y[0][0];

			return rotation_matrix_y * rotation_matrix_x * initial_vector;
		}
		//=================================================================================================//
		Vec2d getVectorChangeRateAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles, const Vec2d &angular_vel)
		{
			/**The derivative of the rotation matrix. */
			Mat2d drotation_matrix_dt(0.0);
			drotation_matrix_dt[0][0] = -sin(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_dt[0][1] = cos(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_dt[1][0] = -drotation_matrix_dt[0][1];
			drotation_matrix_dt[1][1] = drotation_matrix_dt[0][0];

			return drotation_matrix_dt * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorChangeRateAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles, const Vec3d &angular_vel)
		{
			/**The rotation matrix about the X-axis. */
			Mat3d rotation_matrix_x(0.0);
			rotation_matrix_x[0][0] = 1.0;
			rotation_matrix_x[1][1] = cos(rotation_angles[0]);
			rotation_matrix_x[1][2] = -sin(rotation_angles[0]);
			rotation_matrix_x[2][1] = -rotation_matrix_x[1][2];
			rotation_matrix_x[2][2] = rotation_matrix_x[1][1];
			/**The rotation matrix about the Y-axis. */
			Mat3d rotation_matrix_y(0.0);
			rotation_matrix_y[0][0] = cos(rotation_angles[1]);
			rotation_matrix_y[0][2] = sin(rotation_angles[1]);
			rotation_matrix_y[1][1] = 1.0;
			rotation_matrix_y[2][0] = -rotation_matrix_y[0][2];
			rotation_matrix_y[2][2] = rotation_matrix_y[0][0];

			/**The derivative of the rotation matrix of the X-axis. */
			Mat3d drotation_matrix_x_dt(0.0);
			drotation_matrix_x_dt[1][1] = -sin(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_x_dt[1][2] = -cos(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_x_dt[2][1] = -drotation_matrix_x_dt[1][2];
			drotation_matrix_x_dt[2][2] = drotation_matrix_x_dt[1][1];
			/**The derivative of the rotation matrix of the Y-axis. */
			Mat3d drotation_matrix_y_dt(0.0);
			drotation_matrix_y_dt[0][0] = -sin(rotation_angles[1]) * angular_vel[1];
			drotation_matrix_y_dt[0][2] = cos(rotation_angles[1]) * angular_vel[1];
			drotation_matrix_y_dt[2][0] = -drotation_matrix_y_dt[0][2];
			drotation_matrix_y_dt[2][2] = drotation_matrix_y_dt[0][0];

			return (drotation_matrix_y_dt * rotation_matrix_x + rotation_matrix_y * drotation_matrix_x_dt)* initial_vector;
		}
		//=================================================================================================//
		Vec2d getRotationFromPseudoNormalForFiniteDeformation(const Vec2d &dpseudo_n_d2t, const Vec2d &rotation, const Vec2d &angular_vel, Real dt)
		{
			Vec2d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -(dpseudo_n_d2t[0] + sin(rotation[0]) * powerN(angular_vel[0], 2))
								 / (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
		{
			Vec3d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = (dpseudo_n_d2t[1] - sin(rotation[0]) * powerN(angular_vel[0], 2))
								 / (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			dangular_vel_dt[1] = (dpseudo_n_d2t[0] + cos(rotation[0]) * sin(rotation[1])
								  * (powerN(angular_vel[0], 2) + powerN(angular_vel[1], 2))
								  + 2 * sin(rotation[0]) * cos(rotation[1]) * angular_vel[0] * angular_vel[1]
								  + (2 * cos(rotation[0]) * sin(rotation[1]) * angular_vel[0] * dt
								  + 2 * sin(rotation[0]) * cos(rotation[1]) * angular_vel[1] * dt
								  + sin(rotation[0]) * cos(rotation[1])) * dangular_vel_dt[0])
								 / (-2 * sin(rotation[0]) * cos(rotation[1]) * angular_vel[0] * dt
									- 2 * cos(rotation[0]) * sin(rotation[1]) * angular_vel[1] * dt
									+ cos(rotation[0]) * cos(rotation[1]));
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec2d getRotationFromPseudoNormalForSmallDeformation(const Vec2d &dpseudo_n_d2t, const Vec2d &rotation, const Vec2d &angular_vel, Real dt)
		{
			Vec2d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForSmallDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
		{
			Vec3d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -dpseudo_n_d2t[1];
			dangular_vel_dt[1] = dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec2d getNormalFromDeformationGradientTensor(const Mat2d &F)
		{
			Vec2d n = Vec2d(-F.col(0)[1], F.col(0)[0]);
			n = n / (n.norm() + Eps);
			return n;
		}
		//=================================================================================================//
		Vec3d getNormalFromDeformationGradientTensor(const Mat3d &F)
		{
			Vec3d n = F.col(0) % F.col(1);
			n = n / (n.norm() + Eps);
			return n;
		}
		//=================================================================================================//
		Vecd getLinearVariableJump(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			return particle_i_value - particle_j_value
				   + 0.5 * r_ij * (gradient_particle_i_value + gradient_particle_j_value) * e_ij;
		}
		//=================================================================================================//
		Vecd getWENOVariableJump(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			return getWENOLeftState(e_ij, r_ij, particle_i_value,
									gradient_particle_i_value, particle_j_value, gradient_particle_j_value)
				   - getWENORightState(e_ij, r_ij, particle_i_value,
									   gradient_particle_i_value, particle_j_value, gradient_particle_j_value);
		}
		//=================================================================================================//
		Vecd getWENOStateWithStencilPoints(const Vecd &v1, const Vecd &v2, const Vecd &v3, const Vecd &v4)
		{
			Vecd f1 = 0.5 * v2 + 0.5 * v3;
			Vecd f2 = -0.5 * v1 + 1.5 * v2;
			Vecd f3 = v2 / 3.0 + 5.0 * v3 / 6.0 - v4 / 6.0;

			Real epsilon = 1.0e-6;
			Real s1 = dot(v2 - v3, v2 - v3) + epsilon;
			Real s2 = dot(v2 - v1, v2 - v1) + epsilon;
			Real s3 = dot(3.0 * v2 - 4.0 * v3 + v4, 3.0 * v2 - 4.0 * v3 + v4) / 4.0
					  + 13.0 * dot(v2 - 2.0 * v3 + v4, v2 - 2.0 * v3 + v4) / 12.0 + epsilon;
			Real s12 = 13.0 * dot(v1 - 2.0 * v2 + v3, v1 - 2.0 * v2 + v3) / 12.0
					   + dot(v1 - v3, v1 - v3) / 4.0 + epsilon;
			Real s4 = (dot(v1, 6649.0 * v1 - 30414.0 * v2 + 23094.0 * v3 - 5978.0 * v4)
					   + 3.0 * dot(v2, 13667.0 * v2 - 23534.0 * v3 + 6338.0 * v4)
					   + 3.0 * dot(v3, 11147.0 * v3 - 6458.0 * v4)
					   + 3169.0 * dot(v4, v4)) / 2880.0;
			Real tau_4 = s4 - 0.5 * (s1 + s2);

			Real alpha_1 = (1.0 + (tau_4 / s1) * (tau_4 / s12)) / 3.0;
			Real alpha_2 = (1.0 + (tau_4 / s2) * (tau_4 / s12)) / 6.0;
			Real alpha_3 = (1.0 + tau_4 / s3) / 2.0;
			Real w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
			Real w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
			Real w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);

			return w_1 * f1 + w_2 * f2 + w_3 * f3;
		}
		//=================================================================================================//
		Vecd getWENOLeftState(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			Vecd v1 = particle_i_value - gradient_particle_i_value * e_ij * r_ij;
			Vecd v2 = particle_i_value;
			Vecd v3 = particle_j_value;
			Vecd v4 = particle_j_value + gradient_particle_j_value * e_ij * r_ij;

			return getWENOStateWithStencilPoints(v1, v2, v3, v4);
		}
		//=================================================================================================//
		Vecd getWENORightState(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			Vecd v1 = particle_j_value + gradient_particle_j_value * e_ij * r_ij;
			Vecd v2 = particle_j_value;
			Vecd v3 = particle_i_value;
			Vecd v4 = particle_i_value - gradient_particle_i_value * e_ij * r_ij;

			return getWENOStateWithStencilPoints(v1, v2, v3, v4);
		}
		//=================================================================================================//
		Vec2d getRotationJump(const Vec2d &pseudo_n_jump, const Mat2d &transformation_matrix)
		{
			Vec2d local_rotation_jump(0.0);
			Vec2d local_pseuodo_n_jump = transformation_matrix * pseudo_n_jump;
			local_rotation_jump[0] = local_pseuodo_n_jump[0];
			return ~transformation_matrix * local_rotation_jump;
		}
		//=================================================================================================//
		Vec3d getRotationJump(const Vec3d &pseudo_n_jump, const Mat3d &transformation_matrix)
		{
			Vec3d local_rotation_jump(0.0);
			Vec3d local_pseuodo_n_jump = transformation_matrix * pseudo_n_jump;
			local_rotation_jump[0] = local_pseuodo_n_jump[0];
			local_rotation_jump[1] = local_pseuodo_n_jump[1];
			return ~transformation_matrix * local_rotation_jump;
		}
		//=================================================================================================//
		Mat2d getCorrectedAlmansiStrain(const Mat2d &current_local_almansi_strain, const Real &nu_)
		{
			Mat2d corredted_almansi_strain = current_local_almansi_strain;
			corredted_almansi_strain[1][1] = -nu_ * corredted_almansi_strain[0][0] / (1.0 - nu_);
			return corredted_almansi_strain;
		}
		//=================================================================================================//
		Mat3d getCorrectedAlmansiStrain(const Mat3d &current_local_almansi_strain, const Real &nu_)
		{
			Mat3d corredted_almansi_strain = current_local_almansi_strain;
			corredted_almansi_strain[2][2]
				= -nu_ * (corredted_almansi_strain[0][0] + corredted_almansi_strain[1][1]) / (1.0 - nu_);
			return corredted_almansi_strain;
		}
		//=================================================================================================//
	}
}
