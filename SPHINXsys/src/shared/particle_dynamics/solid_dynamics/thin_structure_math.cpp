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
		Vec2d getVectorAfterThinStructureRotation(Vec2d &initial_vector, Vec2d &rotation_angles)
		{
			/**The rotation matrix. */
			Mat2d rotation_matrix(0.0);
			rotation_matrix[0][0] = cos(rotation_angles[0]);
			rotation_matrix[0][1] = -sin(rotation_angles[0]);
			rotation_matrix[1][0] = -rotation_matrix[0][1];
			rotation_matrix[1][1] = rotation_matrix[0][0];

			return rotation_matrix * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorAfterThinStructureRotation(Vec3d &initial_vector, Vec3d &rotation_angles)
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
		Vec2d getVectorChangeRateAfterThinStructureRotation(Vec2d &initial_vector, Vec2d &rotation_angles, Vec2d &angular_vel)
		{
			/**The derivative of the rotation matrix. */
			Mat2d drotation_matrix_dt(0.0);
			drotation_matrix_dt[0][0] = -sin(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_dt[0][1] = -cos(rotation_angles[0]) * angular_vel[0];
			drotation_matrix_dt[1][0] = -drotation_matrix_dt[0][1];
			drotation_matrix_dt[1][1] = drotation_matrix_dt[0][0];

			return drotation_matrix_dt * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorChangeRateAfterThinStructureRotation(Vec3d& initial_vector, Vec3d& rotation_angles, Vec3d& angular_vel)
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
		Matd getTransformationMatrix(Vec2d direction_of_y)
		{
			Matd transformation_matrix(0.0);
			transformation_matrix[0][0] = direction_of_y[1];
			transformation_matrix[0][1] = -direction_of_y[0];
			transformation_matrix[1][0] = direction_of_y[0];
			transformation_matrix[1][1] = direction_of_y[1];
			return transformation_matrix;
		}
		//=================================================================================================//
		Matd getTransformationMatrix(Vec3d direction_of_Z)
		{
			Matd transformation_matrix(0.0);
			transformation_matrix[0][0] = direction_of_Z[2] + powern(direction_of_Z[1], 2) / (1 + direction_of_Z[2] + Eps);
			transformation_matrix[0][1] = -direction_of_Z[0] * direction_of_Z[1] / (1 + direction_of_Z[2] + Eps);
			transformation_matrix[0][2] = -direction_of_Z[0];
			transformation_matrix[1][0] = transformation_matrix[0][1];
			transformation_matrix[1][1] = direction_of_Z[2] + powern(direction_of_Z[0], 2) / (1 + direction_of_Z[2] + Eps);
			transformation_matrix[1][2] = -direction_of_Z[1];
			transformation_matrix[2][0] = direction_of_Z[0];
			transformation_matrix[2][1] = direction_of_Z[1];
			transformation_matrix[2][2] = direction_of_Z[2];
			return transformation_matrix;
		}
		//=================================================================================================//
		Vecd getRotationFromPseudoNormalForFiniteDeformation(Vec2d dpseudo_n_d2t, Vec2d rotation, Vec2d angular_vel, Real dt)
		{
			Vecd dangular_vel_dt(0.0);
			dangular_vel_dt[0] = (dpseudo_n_d2t[0] - sin(rotation[0]) * powern(angular_vel[0], 2))
				/ (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vecd getRotationFromPseudoNormalForFiniteDeformation(Vec3d dpseudo_n_d2t, Vec3d rotation, Vec3d angular_vel, Real dt)
		{
			Vecd dangular_vel_dt(0.0);
			dangular_vel_dt[0] = (dpseudo_n_d2t[1] - sin(rotation[0]) * powern(angular_vel[0], 2))
				/ (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			dangular_vel_dt[1] = (dpseudo_n_d2t[0] + cos(rotation[0]) * sin(rotation[1])
				* (powern(angular_vel[0], 2) + powern(angular_vel[1], 2))
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
		Vecd getRotationFromPseudoNormalForSmallDeformation(Vec2d dpseudo_n_d2t, Vec2d rotation, Vec2d angular_vel, Real dt)
		{
			Vecd dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vecd getRotationFromPseudoNormalForSmallDeformation(Vec3d dpseudo_n_d2t, Vec3d rotation, Vec3d angular_vel, Real dt)
		{
			Vecd dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -dpseudo_n_d2t[1];
			dangular_vel_dt[1] = dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
	}
}
