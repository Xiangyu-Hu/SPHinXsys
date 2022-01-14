/**
 * @file 	thin_structure_math.cpp
 * @author	Dong Wu, Chi ZHang, Massoud Rezavand and Xiangyu Hu
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
		Vec3d getVectorChangeRateAfterThinStructureRotation(const Vec3d& initial_vector, const Vec3d& rotation_angles, const Vec3d& angular_vel)
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
		Vec2d getRotationFromPseudoNormalForFiniteDeformation(const Vec2d& dpseudo_n_d2t, const Vec2d& rotation, const Vec2d& angular_vel, Real dt)
		{
			Vec2d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -(dpseudo_n_d2t[0] + sin(rotation[0]) * powerN(angular_vel[0], 2))
				/ (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d& dpseudo_n_d2t, const Vec3d& rotation, const Vec3d& angular_vel, Real dt)
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
		Vec2d getRotationFromPseudoNormalForSmallDeformation(const Vec2d& dpseudo_n_d2t, const Vec2d& rotation, const Vec2d& angular_vel, Real dt)
		{
			Vec2d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForSmallDeformation(const Vec3d& dpseudo_n_d2t, const Vec3d& rotation, const Vec3d& angular_vel, Real dt)
		{
			Vec3d dangular_vel_dt(0.0);
			dangular_vel_dt[0] = -dpseudo_n_d2t[1];
			dangular_vel_dt[1] = dpseudo_n_d2t[0];
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec2d getNormalFromDeformationGradientTensor(const Mat2d& F)
		{
			Vec2d n = Vec2d(-F.col(0)[1], F.col(0)[0]);
			n = n / (n.norm() + Eps);
			return n;
		}
		//=================================================================================================//
		Vec3d getNormalFromDeformationGradientTensor(const Mat3d& F)
		{
			Vec3d n = F.col(0) % F.col(1);
			n = n / (n.norm() + Eps);
			return n;
		}
		//=================================================================================================//
		Real gaussianQuadratureIntegralTriple(const Vecd& ri, const Vecd& rj, const Real& spacing_ref, Kernel* kernel_) 
		{
			const int number_of_nodes = 3; /** or 5 (typically, 3 nodes are enough) */
			const LegendrePolynomialSet legendreSet(number_of_nodes);
			const std::vector<Real>& abscissa = legendreSet.getAbscissas();
			const std::vector<Real>& weight = legendreSet.getWeights();
			/** number of nodes at each directons x, y and z, respectively. Typically same values */
			int m, n, p;
			m = n = p = number_of_nodes;
			
			const Real a1 = rj[0] - spacing_ref/2.; const Real b1 = rj[0] + spacing_ref/2.;
			const Real c1 = rj[1] - spacing_ref/2.; const Real d1 = rj[1] + spacing_ref/2.;
			const Real e1 = rj[2] - spacing_ref/2.; const Real f1 = rj[2] + spacing_ref/2.;

			Real aj = 0.0;
			const Real h1 = 0.5 * (b1 - a1);
			const Real h2 = 0.5 * (a1 + b1);

			for (int i=0; i != m; i++) {    
				Real jx = 0.0;
				const Real x = h1 * abscissa[i] + h2;     
				const Real k1 = 0.5 * (d1 - c1);
				const Real k2 = 0.5 * (d1 + c1);
				for (int j = 0; j !=n; j++) { 
					Real jy = 0.0;
					const Real y = k1 * abscissa[j] + k2;  
					const Real l1 = 0.5 * (f1 - e1);
					const Real l2 = 0.5 * (f1 + e1);
					for (int k = 0; k != p; k++) { 
						const Real z = l1 * abscissa[k] + l2;
						const Real q = kernel_->W(Vecd(x, y, z).norm(), Vecd(x, y, z));
						jy += weight[k] * q;
					}
					jx += weight[j] * l1 * jy;
				}
				aj += weight[i] * k1 * jx; 
			}
			return aj * h1;
		}
		//=================================================================================================//
		Real gaussianQuadratureIntegralDouble(const Vecd& ri, const Vecd& rj, const Real& spacing_ref, Kernel* kernel_) 
		{
			const int number_of_nodes = 3; /** or 5 (typically, 3 nodes are enough) */
			const LegendrePolynomialSet legendreSet(number_of_nodes);
			const std::vector<Real>& abscissa = legendreSet.getAbscissas();
			const std::vector<Real>& weight = legendreSet.getWeights();
			/** number of nodes at each directons x and y, respectively. Typically same values */
			int m, n;
			m = n = number_of_nodes;

			const Real a1 = rj[0] - spacing_ref/2.; const Real b1 = rj[0] + spacing_ref/2.;
			const Real c1 = rj[1] - spacing_ref/2.; const Real d1 = rj[1] + spacing_ref/2.;

			Real aj = 0.0;
			const Real h1 = 0.5 * (b1 - a1);
			const Real h2 = 0.5 * (a1 + b1);

			for (int i=0; i != m; i++) {    
				Real jx = 0.0;
				const Real x = h1 * abscissa[i] + h2;     
				const Real k1 = 0.5 * (d1 - c1);
				const Real k2 = 0.5 * (d1 + c1);
				for (int j = 0; j !=n; j++) { 
					const Real y = k1 * abscissa[j] + k2;  
					const Real q = kernel_->W(Vecd(x, y).norm(), Vecd(x, y));
					jx += weight[j] * q;
				}
				aj += weight[i] * k1 * jx; 
			}
			printf("aj = %.12f \n", aj);
			return aj * h1;
		}
	}
}
