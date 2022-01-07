/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	thin_structure_math.h
* @brief 	Here, we define the math operation for thin structure dynamics. 
* @author	Dong Wu and Xiangyu Hu
* @version	0.1
*/

#ifndef THIN_STRUCTURE_MATH_H
#define THIN_STRUCTURE_MATH_H


#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"

namespace SPH
{
	namespace thin_structure_dynamics
	{
		/**
		* @function getVectorAfterThinStructureRotation
		* @brief Each of these basic vector rotations appears counterclockwise
		* @brief when the axis about which they occur points toward the observer,
		* @brief and the coordinate system is right-handed.
		*/
		Vec2d getVectorAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles);
		Vec3d getVectorAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles);

		/** Vector change rate after rotation. */
		Vec2d getVectorChangeRateAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles, const Vec2d &angular_vel);
		Vec3d getVectorChangeRateAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles, const Vec3d &angular_vel);

		/** get the rotation from pseudo-normal for finite deformation. */
		Vec2d getRotationFromPseudoNormalForFiniteDeformation(const Vec2d& dpseudo_n_d2t, const Vec2d& rotation, const Vec2d& angular_vel, Real dt);
		Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d& dpseudo_n_d2t, const Vec3d& rotation, const Vec3d& angular_vel, Real dt);

		/** get the rotation from pseudo-normal for small deformation. */
		Vec2d getRotationFromPseudoNormalForSmallDeformation(const Vec2d& dpseudo_n_d2t, const Vec2d& rotation, const Vec2d& angular_vel, Real dt);
		Vec3d getRotationFromPseudoNormalForSmallDeformation(const Vec3d& dpseudo_n_d2t, const Vec3d& rotation, const Vec3d& angular_vel, Real dt);

		/** get the current normal direction from deformation gradient tensor. */
		Vec2d getNormalFromDeformationGradientTensor(const Mat2d& F);
		Vec3d getNormalFromDeformationGradientTensor(const Mat3d& F);

		/** general Gauss-Legendre quadrature for up to 3-D surface integrals*/   
		Real gaussianQuadratureIntegral(Real a, Real b, int n, const std::function<Real (Real)>& F);

		/**
		* @class LegendrePolynomialSet
		* @brief returns abscissas (xi) and weights (wi) for 
		* @brief Gauss-Legendre quadrature integration for 
		* @brief n=3 or n=5 number of integration nodes *.
		*/
		class LegendrePolynomialSet
		{
		public:
			LegendrePolynomialSet(int number_of_nodes)
				: numberOfNodes(number_of_nodes), xi(number_of_nodes), wi(number_of_nodes) {
				setAbscissasAndWeights();
			}

			const std::vector<Real>& getAbscissas() const { return xi; }
			const std::vector<Real>& getWeights() const { return wi; }

		private:
			void setAbscissasAndWeights() 
			{
				if (numberOfNodes == 3)
				{
					xi[0] = -0.774596669241483377035853079956;
					xi[1] = 0.000000000000000000000000000000;
					xi[2] = 0.774596669241483377035853079956;

					wi[0] = 0.555555555555555555555555555556;
					wi[1] = 0.888888888888888888888888888889;
					wi[2] = 0.555555555555555555555555555556;
				}
				else if (numberOfNodes == 5)
				{
					xi[0] = -0.906179845938663992797626878299;
					xi[1] = -0.538469310105683091036314420700;
					xi[2] = 0.000000000000000000000000000000;
					xi[3] = 0.538469310105683091036314420700;
					xi[4] = 0.906179845938663992797626878299;

					wi[0] = 0.236926885056189087514264040720;
					wi[1] = 0.478628670499366468041291514836;
					wi[2] = 0.568888888888888888888888888889;
					wi[3] = 0.478628670499366468041291514836;
					wi[4] = 0.236926885056189087514264040720;
				}
				else
				{
					std::cout << "Illegal number of nodes for Gaussian Integration: n = " << numberOfNodes << std::endl;
					std::cout << "Legal values are 3 or 5" << std::endl;
					exit(1);
				}
			}

			const int numberOfNodes;
			std::vector<Real> xi;
			std::vector<Real> wi;
		};
	}
}
#endif //THIN_STRUCTURE_MATH_H