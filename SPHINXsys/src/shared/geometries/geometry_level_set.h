/**
* @file geometry_level_set.h
* @brief Here, we define geometry based on level set technique.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
*/

#ifndef GEOMETRY_LEVEL_SET_H
#define GEOMETRY_LEVEL_SET_H

#include "geometry.h"

#include <string>

namespace SPH
{

	class SPHBody;
	class BaseLevelSet;
	/**
	 * @class LevelSetComplexShape
	 * @brief the final geomtrical definition of the SPHBody based on a narrow band level set function
	 * generated from the original ComplexShape
	 */
	class LevelSetComplexShape : public ComplexShape
	{
	public:
		LevelSetComplexShape(SPHBody *sph_body, ComplexShape &complex_shape, bool isCleaned = false);
		virtual ~LevelSetComplexShape(){};

		virtual bool checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual bool checkNotFar(const Vecd &input_pnt, Real threshold) override;
		virtual bool checkNearSurface(const Vecd &input_pnt, Real threshold) override;
		virtual Real findSignedDistance(const Vecd &input_pnt) override;
		virtual Vecd findNormalDirection(const Vecd &input_pnt) override;
		virtual Real computeKernelIntegral(const Vecd &input_pnt, Real h_ratio = 1.0);
		virtual Vecd computeKernelGradientIntegral(const Vecd &input_pnt, Real h_ratio = 1.0);

	protected:
		BaseLevelSet *level_set_; /**< narrow bounded levelset mesh. */
	};
}

#endif //GEOMETRY_LEVEL_SET_H