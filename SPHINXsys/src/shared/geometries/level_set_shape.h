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
* @file level_set_shape.h
* @brief Here, we define geometry based on level set technique.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
*/

#ifndef LEVEL_SET_SHAPE_H
#define LEVEL_SET_SHAPE_H

#include "base_geometry.h"
#include "level_set.h"

#include <string>

namespace SPH
{

	class SPHBody;

	/**
	 * @class LevelSetShape
	 * @brief A shape using level set to define geometry
	 */
	class LevelSetShape : public Shape
	{
	private:
		UniquePtrKeeper<BaseLevelSet> level_set_keeper_;

	public:
		LevelSetShape(SPHBody *sph_body, Shape &shape, bool isCleaned = false, bool write_level_set = true);
		virtual ~LevelSetShape(){};

		virtual BoundingBox findBounds() override { return bounding_box_; };
		virtual Vecd findClosestPoint(const Vecd &input_pnt) override;
		virtual bool checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual bool checkNotFar(const Vecd &input_pnt, Real threshold) override;
		virtual bool checkNearSurface(const Vecd &input_pnt, Real threshold) override;
		virtual Real findSignedDistance(const Vecd &input_pnt) override;
		virtual Vecd findNormalDirection(const Vecd &input_pnt) override;
		virtual Vecd findNoneNormalizedNormalDirection(const Vecd& input_pnt);

		virtual Real computeKernelIntegral(const Vecd &input_pnt, Real h_ratio = 1.0);
		virtual Vecd computeKernelGradientIntegral(const Vecd &input_pnt, Real h_ratio = 1.0);

	protected:
		BoundingBox bounding_box_;
		BaseLevelSet *level_set_; /**< narrow bounded levelset mesh. */
	};
}
#endif //LEVEL_SET_SHAPE_H
