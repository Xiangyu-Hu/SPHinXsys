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
