/**
 * @file 	level_set_shape.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "level_set_shape.h"

#include "base_body.h"
#include "in_output.h"
#include "sph_system.h"

namespace SPH
{
	//=================================================================================================//
	LevelSetShape::LevelSetShape(SPHBody *sph_body, Shape &shape, Real refinement_ratio)
		: Shape(shape.getName()),
		  level_set_(level_set_keeper_.movePtr(
			  sph_body->sph_adaptation_->createLevelSet(shape, refinement_ratio)))
	{
		bounding_box_ = shape.getBounds();
		is_bounds_found_ = true;
	}
	//=================================================================================================//
	void LevelSetShape::writeLevelSet(SPHBody &sph_body)
	{
		InOutput *in_output = sph_body.getSPHSystem().in_output_;
		MeshRecordingToPlt write_level_set_to_plt(*in_output, sph_body, level_set_);
		write_level_set_to_plt.writeToFile(0);
	}
	//=================================================================================================//
	LevelSetShape *LevelSetShape::cleanLevelSet(Real small_shift_factor)
	{
		level_set_->cleanInterface(small_shift_factor);
		return this;
	}
	//=================================================================================================//
	bool LevelSetShape::checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED)
	{
		return level_set_->probeSignedDistance(input_pnt) < 0.0 ? true : false;
	}
	//=================================================================================================//
	Vecd LevelSetShape::findClosestPoint(const Vecd &input_pnt)
	{
		Real phi = level_set_->probeSignedDistance(input_pnt);
		Vecd normal = level_set_->probeNormalDirection(input_pnt);
		return input_pnt - phi * normal;
	}
	//=================================================================================================//
	BoundingBox LevelSetShape::findBounds()
	{
		if (!is_bounds_found_)
		{
			std::cout << "\n FAILURE: LevelSetShape bounds should be defined at construction!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		return bounding_box_;
	}
	//=================================================================================================//
	Vecd LevelSetShape::findLevelSetGradient(const Vecd &input_pnt)
	{
		return level_set_->probeLevelSetGradient(input_pnt);
	}
	//=================================================================================================//
	Real LevelSetShape::computeKernelIntegral(const Vecd &input_pnt, Real h_ratio)
	{
		return level_set_->probeKernelIntegral(input_pnt, h_ratio);
	}
	//=================================================================================================//
	Vecd LevelSetShape::computeKernelGradientIntegral(const Vecd &input_pnt, Real h_ratio)
	{
		return level_set_->probeKernelGradientIntegral(input_pnt, h_ratio);
	}
	//=================================================================================================//
}
