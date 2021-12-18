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
	LevelSetShape::
		LevelSetShape(SPHBody *sph_body, Shape &shape, bool isCleaned, bool write_level_set)
		: Shape(sph_body->getBodyName()), level_set_(nullptr)
	{
		name_ = sph_body->getBodyName();
		bounding_box_ = shape.findBounds();
		level_set_ = level_set_keeper_.movePtr(sph_body->sph_adaptation_->createLevelSet(shape));
		if (isCleaned)
			level_set_->cleanInterface();

		if (write_level_set)
		{
			In_Output *in_output = sph_body->getSPHSystem().in_output_;
			MeshRecordingToPlt write_level_set_to_plt(*in_output, *sph_body, level_set_);
			write_level_set_to_plt.writeToFile(0);
		}
	}
	//=================================================================================================//
	bool LevelSetShape::checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED)
	{
		return level_set_->probeSignedDistance(input_pnt) < 0.0 ? true : false;
	}
	//=================================================================================================//
	bool LevelSetShape::checkNearSurface(const Vecd &input_pnt, Real threshold)
	{
		if (!checkNotFar(input_pnt, threshold))
			return false;
		return getMaxAbsoluteElement(findSignedDistance(input_pnt) *
									 findNormalDirection(input_pnt)) < threshold
				   ? true
				   : false;
	}
	//=================================================================================================//
	Real LevelSetShape::findSignedDistance(const Vecd &input_pnt)
	{
		return level_set_->probeSignedDistance(input_pnt);
	}
	//=================================================================================================//
	Vecd LevelSetShape::findNormalDirection(const Vecd &input_pnt)
	{
		//std::cout << "LevelSetComplexShape::findNormalDirection called" << std::endl; //to check if LevelSetComplexShape::findNormalDirection is called
		return level_set_->probeNormalDirection(input_pnt);
	}
	//=================================================================================================//
	bool LevelSetShape::checkNotFar(const Vecd &input_pnt, Real threshold)
	{
		return level_set_->probeIsWithinMeshBound(input_pnt);
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
	Vecd LevelSetShape::findClosestPoint(const Vecd &input_pnt)
	{
		Real phi = level_set_->probeSignedDistance(input_pnt);
		Vecd normal = level_set_->probeNormalDirection(input_pnt);
		return input_pnt - phi * normal;
	}
	//=================================================================================================//
}
