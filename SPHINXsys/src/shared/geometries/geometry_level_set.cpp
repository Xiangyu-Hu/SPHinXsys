/**
 * @file 	geometry_level_set.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "geometry_level_set.h"

#include "level_set.h"
#include "base_body.h"
#include "in_output.h"
#include "sph_system.h"


namespace SPH {
	//=================================================================================================//
	LevelSetComplexShape::
		LevelSetComplexShape(SPHBody* sph_body, ComplexShape& complex_shape, bool isCleaned)
		: ComplexShape(complex_shape), level_set_(NULL)
	{
		name_ = sph_body->getBodyName();
		level_set_ = sph_body->particle_adaptation_->createLevelSet(complex_shape);
		if (isCleaned) level_set_->cleanInterface();

		In_Output* in_output = sph_body->getSPHSystem().in_output_;
		WriteMeshToPlt 	write_level_set(*in_output, sph_body, level_set_);
		write_level_set.WriteToFile(0.0);
	}
	//=================================================================================================//
	bool LevelSetComplexShape::checkContain(Vecd& input_pnt, bool BOUNDARY_INCLUDED)
	{
		return level_set_->probeSignedDistance(input_pnt) < 0.0 ? true : false;
	}
	//=================================================================================================//
	bool LevelSetComplexShape::checkNearSurface(Vecd& input_pnt, Real threshold)
	{
		if(!checkNotFar(input_pnt, threshold)) return false;
		return  getMaxAbsoluteElement(findSignedDistance(input_pnt) 
			* findNormalDirection(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	Real LevelSetComplexShape::findSignedDistance(Vecd& input_pnt)
	{
		return level_set_->probeSignedDistance(input_pnt);
	}
	//=================================================================================================//
	Vecd LevelSetComplexShape::findNormalDirection(Vecd& input_pnt)
	{
		return level_set_->probeNormalDirection(input_pnt);
	}
	//=================================================================================================//
	bool LevelSetComplexShape::checkNotFar(Vecd& input_pnt, Real threshold)
	{
		return level_set_->isWithinMeshBound(input_pnt);
	}
	//=================================================================================================//
	Real LevelSetComplexShape::computeKernelIntegral(Vecd& input_pnt, Real h_ratio)
	{
		return level_set_->probeKernelIntegral(input_pnt, h_ratio);
	}
	//=================================================================================================//
	Vecd LevelSetComplexShape::computeKernelGradientIntegral(Vecd& input_pnt, Real h_ratio)
	{
		return level_set_->probeKernelGradientIntegral(input_pnt,h_ratio);
	}
}
