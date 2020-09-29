/**
 * @file 	geometry_level_set.cpp
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "geometry_level_set.h"
#include "level_set.h"
#include "base_body.h"
#include "in_output.h"

namespace SPH {
	//=================================================================================================//
	LevelSetComplexShape::
		LevelSetComplexShape(SPHBody* sph_body, ComplexShape& complex_shape, bool isCleaned)
		: ComplexShape(complex_shape), level_set_(NULL)
	{
		name_ = sph_body->GetBodyName();
		Vecd lower_bound, upper_bound;
		findBounds(lower_bound, upper_bound);
		Real mesh_spacing = 4.0 * sph_body->particle_spacing_;
		level_set_ = new LevelSet(complex_shape, lower_bound, upper_bound,	mesh_spacing, 4	);
		if (isCleaned) level_set_->cleanInterface();

		In_Output in_output(sph_body->getSPHSystem());
		WriteLevelSetToPlt 	write_level_set(in_output, { sph_body }, level_set_);
		write_level_set.WriteToFile(0.0);
	}
	//=================================================================================================//
	bool LevelSetComplexShape::checkContain(Vecd input_pnt, bool BOUNDARY_INCLUDED)
	{
		return level_set_->probeLevelSet(input_pnt) < 0.0 ? true : false;
	}
	//=================================================================================================//
	Vecd LevelSetComplexShape::findClosestPoint(Vecd input_pnt)
	{
		return  input_pnt - level_set_->probeLevelSet(input_pnt) * level_set_->probeNormalDirection(input_pnt);
	}
	//=================================================================================================//
	Real LevelSetComplexShape::findSignedDistance(Vecd input_pnt)
	{
		return level_set_->probeLevelSet(input_pnt);
	}
	//=================================================================================================//
	Vecd LevelSetComplexShape::findNormalDirection(Vecd input_pnt)
	{
		return level_set_->probeNormalDirection(input_pnt);
	}
	//=================================================================================================//
	bool LevelSetComplexShape::checkNotFar(Vecd input_pnt, Real threshold)
	{
		return level_set_->isWithinMeshBound(input_pnt);
	}
	//=================================================================================================//
	Vecd LevelSetComplexShape::computeKernelIntegral(Vecd input_pnt, Kernel * kernel)
	{
		return level_set_->computeKernelIntegral(input_pnt, kernel);
	}
}
