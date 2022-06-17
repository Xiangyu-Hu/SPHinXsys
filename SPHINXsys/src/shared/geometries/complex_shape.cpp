/**
 * @file 	complex_geometries.cpp
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#include "complex_shape.h"
namespace SPH
{
	//=================================================================================================//
	bool AlignedBoxShape::checkInBounds(int axis, const Vecd &point)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		return point_origin[axis] >= -halfsize_[axis] && point_origin[axis] <= halfsize_[axis]
				   ? true
				   : false;
	}
	//=================================================================================================//
	bool AlignedBoxShape::checkUpperBound(int axis, const Vecd &point)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		return point_origin[axis] > halfsize_[axis] ? true : false;
	}
	//=================================================================================================//
	bool AlignedBoxShape::checkLowerBound(int axis, const Vecd &point)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		return point_origin[axis] < -halfsize_[axis] ? true : false;
	}
	//=================================================================================================//
	bool AlignedBoxShape::checkNearUpperBound(int axis, const Vecd &point, Real threshold)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		return ABS(point_origin[axis] - halfsize_[axis]) <= threshold ? true : false;
	}
	//=================================================================================================//
	bool AlignedBoxShape::checkNearLowerBound(int axis, const Vecd &point, Real threshold)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		return ABS(point_origin[axis] + halfsize_[axis]) <= threshold ? true : false;
	}
	//=================================================================================================//
	Vecd AlignedBoxShape::getUpperPeriodic(int axis, const Vecd &point)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		Vecd shift(0);
		shift[axis] -= 2.0 * halfsize_[axis];
		return transformd_.shiftFrameStationToBase(point_origin + shift);
	}
	//=================================================================================================//
	Vecd AlignedBoxShape::getLowerPeriodic(int axis, const Vecd &point)
	{
		Vecd point_origin = transformd_.shiftBaseStationToFrame(point);
		Vecd shift(0);
		shift[axis] += 2.0 * halfsize_[axis];
		return transformd_.shiftFrameStationToBase(point_origin + shift);
	}
	//=================================================================================================//
}