/**
 * @file 	kernel_cubic_B_spline.cpp
 * @author	Zhentong Wang, Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
 */
#include "kernel_cubic_B_spline.h"

#include <cmath>

namespace SPH
{
	//=================================================================================================//
	KernelCubicBSpline::KernelCubicBSpline()
		: Kernel("CubicBSpline") {}
	//=================================================================================================//
	void KernelCubicBSpline::setBasicParameters()
	{
		factor_W_1D_ = inv_h_ * 2.0 / 3.0;
		factor_W_2D_ = inv_h_ * inv_h_ * 10.0 / (7.0 * Pi);
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ / Pi;
	}
	//=================================================================================================//
	Real KernelCubicBSpline::W_1D(const Real q) const
	{
		if (q < 1.0) 
		{
			return (1.0 - 3.0 * powerN(q, 2) * (1.0 - q / 2.0) / 2.0);
		}
		else 
		{
			return powerN(2.0 - q, 3) / 4.0;
		}
	}
	//=================================================================================================//
	Real KernelCubicBSpline::W_2D(const Real q) const
	{
		return  W_1D(q);
	}
	//=================================================================================================//
	Real KernelCubicBSpline::W_3D(const Real q) const
	{
		return W_2D(q);
	}
	//=================================================================================================//
	Real KernelCubicBSpline::dW_1D(const Real q) const
	{
		if (q < 1.0)
		{
			return (9.0 * powerN(q, 2) / 4.0 - 3.0 * q);
		}
		else
		{
			return (-1.0) * 3.0 * powerN(2.0 - q, 2) / 4.0;
		}
	}
	//=================================================================================================//
	Real KernelCubicBSpline::dW_2D(const Real q) const
	{
		return dW_1D(q);
	}
	//=================================================================================================//
	Real KernelCubicBSpline::dW_3D(const Real q) const
	{
		return dW_2D(q);
	}
	//=================================================================================================//
	Real KernelCubicBSpline::d2W_1D(const Real q) const
	{
		if (q < 1.0) 
		{
			return 9.0 * q / 2.0 - 3.0;
		}
		else {
			return 3.0 * (2.0 - q) / 2.0;
		}
	}
	//=================================================================================================//
	Real KernelCubicBSpline::d2W_2D(const Real q) const
	{
		return d2W_1D(q);
	}
	//=================================================================================================//
	Real KernelCubicBSpline::d2W_3D(const Real q) const
	{
		return d2W_2D(q);
	}
	//=================================================================================================//
}
