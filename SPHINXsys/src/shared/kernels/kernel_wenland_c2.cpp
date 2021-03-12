/**
 * @file 	kernel_wenland.cpp
 * @author	Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
 */
#include "kernel_wenland_c2.h"

#include <cmath>

namespace SPH
{
	//=================================================================================================//
	KernelWendlandC2::KernelWendlandC2()
		: Kernel("Wendland2CKernel") {}
	//=================================================================================================//
	void KernelWendlandC2::setBasicParameters()
	{
		factor_W_1D_ = inv_h_  * 5.0 / 8.0;
		factor_W_2D_ = inv_h_ * inv_h_ * 7.0 / (4.0 * Pi);
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 21.0 / (16.0 * Pi);
	}
	//=================================================================================================//
	Real KernelWendlandC2::W_1D(const Real q) const
	{
			return  powern(1.0 - 0.5*q, 3) * (1.0 + 1.5 * q);
	}
	//=================================================================================================//
	Real KernelWendlandC2::W_2D(const Real q) const
	{
			return  powern(1.0 - 0.5*q, 4) * (1.0 + 2.0 * q);
	}
	//=================================================================================================//
	Real KernelWendlandC2::W_3D(const Real q) const
	{
			return W_2D(q);
	}
	//=================================================================================================//
	Real KernelWendlandC2::dW_1D(const Real q) const
	{
		return -0.75 * powern(q - 2.0, 2) * q;
	}
	//=================================================================================================//
	Real KernelWendlandC2::dW_2D(const Real q) const
	{
		return 0.625 * powern(q - 2.0, 3) * q;
	}
	//=================================================================================================//
	Real KernelWendlandC2::dW_3D(const Real q) const
	{
		return dW_2D(q);
	}
	//=================================================================================================//
	Real KernelWendlandC2::d2W_1D(const Real q) const
	{
		return -0.75 * (q - 2.0) * (3.0 * q - 2.0);
	}
	//=================================================================================================//
	Real KernelWendlandC2::d2W_2D(const Real q) const
	{
		return 1.25 * powern(q - 2.0, 2) * (2.0 * q - 1.0);
	}
	//=================================================================================================//
	Real KernelWendlandC2::d2W_3D(const Real q) const
	{
		return d2W_2D(q);
	}
	//=================================================================================================//
}
