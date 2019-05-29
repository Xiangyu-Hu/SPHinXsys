#include "kernel_wenland_c2.h"

#include <cmath>

namespace SPH
{
	///constructor to initialize the data members (auxiliary factors for kernel calculation)
	KernelWendlandC2::KernelWendlandC2(Real h)
		: Kernel(h, 2.0, "Wendland2C")
	{
		factor_W_1D_ = inv_h_  * 5.0 / 8.0;
		factor_W_2D_ = inv_h_ * inv_h_ * 7.0 / (4.0 * pi);
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 21.0 / (16.0 * pi);
		factor_dW_1D_ = inv_h_ * factor_W_1D_;
		factor_dW_2D_ = inv_h_ * factor_W_2D_;
		factor_dW_3D_ = inv_h_ * factor_W_3D_;
	}
	//===========================================================//
	Real KernelWendlandC2::W(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return GetFactorW1D(rt_h) * powern(1.0 - 0.5*q, 3) * (1.0 + 1.5 * q);
		}
		else
		{
			return 0.0;
		}

	}
	//===========================================================//
	Real KernelWendlandC2::W(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		
		if ((q >= 0.0) && (q < 2.0)) {
			return factor_W_2D_ * powern(1.0 - 0.5*q, 4) * (1.0 + 2.0 * q);
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::W(const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 2.0)) {
			return factor_W_3D_ * powern(1.0 - 0.5*q, 4) * (1.0 + 2.0 * q);
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::dW(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;

		if ((q >= 0.0) && (q < 2.0)) {
			return -0.75 * factor_dW_1D_  * powern(q - 2.0, 2) * q;
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::dW(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 2.0)) {
			return 0.625 * factor_dW_2D_ * powern(q - 2.0, 3) * q;
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::dW(const Vec3d& r_ij)  const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 2.0)) {
			return 0.625 * factor_dW_3D_ * powern(q - 2.0, 3) * q;
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::W(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;

		if ((q >= 0.0) && (q < 2.0)) {
			return factor_W_1D_ * powern(1.0 - 0.5*q, 3) * (1.0 + 1.5 * q);
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::W(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return GetFactorW2D(rt_h) * powern(1.0 - 0.5*q, 4) * (1.0 + 2.0 * q);
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::W(Real rt_h, const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return GetFactorW3D(rt_h) * powern(1.0 - 0.5*q, 4) * (1.0 + 2.0 * q);
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::dW(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return -0.75 * GetFactorDW1D(rt_h)  * powern(q - 2.0, 2) * q;
		}
		else
		{
			return 0.0;
		}

	}
	//===========================================================//
	Real KernelWendlandC2::dW(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return 0.625 * GetFactorDW2D(rt_h) * powern(q - 2.0, 3) * q;
		}
		else
		{
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelWendlandC2::dW(Real rt_h, const Vec3d& r_ij)  const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 2.0)) {
			return 0.625 * GetFactorDW3D(rt_h) * powern(q - 2.0, 3) * q;
		}
		else
		{
			return 0.0;
		}
	}
}
