#include "kernel_hyperbolic.h"

#include <cmath>

namespace SPH
{
	///constructor to initialize the data members (auxiliary factors for kernel calculation)
	KernelHyperbolic::KernelHyperbolic(Real h)
		: Kernel(h, 2.0, "Hyperbolic")
	{
		factor_W_1D_ = inv_h_ / 7.0;
		factor_W_2D_ = inv_h_ * inv_h_ / (3.0 * pi);
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 15.0 / (62.0 * pi);
		factor_dW_1D_ = inv_h_ * factor_W_1D_;
		factor_dW_2D_ = inv_h_ * factor_W_2D_;
		factor_dW_3D_ = inv_h_ * factor_W_3D_;
	}
	//===========================================================//
	Real KernelHyperbolic::W(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_W_1D_ * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_W_1D_ * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_dW_1D_ * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_dW_1D_ * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_W_2D_ * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_W_2D_ * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_dW_2D_ * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_dW_2D_ * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W(const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_W_3D_ * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_W_3D_ * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(const Vec3d& r_ij)  const
	{
		Real q = r_ij.norm() * inv_h_;

		if ((q >= 0.0) && (q < 1.0)) {
			return factor_dW_3D_ * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return factor_dW_3D_ * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorW1D(rt_h) * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorW1D(rt_h) * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorDW1D(rt_h) * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorDW1D(rt_h) * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorW2D(rt_h) * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorW2D(rt_h) * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorDW2D(rt_h) * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorDW2D(rt_h) * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W(Real rt_h, const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorW3D(rt_h) * (6.0 - 6.0 * q + powern(q, 3));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorW3D(rt_h) * powern(2.0 - q, 3);
		}
		else {
			return 0.0;
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW(Real rt_h, const Vec3d& r_ij)  const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;

		if ((q >= 0.0) && (q < 1.0)) {
			return GetFactorDW3D(rt_h) * (-6.0 + 3.0 * powern(q, 2));
		}
		else if ((q >= 1.0) && (q < 2.0)) {
			return GetFactorDW3D(rt_h) * 3.0 * powern(2.0 - q, 2) * (-1.0);
		}
		else {
			return 0.0;
		}
	}
}
