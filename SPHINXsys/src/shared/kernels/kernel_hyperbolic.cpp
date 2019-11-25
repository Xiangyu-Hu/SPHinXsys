/**
 * @file 	kerel_hyperbolic.cpp
 * @author	Luhui Han, Chi ZHang, Yongchuan Yu and Xiangyu Hu
 * @version	0.1
 */

#include "kernel_hyperbolic.h"

#include <cmath>

namespace SPH
{
	///constructor to initialize the data members (auxiliary factors for kernel calculation)
	KernelHyperbolic::KernelHyperbolic(Real h)
		: Kernel(h, "Hyperbolic")
	{
		factor_W_1D_ = inv_h_ / 7.0;
		factor_W_2D_ = inv_h_ * inv_h_ / (3.0 * pi);
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 15.0 / (62.0 * pi);
		SetDerivativeFactors();
	}
	//===========================================================//
	Real KernelHyperbolic::W_1D(const Real q) const
	{
		if (q < 1.0) {
			return (6.0 - 6.0 * q + powern(q, 3));
		}
		else  {
			return powern(2.0 - q, 3);
		}
	}
	//===========================================================//
	Real KernelHyperbolic::W_2D(const Real q) const
	{
		return W_1D(q);
	}
	//===========================================================//
	Real KernelHyperbolic::W_3D(const Real q) const
	{
		return W_1D(q);
	}
	//===========================================================//
	Real KernelHyperbolic::dW_1D(const Real q) const
	{
		if (q < 1.0) {
			return (-6.0 + 3.0 * powern(q, 2));
		}
		else {
			return powern(2.0 - q, 2) * (-1.0);
		}
	}
	//===========================================================//
	Real KernelHyperbolic::dW_2D(const Real q) const
	{
		return dW_1D(q);
	}
	//===========================================================//
	Real KernelHyperbolic::dW_3D(const Real q) const
	{
		return dW_1D(q);
	}
	//===========================================================//
	Real KernelHyperbolic::d2W_1D(const Real q) const
	{
		if (q < 1.0) {
			return 6.0 * q;
		}
		else {
			return 2.0 * (2.0 - q);
		}
	}
	//===========================================================//
	Real KernelHyperbolic::d2W_2D(const Real q) const
	{
		return d2W_1D(q);
	}
	//===========================================================//
	Real KernelHyperbolic::d2W_3D(const Real q) const
	{
		return d2W_1D(q);
	}
	//===========================================================//
}
	

	
	

