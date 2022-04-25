/**
 * @file 	kernel_laguerre_gauss.cpp
 * @author	Zhentong Wang,Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */
#include "kernel_laguerre_gauss.h"

#include <cmath>

namespace SPH
{
	//=================================================================================================//
	KernelLaguerreGauss::KernelLaguerreGauss()
		: Kernel("LaguerreGauss") {}
	//=================================================================================================//
	void KernelLaguerreGauss::setBasicParameters()
	{
		factor_W_1D_ = inv_h_ * 8.0 / (5.0 * pow(Pi, 0.5));
		factor_W_2D_ = inv_h_ * inv_h_ * 3.0 / Pi;
		factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 8.0 / pow(Pi, 1.5);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::W_1D(const Real q) const
	{
		return (1.0 - powerN(q, 2) + powerN(q, 4) / 6.0)*exp(-powerN(q, 2));
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::W_2D(const Real q) const
	{
		return  W_1D(q);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::W_3D(const Real q) const
	{
		return W_2D(q);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::dW_1D(const Real q) const
	{
		return (-powerN(q, 5) / 3.0 + 8.0*powerN(q, 3) / 3.0 - 4.0*q)*exp(-powerN(q, 2));
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::dW_2D(const Real q) const
	{
		return dW_1D(q);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::dW_3D(const Real q) const
	{
		return dW_2D(q);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::d2W_1D(const Real q) const
	{
		return (2.0*powerN(q, 6) / 3.0 - 7.0*powerN(q, 4) + 16.0*powerN(q, 2) - 4.0)*exp(-powerN(q, 2));
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::d2W_2D(const Real q) const
	{
		return d2W_1D(q);
	}
	//=================================================================================================//
	Real KernelLaguerreGauss::d2W_3D(const Real q) const
	{
		return d2W_2D(q);
	}
	//=================================================================================================//
}
