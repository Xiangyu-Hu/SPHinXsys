/**
* @file kernel_wenland_c2.h
* @brief This is the classe Wenland kernel..
* @details  NThis kernel has compact support of 2h.
* The smoothing length h can be variable when variable h functions are applied.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_kernel.h"

namespace SPH
{
	/**
	 * @class KernelWendlandC2
	 * @brief Kernel WendlandC2
	 */
	class KernelWendlandC2 : public Kernel
	{
	public:

		/** constructor to initialize the data members 
		(auxiliary factors for kernel calculation) */
		KernelWendlandC2(Real h);

		/** Calculates the kernel value for 
		the given distance of two particles */
		virtual Real W(const Real& r_ij) const override;
		virtual Real W(const Vec2d& r_ij) const override;
		virtual Real W(const Vec3d& r_ij) const override;

		/** calculate derivative of W */
		virtual Real dW(const Real& r_ij) const override;
		virtual Real dW(const Vec2d& r_ij) const override;
		virtual Real dW(const Vec3d& r_ij) const override;

		//---------------------------------------------------------------------
		//for variable smoothing lenght
		//---------------------------------------------------------------------
		/** Calculates the kernel value for 
		the given distance of two particles */
		virtual Real W(Real rt_h, const Real& r_ij) const override;
		virtual Real W(Real rt_h, const Vec2d& r_ij) const override;
		virtual Real W(Real rt_h, const Vec3d& r_ij) const override;

		/** calculate derivative of W */
		virtual Real dW(Real rt_h, const Real& r_ij) const override;
		virtual Real dW(Real rt_h, const Vec2d& r_ij) const override;
		virtual Real dW(Real rt_h, const Vec3d& r_ij) const override;
	};
}
