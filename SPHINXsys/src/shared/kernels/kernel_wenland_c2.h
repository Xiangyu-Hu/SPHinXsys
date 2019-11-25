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
		virtual Real W_1D(const Real q) const override;
		virtual Real W_2D(const Real q) const override;
		virtual Real W_3D(const Real q) const override;

		virtual Real dW_1D(const Real q) const override;
		virtual Real dW_2D(const Real q) const override;
		virtual Real dW_3D(const Real q) const override;

		virtual Real d2W_1D(const Real q) const override;
		virtual Real d2W_2D(const Real q) const override;
		virtual Real d2W_3D(const Real q) const override;
	};
}
