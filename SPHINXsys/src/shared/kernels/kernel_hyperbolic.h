/**
* @file kernel_hyperbolic.h
* @brief Here, we define hperbolic kernel functions.
* @details  Numerical experiments suggests
* the this kernel is more stable than gaussian like kernel due to its spike
* at the origin. However, it is also found that such kernels give bad density
* predictions. Therefore, the application of this kernel should be clarified.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_kernel.h"

namespace SPH
{
	/**
	 * @class KernelHyperbolic
	 * @brief Kernel from Yang el al.
	 */
	class KernelHyperbolic : public Kernel
	{
	public:

		/** constructor to initialize the data members 
		(auxiliary factors for kernel calculation) */
		KernelHyperbolic(Real h);

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
