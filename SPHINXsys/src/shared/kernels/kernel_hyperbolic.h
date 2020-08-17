/**
* @file kernel_hyperbolic.h
* @brief Here, we define hyperbolic kernel functions.
* @details  Numerical experiments suggests
* the this kernel is more stable than gaussian like kernel due to its spike
* at the origin. However, it is also found that such kernels give bad density
* predictions. Therefore, the application of this kernel should be clarified.
* @author	Chi ZHang and Xiangyu Hu
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
