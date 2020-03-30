/**
* @file 	base_kernel.h
* @brief 	This is the base classes of kernel functions.  Implementation will be
*			implemented in derived classes. The kernal function define the relevance
* 			between two neighboring particles. Basically, the further the two
*			particles, the less relevance they have.     
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_data_package.h"

#include <string>

using namespace std;

namespace SPH
{
	/**
	 * @class Kernel
	 * @brief Abstract base class of a general SPH kernel function which
	 * is a smoothed Dirac delta function,
	 * a kernel function is radial symmetric, and has a scaling factor.
	 * Based on difference data type in 2d or 3d buildings,
	 * the kernel is defined for 2 and 3 dimensions.
	 * The kernel gives value one at the origin.
	 * The naming of kernel function follows the stand SPH literature.
	 * Currently, only constant smoothing length is applied.
	 * Basically, one can assign different kernel for different particle interactions.
	 */
	class Kernel
	{
	protected:
		/** name of the kernel **/
		const string kernel_name_;
		/** reference smoothing length **/
		const Real h_, inv_h_;
		/** non-dimensional size of the kernel **/
		Real kernel_size_;
		/** Normalization factor for the kernel function  **/
		Real factor_W_1D_, factor_W_2D_, factor_W_3D_;
		/** Auxiliary factors for the derivative of kernel function  **/
		Real factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;
		/** Auxiliary factors for the second order derivative of kernel function  **/
		Real factor_d2W_1D_, factor_d2W_2D_, factor_d2W_3D_;

		/** Fractors for derivatives*/
		virtual void SetDerivativeFactors();
		Real getSmoothingLengthFactor1D(Real inv_h_in) const { 
			return h_ * inv_h_in; 
		};
		Real getSmoothingLengthFactor2D(Real inv_h_in) const {
			return h_ * h_ * inv_h_in * inv_h_in; 
		};
		Real getSmoothingLengthFactor3D(Real inv_h_in)  const {
			return h_ * h_ * h_ * inv_h_in * inv_h_in * inv_h_in;
		};
	public:
		/** Constructor **/
		Kernel(Real h, string kernel_name = "kernel");
		/** Base classes with virtual member functions should have a virtual destructor **/
		virtual ~Kernel() {};
		
		/** access esstential information of the kernel **/
		string GetKerenlName() const { return kernel_name_; };
		Real GetSmoothingLength() const { return h_; };
		Real GetKernelSize() const { return kernel_size_; };
		Real GetCutOffRadius() const { return h_ * kernel_size_; };
		Real GetFactorW1D() const { return factor_W_1D_; };
		Real GetFactorW2D() const { return factor_W_2D_; };
		Real GetFactorW3D() const { return factor_W_3D_; };

		/** Calculates the kernel value for the given displacement of two particles
		  * r_ij pointing from particle j to particle i **/
		virtual Real W(const Real& r_ij) const;
		virtual Real W(const Vec2d& r_ij) const;
		virtual Real W(const Vec3d& r_ij) const;

		/** this value could be use to calculate the value of W **/
		virtual Real W_1D(const Real q) const = 0;
		virtual Real W_2D(const Real q) const = 0;
		virtual Real W_3D(const Real q) const = 0;

		/** Calculates the kernel value at the origin **/
		virtual Real W0(const Real& r_i)  const { return factor_W_1D_; };
		virtual Real W0(const Vec2d& r_i) const { return factor_W_2D_; };
		virtual Real W0(const Vec3d& r_i) const { return factor_W_3D_; };

		/** Calculates the kernel derivation for
		  * the given distance of two particles **/
		virtual Real dW(const Real& r_ij) const;
		virtual Real dW(const Vec2d& r_ij) const;
		virtual Real dW(const Vec3d& r_ij) const;

		/** this value could be use to calculate the value of dW **/
		virtual Real dW_1D(const Real q) const = 0;
		virtual Real dW_2D(const Real q) const = 0;
		virtual Real dW_3D(const Real q) const = 0;

		/** Calculates the kernel second order derivation for
		  * the given distance of two particles **/
		virtual Real d2W(const Real& r_ij) const;
		virtual Real d2W(const Vec2d& r_ij) const;
		virtual Real d2W(const Vec3d& r_ij) const;

		/** this value could be use to calculate the value of d2W **/
		virtual Real d2W_1D(const Real q) const = 0;
		virtual Real d2W_2D(const Real q) const = 0;
		virtual Real d2W_3D(const Real q) const = 0;

		/** for variable smoothing lenght
		  * note that we input the inverse of 
		  * the variable smoothing length.
		  */
		Real GetCutOffRadius(Real smoothing_length) const { 
			return smoothing_length * kernel_size_; 
		};

		Real W(Real inv_h_in, const Real& r_ij) const;
		Real W(Real inv_h_in, const Vec2d& r_ij) const;
		Real W(Real inv_h_in, const Vec3d& r_ij) const;

		/** Calculates the kernel value at the origin **/
		Real W0(Real inv_h_in, const Real& r_i)  const;
		Real W0(Real inv_h_in, const Vec2d& r_i) const;
		Real W0(Real inv_h_in, const Vec3d& r_i) const;

		/** Calculates the kernel derivation for
		  * the given distance of two particles **/
		virtual Real dW(Real inv_h_in, const Real& r_ij) const;
		virtual Real dW(Real inv_h_in, const Vec2d& r_ij) const;
		virtual Real dW(Real inv_h_in, const Vec3d& r_ij) const;

		/** Calculates the kernel second order derivation for
		  * the given distance of two particles **/
		virtual Real d2W(Real inv_h_in, const Real& r_ij) const;
		virtual Real d2W(Real inv_h_in, const Vec2d& r_ij) const;
		virtual Real d2W(Real inv_h_in, const Vec3d& r_ij) const;
	};
}
