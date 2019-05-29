/**
* @file 	base_kernel.h
* @brief 	This is the base classes of kernel functions.  Implementation will be
*			implemented in derived classes. The kernal function define the relevance
* 			between two neigghboring particles. basically, the further the two
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
	 * Based on difference data_type.h in 2d or 3d buildings
	 * the kernel is defined for 2 and 3 dimensional cases
	 * The kernel gives one at the origin.
	 * The naming of kernel data and function follow the stand SPH literature.
	 * Currently, only constant smoothing length is applied.
	 */
	class Kernel
	{

	protected:
		/** name of the kernel **/
		const string kernel_name_;
		/** reference smoothing length **/
		const Real h_, inv_h_;
		/** non-dimensional size of the kernel **/
		const Real kernel_size_;
		/** Normalization factor for the kernel function  **/
		Real factor_W_1D_, factor_W_2D_, factor_W_3D_;
		/** Auxiliary factors for the derivative of kernel function  **/
		Real factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;

		//---------------------------------------------------------------------
		//for variable smoothing lenght
		//using ratio to the reference smoothing length
		//---------------------------------------------------------------------
		virtual Real GetFactorW1D(Real rt_h) const { return factor_W_1D_ / rt_h; };
		virtual Real GetFactorW2D(Real rt_h) const { return factor_W_2D_ / rt_h / rt_h; };
		virtual Real GetFactorW3D(Real rt_h) const { return factor_W_2D_ / rt_h / rt_h / rt_h; };
		virtual Real GetFactorDW1D(Real rt_h) const { return GetFactorW1D(rt_h) / rt_h; };
		virtual Real GetFactorDW2D(Real rt_h) const { return GetFactorW2D(rt_h) / rt_h; };
		virtual Real GetFactorDW3D(Real rt_h) const { return GetFactorW3D(rt_h) / rt_h; };

	public:
		/** Constructor **/
		Kernel(Real h, Real kernel_size, string kernel_name = "kernel");
		/**Base classes with virtual member functions should have a virtual destructor **/
		virtual ~Kernel() {};

		/** access esstential information of the kernel **/
		string GetKerenlName() const;
		Real GetSmoothingLength() const;
		Real GetKernelSize() const;
		Real GetCutOffRadius() const;

		/** Calculates the kernel value for the given distance of two particles
		  * r_ij pointing from particle j to particle i **/
		virtual Real W(const Real& r_ij) const = 0;
		virtual Real W(const Vec2d& r_ij) const = 0;
		virtual Real W(const Vec3d& r_ij) const = 0;

		/** Calculates the kernel value at the origin **/
		virtual Real W0(const Real& r_i) const { return factor_W_1D_; };
		virtual Real W0(const Vec2d& r_i) const { return factor_W_2D_; };
		virtual Real W0(const Vec3d& r_i) const { return factor_W_3D_; };

		/** Calculates the kernel derivation for
		  * the given distance of two particles **/
		virtual Real dW(const Real& r_ij) const = 0;
		virtual Real dW(const Vec2d& r_ij) const = 0;
		virtual Real dW(const Vec3d& r_ij) const = 0;

		/** Calculates the kernel gradient pointing from r to origin
		  * for the given distance of two particles
		  * same form for 2D and 3D comuptation **/
		Vecd GradW(const Vecd& r_ij) const {
			return dW(r_ij) * normalize(r_ij);
		};

		//---------------------------------------------------------------------
		//for variable smoothing lenght
		//---------------------------------------------------------------------
		/** Calculates the kernel value for the given distance of two particles
		  * r_ij pointing from particle j to particle i **/
		virtual Real W(Real rt_h, const Real& r_ij) const = 0;
		virtual Real W(Real rt_h, const Vec2d& r_ij) const = 0;
		virtual Real W(Real rt_h, const Vec3d& r_ij) const = 0;

		/** Calculates the kernel value at the origin **/
		virtual Real W0(Real rt_h, const Real& r_i) const { return GetFactorW1D(rt_h); };
		virtual Real W0(Real rt_h, const Vec2d& r_i) const { return GetFactorW2D(rt_h); };
		virtual Real W0(Real rt_h, const Vec3d& r_i) const { return GetFactorW3D(rt_h); };

		/** Calculates the kernel derivation for
		  * the given distance of two particles **/
		virtual Real dW(Real rt_h, const Real& r_ij) const = 0;
		virtual Real dW(Real rt_h, const Vec2d& r_ij) const = 0;
		virtual Real dW(Real rt_h, const Vec3d& r_ij) const = 0;

		/** Calculates the kernel gradient pointing from r to origin
		  * for the given distance of two particles
		  * same form for 2D and 3D comuptation **/
		Vecd GradW(Real rt_h, const Vecd& r_ij) const {
			return dW(rt_h, r_ij) * normalize(r_ij);
		};
	};
}
