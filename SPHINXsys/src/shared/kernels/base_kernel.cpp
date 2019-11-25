/**
 * @file 	base_kernel.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_kernel.h"

namespace SPH
{
	//===========================================================//
	Kernel::Kernel(Real h, string kernel_name)
		: h_(h), inv_h_(1.0/h), kernel_size_(2.0), 
		kernel_name_(kernel_name)
	{
		if (h <= 0.0)
		{
			std::cout << "\n FAILURE: The Kernel gets a non-positive smoothing length \""
				<< h << "\"!\n";
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	};
	//===========================================================//
	void Kernel::SetDerivativeFactors()
	{
		factor_dW_1D_ = inv_h_ * factor_W_1D_;
		factor_dW_2D_ = inv_h_ * factor_W_2D_;
		factor_dW_3D_ = inv_h_ * factor_W_3D_;
		factor_d2W_1D_ = inv_h_ * factor_dW_1D_;
		factor_d2W_2D_ = inv_h_ * factor_dW_2D_;
		factor_d2W_3D_ = inv_h_ * factor_dW_3D_;
	}
	//===========================================================//
	Real Kernel::W(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;

		return factor_W_1D_ * W_1D(q);
	}
	//===========================================================//
	Real Kernel::W(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		return factor_W_2D_ * W_2D(q);
	}
	//===========================================================//
	Real Kernel::W(const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;

		return factor_W_3D_ * W_3D(q);
	}
	//===========================================================//
	Real Kernel::dW(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;
		return factor_dW_1D_ * dW_1D(q);
	}
	//===========================================================//
	Real Kernel::dW(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		return factor_dW_2D_ * dW_2D(q);
	}
	//===========================================================//
	Real Kernel::dW(const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		return factor_dW_3D_ * dW_3D(q);
	}
	//===========================================================//
	Real Kernel::d2W(const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;
		return factor_d2W_1D_ * d2W_1D(q);
	}
	//===========================================================//
	Real Kernel::d2W(const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		return factor_d2W_2D_ * d2W_2D(q);
	}
	//===========================================================//
	Real Kernel::d2W(const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		return factor_d2W_3D_ * d2W_3D(q);
	}
	//===========================================================//
	Real Kernel::W(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_ * rt_h;
		cout << "\n This function Kernel::W(Real rt_h, const Real& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return factor_W_1D_ * rt_h * W_1D(q);

	}
	//===========================================================//
	Real Kernel::W(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ * rt_h;
		cout << "\n This function Kernel::W(Real rt_h, const Vec2d& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return factor_W_2D_ * rt_h * W_2D(q);

	}
	//===========================================================//
	Real Kernel::W(Real rt_h, const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_ / rt_h;
		cout << "\n This function Kernel::W(Real rt_h, const Vec3d& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return W_3D(q);
	}
	//===========================================================//
	Real Kernel::dW(Real rt_h, const Real& r_ij) const
	{
		Real q = abs(r_ij) * inv_h_;
		cout << "\n This function Kernel::dW(Real rt_h, const Real& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return dW_1D(q); // missing rt_h
	}
	//===========================================================//
	Real Kernel::dW(Real rt_h, const Vec2d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		cout << "\n This function ::dW(Real rt_h, const Vec2d& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return dW_2D(q); // missing rt_h
	}
	//===========================================================//
	Real Kernel::dW(Real rt_h, const Vec3d& r_ij) const
	{
		Real q = r_ij.norm() * inv_h_;
		cout << "\n This function Kernel::dW(Real rt_h, const Vec3d& r_ij) is not done in 3D. Exit the program! \n";
		exit(0);
		return dW_3D(q); // missing rt_h
	}
	//===========================================================//
}

