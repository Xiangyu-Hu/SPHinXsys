#include "base_kernel.h"

namespace SPH
{
	/// Kernel abstract base
	///Constructor to initialize data members
	Kernel::Kernel(Real h, Real kernel_size, string kernel_name)
		: h_(h), inv_h_(1.0/h), kernel_size_(kernel_size), 
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
	string Kernel::GetKerenlName() const
	{
		return kernel_name_;
	}
	//===========================================================//
	Real Kernel::GetSmoothingLength() const
	{
		return h_;
	}
	//===========================================================//
	Real Kernel::GetKernelSize() const
	{
		return kernel_size_;
	}
	//===========================================================//
	Real Kernel::GetCutOffRadius() const
	{
		return h_ * kernel_size_;
	}
}

