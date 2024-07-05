#pragma once
#include "custom_io_plot.h"

namespace SPH
{
	void CustomPltEngine::writeAQuantityHeader(std::ofstream& out_file, const Matd& quantity, const std::string& quantity_name)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				out_file << "\"" << quantity_name << "[" << i << "][" << j << "]\"" << " ";
			}
		}
	}
	
	
	
	void CustomPltEngine::writeAQuantity(std::ofstream& out_file, const Matd& quantity)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				out_file << std::fixed << std::setprecision(9) << quantity(i, j) << " ";
			}
		}
	}
}
