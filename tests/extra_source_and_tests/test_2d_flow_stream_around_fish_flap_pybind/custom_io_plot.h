#pragma once
#ifndef CUSTOM_IO_PLT_H
#define CUSTOM_IO_PLT_H
#include "io_base.h"

namespace SPH
{
	class CustomPltEngine {
	public:
		CustomPltEngine() {};
		virtual ~CustomPltEngine() {};

		void writeAQuantityHeader(std::ofstream& out_file, const Matd& quantity, const std::string& quantity_name);
		void writeAQuantity(std::ofstream& out_file, const Matd& quantity);
	};
}

#endif