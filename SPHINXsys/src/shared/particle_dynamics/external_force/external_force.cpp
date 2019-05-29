/**
 * @file 	sexternal_froce.cpp
 * @brief 	Defination of functions decleared in external_foce.h
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "external_force.h"

namespace SPH {
	//===============================================================//
	ExternalForce::ExternalForce() : exteranl_acceleration_(0)
	{	
	}
	//===============================================================//
	void ExternalForce::UpdateAcceleration()
	{

	}
	//===============================================================//
	Vecd ExternalForce
		::InducedAcceleration(Vecd position, Vecd velocity, Real t)
	{
		return exteranl_acceleration_;
	}
	//===============================================================//
	Gravity::Gravity(Vecd gravity_vector)
		: ExternalForce() 
	{
		exteranl_acceleration_ = gravity_vector;
	}
}