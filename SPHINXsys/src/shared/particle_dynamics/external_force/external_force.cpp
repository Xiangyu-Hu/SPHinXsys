/**
 * @file 	sexternal_froce.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "external_force.h"
//=================================================================================================//
namespace SPH {
//=================================================================================================//
	ExternalForce::ExternalForce()
	{	
	}
//=================================================================================================//
	Gravity::Gravity(Vecd global_acceleration, Vecd reference_position)
		: ExternalForce(), global_acceleration_(global_acceleration),
		reference_position_(reference_position)
	{
	}
	//=================================================================================================//
	Vecd Gravity::InducedAcceleration(Vecd& position)
	{
		return global_acceleration_;
	}
	//=================================================================================================//
	Real Gravity::getPotential(Vecd& position)
	{
		return dot(InducedAcceleration(position), reference_position_ - position);
	}
	//=================================================================================================//
}
//=================================================================================================//