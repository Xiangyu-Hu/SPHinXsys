/**
 * @file external_force.h
 * @brief Here, we define the base external foce class.
 * @details The simple derived classes, such as gravity will be defined in applicaitons.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"

namespace SPH {
	/**
	 * @class ExternalForce
	 * @brief This class define external forces.
	 */
	class ExternalForce
	{
	protected:
		/** external force. */
		Vecd exteranl_acceleration_;	

	public:
		ExternalForce();
		virtual ~ExternalForce() {};
		/** The function to define acceleration 
		 * due to external force with default arguments for simple functions. */
		virtual void UpdateAcceleration();
		/** This function can be used for runtime control of external force. */
		virtual Vecd InducedAcceleration(Vecd position, Vecd velocity, Real t = 0.0);
		/** simple constant acceleration */
		virtual Vecd InducedAcceleration() { return exteranl_acceleration_; };
	};

	/**
	 * @class Gravity
	 * @brief The gravity force, derived class of External force.
	 */
	class Gravity : public ExternalForce
	{
	public:
		Gravity(Vecd gravity_vector);
		virtual ~Gravity() {};
	};
}