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
		Vecd exteranl_acceleration_;	/**< external force. */

	public:
		/**
		 * @brief Default constructor.
		 */
		ExternalForce();
		/**
		 * @brief Default destructor.
		 */
		virtual ~ExternalForce() {};
		/**
		 * @brief The function to define acceleration due to external force with default arguments for simple functions.
		 */
		virtual void UpdateAcceleration();
		/**
		 * @brief This function can be used for runtime control of external force.
		 * @param[in] position The input position of a particle.
		 * @param[in] velocity The input velocity of a particle.
		 * @param[in] t The physical time of simulation.
		 */
		virtual Vecd InducedAcceleration(Vecd position, Vecd velocity, Real t = 0.0);
		virtual Vecd InducedAcceleration() { return exteranl_acceleration_; };
	};

	/**
	 * @class Gravity
	 * @brief The gravity force, derived class of External force.
	 */
	class Gravity : public ExternalForce
	{
	public:
		/**
		 * @brief Default constructor.
		 * @param[in] gravity_vector The input vector of gravity force.
		 */
		Gravity(Vecd gravity_vector);
		/**
		 * @brief Default destructor.
		 */
		virtual ~Gravity() {};
	};
}