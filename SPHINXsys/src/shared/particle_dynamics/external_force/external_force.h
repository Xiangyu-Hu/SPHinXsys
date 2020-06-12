/**
 * @file 	external_force.h
 * @brief 	Here, we define the base external foce class.
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
	public:
		ExternalForce();
		virtual ~ExternalForce() {};
		/** This function can be used for runtime control of external force. */
		virtual Vecd InducedAcceleration(Vecd& position) = 0;
	};

	/**
	 * @class Gravity
	 * @brief The gravity force, derived class of External force.
	 */
	class Gravity : public ExternalForce
	{
	protected:
		/** global accerlaeration. */
		Vecd global_acceleration_;
		/** global reference position for zero potential. */
		Vecd reference_position_;
	public:
		Gravity(Vecd gravity_vector, Vecd reference_position = Vecd(0));
		virtual ~Gravity() {};

		/** This function can be used for runtime control of external force. */
		virtual Vecd InducedAcceleration(Vecd& position) override;
		/** Compute potential. */
		Real getPotential(Vecd& position);
	};
}