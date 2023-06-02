/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	external_force.h
 * @brief 	Here, we define the base external force class.
 * @details The simple derived classes, such as gravity will be defined in applications.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef EXTERNAL_FORCE_H
#define EXTERNAL_FORCE_H

#include <memory>

#include "base_data_package.h"
#include "execution_unit/execution_proxy.hpp"
#include "execution_unit/execution_argument.hpp"

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

    using namespace execution;

    class GravityKernel {
    public:
        template<class GlobalAccelerationType>
        static Vecd InducedAcceleration(Vecd& position, GlobalAccelerationType&& global_acceleration) {
            return global_acceleration;
        }

        Vecd InducedAcceleration(Vecd& position) {
            return InducedAcceleration(position, global_acceleration_accessor[0]);
        }

        void setAccessors(const std::tuple<sycl::accessor<Vecd, 1, sycl::access_mode::read>> &accessors) {
            global_acceleration_accessor = std::get<0>(accessors);
        }

    private:
        sycl::accessor<Vecd, 1, sycl::access_mode::read> global_acceleration_accessor;
    };

	/**
	 * @class Gravity
	 * @brief The gravity force, derived class of External force.
	 */
	class Gravity : public ExternalForce
	{
	protected:
		Vecd global_acceleration_;
		Vecd zero_potential_reference_;

        DeviceVariable<Vecd, sycl::access_mode::read> global_acceleration_device;
        DeviceProxy<Gravity, GravityKernel, decltype(global_acceleration_device)> device_proxy;

	public:
		Gravity(Vecd gravity_vector, Vecd reference_position = Vecd::Zero());
		virtual ~Gravity() {};

		/** This function can be used for runtime control of external force. */
		virtual Vecd InducedAcceleration(Vecd& position) override;
		Real getPotential(Vecd& position);

        auto& getDeviceProxy() {
            return device_proxy;
        }

        using Proxy = decltype(device_proxy);
	};
}
#endif //EXTERNAL_FORCE_H