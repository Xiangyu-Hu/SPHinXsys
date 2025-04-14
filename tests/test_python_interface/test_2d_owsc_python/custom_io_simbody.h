#ifndef CUSTOM_IO_SIMBODY_H
#define CUSTOM_IO_SIMBODY_H

#include "io_simbody.h"
#include <SimTKsimbody.h>

namespace SPH {
/**
 * @class WriteSimBodyPinDataExtended
 * @brief Extended class to write total force acting on a solid body and get angles to Python.
 */
class WriteSimBodyPinDataExtended : public WriteSimBodyPinData
{
public:
    WriteSimBodyPinDataExtended(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                                SimTK::MobilizedBody::Pin &pinbody)
        : WriteSimBodyPinData(sph_system, integ, pinbody){};

    // Function to get angle
    Real getAngleToPython(size_t iteration_step = 0) 
    {
        const SimTK::State& state = integ_.getState();
        Real angle = mobody_.getAngle(state);
        return angle;
    }

    // Function to get angle rate
    Real getAngleRateToPython(size_t iteration_step = 0) 
    {
        const SimTK::State& state = integ_.getState();
        Real angle_rate = mobody_.getRate(state);
        return angle_rate;
    }
};
} // namespace SPH
#endif // CUSTOM_IO_SIMBODY_H
