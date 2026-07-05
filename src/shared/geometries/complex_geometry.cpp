#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
bool OrientedBox::checkNotFar(const Vecd &probe_point, Real threshold)
{
    return checkContain(probe_point) || checkNearSurface(probe_point, threshold);
}
//=================================================================================================//
bool OrientedBox::checkNearSurface(const Vecd &probe_point, Real threshold)
{
    Vecd distance = probe_point - findClosestPoint(probe_point);
    return distance.cwiseAbs().maxCoeff() < threshold;
}
//=================================================================================================//
bool OrientedBox::checkNearUpperBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[axis_ref_] - halfsize_[axis_ref_]) <= threshold;
}
//=================================================================================================//
bool OrientedBox::checkNearLowerBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[axis_ref_] + halfsize_[axis_ref_]) <= threshold;
}
//=================================================================================================//
Vecd OrientedBox::getLowerPeriodic(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[axis_ref_] += 2.0 * halfsize_[axis_ref_];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
} // namespace SPH