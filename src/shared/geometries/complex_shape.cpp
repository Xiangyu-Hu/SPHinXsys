#include "complex_shape.h"

namespace SPH
{
//=================================================================================================//
bool AlignedBoxShape::checkInBounds(int axis, const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[axis] >= -halfsize_[axis] && position_in_frame[axis] <= halfsize_[axis]
               ? true
               : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkUpperBound(int axis, const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[axis] > halfsize_[axis] ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkLowerBound(int axis, const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[axis] < -halfsize_[axis] ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkNearUpperBound(int axis, const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[axis] - halfsize_[axis]) <= threshold ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkNearLowerBound(int axis, const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[axis] + halfsize_[axis]) <= threshold ? true : false;
}
//=================================================================================================//
Vecd AlignedBoxShape::getUpperPeriodic(int axis, const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[axis] -= 2.0 * halfsize_[axis];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
Vecd AlignedBoxShape::getLowerPeriodic(int axis, const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[axis] += 2.0 * halfsize_[axis];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
} // namespace SPH