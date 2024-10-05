#include "complex_shape.h"

namespace SPH
{
//=================================================================================================//
bool AlignedBoxShape::checkInBounds(const Vecd &probe_point, Real lower_bound_fringe, Real upper_bound_fringe)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[alignment_axis_] >= -halfsize_[alignment_axis_] - lower_bound_fringe &&
                   position_in_frame[alignment_axis_] <= halfsize_[alignment_axis_] + upper_bound_fringe
               ? true
               : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkUpperBound(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[alignment_axis_] > halfsize_[alignment_axis_] ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkLowerBound(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return position_in_frame[alignment_axis_] < -halfsize_[alignment_axis_] ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkNearUpperBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[alignment_axis_] - halfsize_[alignment_axis_]) <= threshold ? true : false;
}
//=================================================================================================//
bool AlignedBoxShape::checkNearLowerBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[alignment_axis_] + halfsize_[alignment_axis_]) <= threshold ? true : false;
}
//=================================================================================================//
Vecd AlignedBoxShape::getUpperPeriodic(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[alignment_axis_] -= 2.0 * halfsize_[alignment_axis_];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
Vecd AlignedBoxShape::getLowerPeriodic(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[alignment_axis_] += 2.0 * halfsize_[alignment_axis_];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
} // namespace SPH