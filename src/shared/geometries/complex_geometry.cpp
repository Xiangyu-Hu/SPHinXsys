#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
bool AlignedBox::checkNotFar(const Vecd &probe_point, Real threshold)
{
    return checkContain(probe_point) || checkNearSurface(probe_point, threshold) ? true : false;
}
//=================================================================================================//
bool AlignedBox::checkNearSurface(const Vecd &probe_point, Real threshold)
{
    Vecd distance = probe_point - findClosestPoint(probe_point);
    return distance.cwiseAbs().maxCoeff() < threshold ? true : false;
}
//=================================================================================================//
bool AlignedBox::checkNearUpperBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[alignment_axis_] - halfsize_[alignment_axis_]) <= threshold ? true : false;
}
//=================================================================================================//
bool AlignedBox::checkNearLowerBound(const Vecd &probe_point, Real threshold)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    return ABS(position_in_frame[alignment_axis_] + halfsize_[alignment_axis_]) <= threshold ? true : false;
}
//=================================================================================================//
Vecd AlignedBox::getLowerPeriodic(const Vecd &probe_point)
{
    Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
    Vecd shift = Vecd::Zero();
    shift[alignment_axis_] += 2.0 * halfsize_[alignment_axis_];
    return transform_.shiftFrameStationToBase(position_in_frame + shift);
}
//=================================================================================================//
} // namespace SPH