#include "sdf_primitive.h"

namespace SPH
{
//=================================================================================================//
Real SDFBox::operator()(const Vec3d &point) const
{
    Vec3d d = point.cwiseAbs() - halfsize_;
    return d.cwiseMax(Vec3d::Zero()).norm() + SMIN(d.maxCoeff(), Real(0));
}
//=================================================================================================//
void SDFCylinder::setParameters(Real halflength, Real radius)
{
    halflength_ = halflength;
    radius_ = radius;
}
//=================================================================================================//
Real SDFCylinder::operator()(const Vec3d &point) const
{
    Real axial = point[0];
    Real radial_distance = point.tail(2).norm();
    Real dh = ABS(axial) - halflength_;
    Real dr = radial_distance - radius_;
    return SMAX(dh, dr);
}
//=================================================================================================//
BoundingBox3d SDFCylinder::findBounds() const
{
    return BoundingBox3d(Vec3d(halflength_, radius_, radius_));
}
//=================================================================================================//
void SDFCapsule::setParameters(Real halflength, Real radius)
{
    halflength_ = halflength;
    radius_ = radius;
}
//=================================================================================================//
Real SDFCapsule::operator()(const Vec3d &point) const
{
    Real axial = point[0];
    Real radial_distance = point.tail(2).norm();
    if (axial < 0.0)
        return (point - Vec3d(0.0, 0.0, 0.0)).norm() - radius_; // bottom hemisphere
    else if (axial > halflength_)
        return (point - Vec3d(halflength_, 0.0, 0.0)).norm() - radius_; // top hemisphere
    else
        return radial_distance - radius_; // cylindrical part
}
//=================================================================================================//
BoundingBox3d SDFCapsule::findBounds() const
{
    return BoundingBox3d(Vec3d(halflength_ + radius_, radius_, radius_));
}
//=================================================================================================//
void SDFCone::setParameters(Real halfheight, Real radius)
{
    halfheight_ = halfheight;
    radius_ = radius;
}
//=================================================================================================//
Real SDFCone::operator()(const Vec3d &point) const
{
    Vec2d p = Vec2d(Vec2d(point.x(), point.z()).norm() - radius_, point.y() + halfheight_);
    Vec2d e = Vec2d(-radius_, 2.0 * halfheight_);
    Vec2d q = p - e * std::clamp(p.dot(e) / e.squaredNorm(), Real(0), Real(1));
    Real d = q.norm();
    return SMAX(q.x(), q.y()) > 0.0 ? d : -SMIN(d, p.y());
}
//=================================================================================================//
BoundingBox3d SDFCone::findBounds() const
{
    return BoundingBox3d(Vec3d(halfheight_, radius_, radius_));
}
//=================================================================================================//
void SDFCappedCone::setParameters(Real height_, Real radius1, Real radius2)
{
    height_ = height_;
    radius1_ = radius1;
    radius2_ = radius2;
}
//=================================================================================================//
Real SDFCappedCone::operator()(const Vec3d &point) const
{
    Vec2d q = Vec2d(Vec2d(point.x(), point.z()).norm(), point.y());
    Vec2d k1 = Vec2d(radius2_, halfheight_);
    Vec2d k2 = Vec2d(radius2_ - radius1_, 2.0 * halfheight_);
    Vec2d ca = Vec2d(q.x() - SMIN(q.x(), (q.y() < 0.0) ? radius1_ : radius2_), ABS(q.y()) - halfheight_);
    Vec2d cb = q - k1 + k2 * std::clamp((k1 - q).dot(k2) / k2.squaredNorm(), Real(0), Real(1));
    float s = (cb.x() < 0.0 && ca.y() < 0.0) ? -1.0 : 1.0;
    return s * std::sqrt(SMIN(ca.squaredNorm(), cb.squaredNorm()));
}
//=================================================================================================//
BoundingBox3d SDFCappedCone::findBounds() const
{
    Real max_radius = SMAX(radius1_, radius2_);
    return BoundingBox3d(Vec3d(halfheight_, max_radius, max_radius));
}
//=================================================================================================//
} // namespace SPH
