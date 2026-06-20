#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_GEOMETRY_H

#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * TEAM7-like scaffold on a 1 m cube (same scale as AphiTeam7PhysicalDimensions in the A–phi branch).
 * Plate = conducting slab (aluminum-like); coil = thick annulus around z; remainder = air (no particles).
 */
struct OphelieTeam7Geometry
{
    Real body_length_ = 1.0;
    Real body_height_ = 1.0;
    Real body_width_ = 1.0;

    /** Thin conducting plate (cylinder along z). */
    Vecd plate_center_ = Vecd(0.5, 0.5, 0.08);
    Real plate_radius_ = 0.36;
    Real plate_half_thickness_ = 0.02;

    /** Annular coil (outer minus inner cylinder, axis z). */
    Vecd coil_center_ = Vecd(0.5, 0.5, 0.5);
    Real coil_inner_radius_ = 0.12;
    Real coil_outer_radius_ = 0.18;
    Real coil_half_height_ = 0.12;

    Vecd domainMin(Real boundary_width) const
    {
        return Vecd(-boundary_width, -boundary_width, -boundary_width);
    }

    Vecd domainMax(Real boundary_width) const
    {
        return Vecd(body_length_ + boundary_width, body_height_ + boundary_width, body_width_ + boundary_width);
    }

    /** Annular coil volume (m^3); do not use for J0 = N*I/area. */
    Real coilAnnulusVolume() const
    {
        return Pi * (coil_outer_radius_ * coil_outer_radius_ - coil_inner_radius_ * coil_inner_radius_) *
               (2.0 * coil_half_height_);
    }

    /** Radial–axial cross section for azimuthal current: N*I = J_theta * area. */
    Real coilCurrentCrossSectionArea() const
    {
        return (coil_outer_radius_ - coil_inner_radius_) * (2.0 * coil_half_height_);
    }
};

class OphelieTeam7PlateShape : public ComplexShape
{
  public:
    OphelieTeam7PlateShape(const std::string &shape_name, const OphelieTeam7Geometry &geom)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeCylinder>(Transform(geom.plate_center_), geom.plate_radius_, geom.plate_half_thickness_);
    }
};

class OphelieTeam7AnnularCoilShape : public ComplexShape
{
  public:
    OphelieTeam7AnnularCoilShape(const std::string &shape_name, const OphelieTeam7Geometry &geom)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeCylinder>(Transform(geom.coil_center_), geom.coil_outer_radius_, geom.coil_half_height_);
        subtract<GeometricShapeCylinder>(Transform(geom.coil_center_), geom.coil_inner_radius_, geom.coil_half_height_);
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_GEOMETRY_H
