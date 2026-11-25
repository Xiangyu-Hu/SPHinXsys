/**
 * @file 	spdb_geometry.h
 * @brief 	Parametric SPDB workpiece geometry (SI units).
 * @details Workpiece shapes for the SPDB case. The executable is
 *          tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry.
 *          Synthetic elastic material only — NOT 950 C 20Cr13 prediction.
 */
#ifndef SPDB_GEOMETRY_H
#define SPDB_GEOMETRY_H

#include "sphinxsys.h"

namespace SPH
{
//----------------------------------------------------------------------
// Geometry parameters (SI). Loading axis = z.
//----------------------------------------------------------------------
inline Real cylinder_radius = 0.020;          // 40 mm diameter
inline Real cylinder_height = 0.020;          // 20 mm height
inline Real cylinder_halflength = 0.5 * cylinder_height;
inline Real slot_half_xy = 0.005;             // 10 mm square slot half-width
inline Real slot_depth = 0.010;               // 10 mm blind depth from join face
inline Real particle_spacing_ref = 0.001;     // 1 mm coarse debug spacing
inline Real resolution_ref = particle_spacing_ref;

inline Vec3d lower_center(0.0, 0.0, -cylinder_halflength); // z in [-0.020, 0]
inline Vec3d upper_center(0.0, 0.0, +cylinder_halflength); // z in [0, +0.020]
inline Vec3d slot_center(0.0, 0.0, -0.5 * slot_depth);     // z in [-0.010, 0]
inline Vec3d slot_halfsize(slot_half_xy, slot_half_xy, 0.5 * slot_depth);

inline Real theoretical_full_cylinder_volume()
{
    return Pi * cylinder_radius * cylinder_radius * cylinder_height;
}
inline Real theoretical_slot_volume()
{
    return (2.0 * slot_half_xy) * (2.0 * slot_half_xy) * slot_depth; // 1000 mm^3
}
inline Real theoretical_lower_volume(bool with_slot)
{
    return with_slot ? theoretical_full_cylinder_volume() - theoretical_slot_volume()
                     : theoretical_full_cylinder_volume();
}
inline Real theoretical_upper_volume()
{
    return theoretical_full_cylinder_volume();
}

inline BoundingBoxd system_bounds_for(Real dp)
{
    return BoundingBoxd(
        Vec3d(-cylinder_radius - 2.0 * dp,
              -cylinder_radius - 2.0 * dp,
              -cylinder_height - 2.0 * dp),
        Vec3d(+cylinder_radius + 2.0 * dp,
              +cylinder_radius + 2.0 * dp,
              +cylinder_height + 2.0 * dp));
}

inline BoundingBoxd system_domain_bounds = system_bounds_for(particle_spacing_ref);

inline Transform axis_z_transform(const Vec3d &center)
{
    Rotation3d rotation(Eigen::Quaternion<Real>::FromTwoVectors(Vec3d::UnitX(), Vec3d::UnitZ()));
    return Transform(rotation, center);
}

//----------------------------------------------------------------------
// Slot configuration (first round: none / lower only).
//----------------------------------------------------------------------
enum class SlotMode
{
    None = 0,
    LowerOnly = 1,
    Both = 2
};

//----------------------------------------------------------------------
// Material: synthetic elastic verification (source = synthetic).
// NOT claimed as 950 C 20Cr13.
//----------------------------------------------------------------------
inline Real rho0_s = 7750.0;          // proxy steel density order, demonstration
inline Real Youngs_modulus = 2.0e11;  // synthetic room-temp-like E for elastic check
inline Real poisson = 0.30;
inline Real physical_viscosity = 5.0e6;

// Nominal bonding pressure on the initial full circular section.
// F = 10 MPa * π * R². Constant nominal force, not constant true pressure.
inline Real nominal_pressure = 10.0e6;
inline Real nominal_axial_force()
{
    return nominal_pressure * Pi * cylinder_radius * cylinder_radius;
}
// Softer than Youngs_modulus so |Δh| is observable at dp = 1 mm.
// |ε_z| = p/E = 0.005. This is still synthetic — not 20Cr13 stiffness.
// A 200 GPa bar at the same force would shorten by ~1 µm, below this spacing.
inline Real Youngs_modulus_axial = 2.0e9;

class LowerWorkpiece : public ComplexShape
{
  public:
    explicit LowerWorkpiece(SlotMode slot_mode, const std::string &shape_name = "LowerWorkpiece")
        : ComplexShape(shape_name)
    {
        add<GeometricShapeCylinder>(axis_z_transform(lower_center), cylinder_radius, cylinder_halflength,
                                    "LowerCylinder");
        if (slot_mode == SlotMode::LowerOnly || slot_mode == SlotMode::Both)
        {
            subtract<GeometricShapeBox>(Transform(slot_center), slot_halfsize, "BlindSlotLower");
        }
    }
};

class UpperWorkpiece : public ComplexShape
{
  public:
    explicit UpperWorkpiece(SlotMode slot_mode, const std::string &shape_name = "UpperWorkpiece")
        : ComplexShape(shape_name)
    {
        add<GeometricShapeCylinder>(axis_z_transform(upper_center), cylinder_radius, cylinder_halflength,
                                    "UpperCylinder");
        if (slot_mode == SlotMode::Both)
        {
            // Upper join-face slot opening downward: z in [0, +0.010]
            Vec3d upper_slot_center(0.0, 0.0, 0.5 * slot_depth);
            subtract<GeometricShapeBox>(Transform(upper_slot_center), slot_halfsize, "BlindSlotUpper");
        }
    }
};

inline Real particle_volume_sum(BaseParticles &particles)
{
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    size_t n = particles.TotalRealParticles();
    Real sum = 0.0;
    for (size_t i = 0; i < n; ++i)
        sum += Vol[i];
    return sum;
}

inline void report_body_metrics(const std::string &name, SolidBody &body, Real theoretical_volume)
{
    BaseParticles &particles = body.getBaseParticles();
    size_t n = particles.TotalRealParticles();
    Real discrete_volume = particle_volume_sum(particles);
    Real err = theoretical_volume > 0.0
                   ? std::abs(discrete_volume - theoretical_volume) / theoretical_volume
                   : 0.0;

    Vecd *pos = particles.ParticlePositions();
    Real z_min = pos[0][2], z_max = pos[0][2];
    Real r_max = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        z_min = std::min(z_min, pos[i][2]);
        z_max = std::max(z_max, pos[i][2]);
        r_max = std::max(r_max, std::sqrt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]));
    }

    std::cout << "=== " << name << " ===\n"
              << "  particles: " << n << "\n"
              << "  theoretical_volume_m3: " << theoretical_volume << "\n"
              << "  discrete_volume_m3: " << discrete_volume << "\n"
              << "  relative_volume_error: " << err << "\n"
              << "  z_min_m: " << z_min << "  z_max_m: " << z_max
              << "  height_m: " << (z_max - z_min) << "\n"
              << "  r_max_m: " << r_max << "\n";
}

} // namespace SPH

#endif // SPDB_GEOMETRY_H
