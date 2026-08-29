#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_STIRRING_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_STIRRING_GEOMETRY_H

#include "electromagnetic_ophelie_french_reduced_geometry.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * French stirring case from CAD STLs (mm units, geometry_scale=0.001 → meters).
 *
 * The CAD was delivered with the melt axis along -y, but every shared French helper
 * (coil stacking, cylinder level sets, Robin/radiation face tagging, phi boundary normals)
 * assumes a +z cylinder. data/rotate_cad_to_z_up.py rewrites the meshes once into the z-up
 * frame, so the defaults below point at the rotated files:
 *
 *   glass-z.stl             — molten glass, z-cylinder R=0.25 m, z in [0, 0.21] m
 *   stirring_paddle_2_z.stl — paddle solid, hanging from the free surface down to z=0.04 m
 *
 * Reload bodies after relax:
 *   GlassBody    — glass-z minus paddle
 *   Rotor        — paddle
 *   WallBoundary — crucible shell (side annulus + bottom slab, open top)
 */
struct OphelieFrenchStirringCaseParams
{
    OphelieFrenchReducedCaseParams french;

    std::string glass_stl_path = "./input/glass-z.stl";
    std::string rotor_stl_path = "./input/stirring_paddle_2_z.stl";
    std::string rotor_bbox_stl_path;

    Real geometry_scale = 0.001;
    Vecd glass_translation = Vecd::Zero();
    Vecd rotor_translation = Vecd::Zero();
    bool auto_center_rotor_in_glass = false;
    bool rotor_translation_user_set = false;

    Vecd rotation_axis = Vecd(0.0, 0.0, 1.0);
    Vecd rotation_center = Vecd(0.0, 0.0, 0.0);
    /** Jacoutot stirring case: 10 rev/min. */
    Real rotation_speed_rad_s = 1.0471975512;

    Real rho_rotor = 8000.0;

    /** Crucible shell used for the no-slip wall; open at the free surface. */
    Real wall_thickness_factor = 4.0; // multiples of dp
    int wall_mesh_resolution = 40;

    /**
     * Inductor stack height, centred on the melt. EREBUS reports 230 mm, i.e. taller than the
     * 210 mm melt, so the default coil_z_inset rule of the reduced case does not apply here.
     */
    Real coil_stack_height = 0.230;

    /**
     * Level-set mesh refinement relative to dp (SPHinXsys `defineBodyLevelSetShape(par_ck, ratio)`).
     * Thin paddle blades need a finer level set than the particle spacing, following
     * sphinxsys-einsimo test_3d_stirring: ratio = dp / (0.1 * blade_thickness).
     */
    Real rotor_level_set_refinement = 5.0;
    Real glass_level_set_refinement = 1.0;

    /** Comma-separated subset of {glass, rotor, wall}, or "all". */
    std::string relax_bodies = "all";

    int relax_steps = 1000;
    int relax_vtp_every = 100;
};

/**
 * Molten-glass physics for the French stirring case.
 *
 * Sources (see docs/ophelie): Jacoutot et al., Chem. Eng. Process. 47 (2008) 449-455 Table 1
 * for the shared glass properties, and the stirred-run values recorded in the project docs
 * (rho0 = 2800, beta = 3e-5, P = 60 kW, N = 10 rev/min). The electrical conductivity follows
 * the Jacoutot thesis (HAL tel-00121350) figure IV.10 for the mechanically stirred Uox2 glass
 * with 0.96% RuO2, which differs from the figure III.12 law used by the natural-convection case.
 */
struct OphelieFrenchStirringPhysicsParams
{
    Real rho0_glass = 2800.0;   // kg/m^3
    Real mu_glass = 4.0;        // Pa.s
    Real cp_glass = 1150.0;     // J/(kg.K)
    Real k_glass = 4.0;         // W/(m.K)
    Real beta_glass = 3.0e-5;   // 1/K
    Real t_initial = 1473.0;    // K
    Real gravity_g = 9.81;      // m/s^2

    /** log10(sigma) = a - b/T, thesis figure IV.10; sigma(1473 K) = 24.7 S/m. */
    Real sigma_log10_a = 4.05726;
    Real sigma_log10_b = 3923.73;
    Real sigma_min = 1.0e-6;
    Real sigma_max = 30.0;

    /** Robin/radiation boundary condition on the melt, French natural-production convention. */
    Real h_side = 300.0;   // W/(m^2.K)
    Real h_bottom = 35.0;  // W/(m^2.K)
    Real h_free = 20.0;    // W/(m^2.K)
    Real emissivity = 0.8;
    Real t_cool = 300.0;        // crucible coolant [K]
    Real t_ambient = 300.0;     // free-surface convection sink [K]
    Real t_rad_ambient = 300.0; // radiative sink [K]
};

inline Real frenchStirringSigmaAtTemperature(const OphelieFrenchStirringPhysicsParams &physics, Real temperature)
{
    const Real sigma = std::pow(Real(10), physics.sigma_log10_a - physics.sigma_log10_b / temperature);
    return std::min(physics.sigma_max, std::max(physics.sigma_min, sigma));
}

inline bool frenchStirringBodySelected(const OphelieFrenchStirringCaseParams &stirring, const std::string &body)
{
    const std::string &selection = stirring.relax_bodies;
    if (selection == "all" || selection == "both")
    {
        // "both" kept for the pre-wall command lines: glass + rotor only.
        return selection == "all" ? true : body != "wall";
    }
    return selection.find(body) != std::string::npos;
}

inline bool frenchStirringRelaxRotor(const OphelieFrenchStirringCaseParams &stirring)
{
    return frenchStirringBodySelected(stirring, "rotor");
}

inline bool frenchStirringRelaxGlass(const OphelieFrenchStirringCaseParams &stirring)
{
    return frenchStirringBodySelected(stirring, "glass");
}

inline bool frenchStirringRelaxWall(const OphelieFrenchStirringCaseParams &stirring)
{
    return frenchStirringBodySelected(stirring, "wall");
}

inline std::string frenchStirringRotorBboxStlPath(const OphelieFrenchStirringCaseParams &stirring)
{
    return stirring.rotor_bbox_stl_path.empty() ? stirring.rotor_stl_path : stirring.rotor_bbox_stl_path;
}

inline Vecd frenchStirringScaledStlCenter(const std::string &stl_path, const Vecd &translation, Real geometry_scale)
{
    TriangleMeshShapeSTL shape(stl_path, translation, geometry_scale, "StlCenterProbe");
    const BoundingBoxd bounds = shape.findBounds();
    return Real(0.5) * (bounds.lower_ + bounds.upper_);
}

inline BoundingBoxd frenchStirringScaledStlBounds(const std::string &stl_path, const Vecd &translation,
                                                  Real geometry_scale, Real pad = Real(0))
{
    TriangleMeshShapeSTL shape(stl_path, translation, geometry_scale, "StlBoundsProbe");
    const BoundingBoxd bounds = shape.findBounds();
    if (pad <= TinyReal)
    {
        return bounds;
    }
    const Vec3d delta = Vec3d::Constant(pad);
    return BoundingBoxd(bounds.lower_ - delta, bounds.upper_ + delta);
}

/**
 * Feed the CAD melt extent into OphelieFrenchReducedCaseParams.
 *
 * The shared EM pipeline and the thermal Robin/radiation face tagging both describe the melt
 * analytically as a z-cylinder (glass_center, glass_radius, glass_half_height), so those fields
 * have to match the STL the fluid particles actually came from.
 */
inline void syncFrenchStirringGlassToFrenchParams(OphelieFrenchStirringCaseParams &stirring)
{
    const BoundingBoxd bounds =
        frenchStirringScaledStlBounds(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale);
    const Vecd center = Real(0.5) * (bounds.lower_ + bounds.upper_);
    stirring.french.glass_radius =
        Real(0.25) * ((bounds.upper_[0] - bounds.lower_[0]) + (bounds.upper_[1] - bounds.lower_[1]));
    stirring.french.glass_half_height = Real(0.5) * (bounds.upper_[2] - bounds.lower_[2]);
    // Faceting makes the STL bbox centre miss the axis by ~1e-5 m; snap so the coil stays coaxial.
    const Real snap = Real(0.01) * stirring.french.glass_radius;
    stirring.french.glass_center =
        Vecd(std::abs(center[0]) < snap ? Real(0) : center[0], std::abs(center[1]) < snap ? Real(0) : center[1],
             center[2]);
    refreshFrenchReducedCoilStack(stirring.french);
    if (!stirring.french.coil_z_user_set && stirring.coil_stack_height > TinyReal)
    {
        const Real half_stack = Real(0.5) * stirring.coil_stack_height;
        stirring.french.coil.z_min = stirring.french.glass_center[2] - half_stack;
        stirring.french.coil.z_max = stirring.french.glass_center[2] + half_stack;
        stirring.french.coil.stack_center =
            Vecd(stirring.french.glass_center[0], stirring.french.glass_center[1], Real(0));
        syncFrenchReducedCoilCurrentFromAmpereTurns(stirring.french);
    }
}

inline void resolveFrenchStirringGeometryPlacement(OphelieFrenchStirringCaseParams &stirring)
{
    if (stirring.rotor_bbox_stl_path.empty())
    {
        stirring.rotor_bbox_stl_path = stirring.rotor_stl_path;
    }
    if (stirring.auto_center_rotor_in_glass && !stirring.rotor_translation_user_set)
    {
        const Vecd glass_center =
            frenchStirringScaledStlCenter(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale);
        const Vecd rotor_center =
            frenchStirringScaledStlCenter(stirring.rotor_stl_path, Vecd::Zero(), stirring.geometry_scale);
        stirring.rotor_translation = glass_center - rotor_center;
    }
    syncFrenchStirringGlassToFrenchParams(stirring);
}

inline void applyFrenchReducedStirringDefaults(OphelieFrenchStirringCaseParams &stirring,
                                               const OphelieFrenchStirringPhysicsParams &physics =
                                                   OphelieFrenchStirringPhysicsParams())
{
    OphelieParameters params;
    applyFrenchReducedDefaults(params, stirring.french);
    // Thinnest paddle features are ~8 mm, so dp must stay well below that.
    stirring.french.dp = 0.005;
    stirring.glass_stl_path = "./input/glass-z.stl";
    stirring.rotor_stl_path = "./input/stirring_paddle_2_z.stl";
    stirring.geometry_scale = 0.001;
    stirring.glass_translation = Vecd::Zero();
    stirring.rotor_translation = Vecd::Zero();
    stirring.auto_center_rotor_in_glass = false;
    stirring.rotor_translation_user_set = false;
    stirring.rotation_axis = Vecd(0.0, 0.0, 1.0);
    stirring.rotation_center = Vecd::Zero();
    stirring.rotation_speed_rad_s = Real(10.0) * Real(2.0 * M_PI) / Real(60.0); // 10 rev/min

    // Induction settings for the stirred run: EREBUS inductor (7 turns, D=570 mm, stack 230 mm),
    // 300 kHz, 60 kW absorbed by the melt, sigma from thesis figure IV.10 at the initial temperature.
    stirring.french.frequency_hz = 300000.0;
    stirring.french.target_joule_power = 60000.0;
    stirring.french.sigma_glass = frenchStirringSigmaAtTemperature(physics, physics.t_initial);
    stirring.french.coil.loop_radius = 0.285;
    stirring.french.coil.num_loops = 7;
    stirring.french.coil.segments_per_loop = 256;
    stirring.french.coil.use_cell_centered_loops = true;
    stirring.french.ampere_turns = static_cast<Real>(stirring.french.coil.num_loops);
    stirring.french.auto_coil_z = true;
    stirring.french.coil_z_user_set = false;

    // The glass level set must resolve the paddle-shaped cavity, not just the crucible wall,
    // so it needs the same refinement as the paddle itself.
    stirring.rotor_level_set_refinement = 5.0;
    stirring.glass_level_set_refinement = 5.0;
    stirring.relax_bodies = "all";
    stirring.relax_steps = 1000;
    stirring.relax_vtp_every = 100;
}

inline void applyFrenchNaturalStirringDefaults(OphelieFrenchStirringCaseParams &stirring)
{
    applyFrenchReducedStirringDefaults(stirring);
}

inline Real frenchStirringWallThickness(const OphelieFrenchStirringCaseParams &stirring)
{
    return stirring.wall_thickness_factor * stirring.french.dp;
}


/** Molten glass = glass-y.stl minus paddle STL. */
class OphelieFrenchStirringGlassFluidShape : public ComplexShape
{
  public:
    OphelieFrenchStirringGlassFluidShape(const std::string &shape_name, const OphelieFrenchStirringCaseParams &stirring)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale,
                                  "GlassFluid");
        subtract<TriangleMeshShapeSTL>(stirring.rotor_stl_path, stirring.rotor_translation, stirring.geometry_scale,
                                       "RotorCutout");
    }
};

class OphelieFrenchStirringRotorShape : public ComplexShape
{
  public:
    explicit OphelieFrenchStirringRotorShape(const std::string &shape_name,
                                             const OphelieFrenchStirringCaseParams &stirring)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stirring.rotor_stl_path, stirring.rotor_translation, stirring.geometry_scale, "Rotor");
    }
};

class OphelieFrenchStirringBoundingBoxCalculator : public TriangleMeshShapeSTL
{
  public:
    OphelieFrenchStirringBoundingBoxCalculator(const std::string &stl_path, const Vec3d &translation, Real scale,
                                               Real resolution_ref)
        : TriangleMeshShapeSTL(stl_path, translation, scale, "BoundingBoxCalc"), resolution_ref_(resolution_ref)
    {
        initializeFromSTLMesh(stl_path, translation, scale);
    }

    using TriangleMeshShapeSTL::findBounds;

    BoundingBoxd calculate()
    {
        const BoundingBoxd b = this->findBounds();
        const Vec3d delta = Vec3d::Constant(8.0 * resolution_ref_);
        return BoundingBoxd(b.lower_ - delta, b.upper_ + delta);
    }

  private:
    Real resolution_ref_;
};

inline BoundingBoxd frenchStirringRelaxDomainBounds(const OphelieFrenchStirringCaseParams &stirring)
{
    const Real pad = 4.0 * stirring.french.dp;
    const BoundingBoxd glass_bounds =
        frenchStirringScaledStlBounds(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale, pad);
    const BoundingBoxd rotor_bounds = frenchStirringScaledStlBounds(
        frenchStirringRotorBboxStlPath(stirring), stirring.rotor_translation, stirring.geometry_scale, pad);

    Vec3d lower = glass_bounds.lower_;
    Vec3d upper = glass_bounds.upper_;
    for (int d = 0; d < 3; ++d)
    {
        lower[d] = std::min(lower[d], rotor_bounds.lower_[d]);
        upper[d] = std::max(upper[d], rotor_bounds.upper_[d]);
    }
    // The crucible shell sits outside the melt radially and below its floor.
    const Real wall = frenchStirringWallThickness(stirring);
    lower -= Vec3d(wall, wall, wall);
    upper += Vec3d(wall, wall, Real(0));
    return BoundingBoxd(lower, upper);
}

inline BoundingBoxd frenchStirringFlowDomainBounds(const OphelieFrenchStirringCaseParams &stirring)
{
    return frenchStirringRelaxDomainBounds(stirring);
}

inline Real frenchStirringReferenceSpeed(const OphelieFrenchStirringCaseParams &stirring)
{
    const BoundingBoxd glass_bounds =
        frenchStirringScaledStlBounds(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale);
    const Real extent_xy =
        std::max(glass_bounds.upper_[0] - glass_bounds.lower_[0], glass_bounds.upper_[1] - glass_bounds.lower_[1]);
    return Real(1.2) * stirring.rotation_speed_rad_s * Real(0.25) * extent_xy;
}

class UpdateSpinningParticlePosition : public LocalDynamics
{
  public:
    explicit UpdateSpinningParticlePosition(SPHBody &sph_body, const Vecd &rotation_center, const Vecd &spin_axis,
                                            Real omega_spin)
        : LocalDynamics(sph_body), rotation_center0_(rotation_center), spin_axis_(spin_axis.normalized()),
          omega_spin_(omega_spin), dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, UpdateSpinningParticlePosition &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)), rotation_center0_(encloser.rotation_center0_),
              spin_axis_(encloser.spin_axis_), omega_spin_(encloser.omega_spin_)
        {
        }

        void update(size_t index_i, Real physical_time)
        {
            const Real theta = omega_spin_ * physical_time;
            const Real c = std::cos(theta);
            const Real s = std::sin(theta);
            const Vecd rel = pos0_[index_i] - rotation_center0_;
            const Vecd rotated = rel * c + spin_axis_.cross(rel) * s + spin_axis_ * (spin_axis_.dot(rel) * (1.0 - c));
            pos_[index_i] = rotation_center0_ + rotated;
        }

      protected:
        Vecd *pos_, *pos0_;
        Vecd rotation_center0_;
        Vecd spin_axis_;
        Real omega_spin_;
    };

  protected:
    Vecd rotation_center0_;
    Vecd spin_axis_;
    Real omega_spin_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos0_;
};

inline bool isFrenchStirringCommandLineOption(const char *arg)
{
    return std::strncmp(arg, "--glass-stl=", 12) == 0 || std::strncmp(arg, "--rotor-stl=", 12) == 0 ||
           std::strncmp(arg, "--rotor-bbox-stl=", 17) == 0 || std::strncmp(arg, "--geometry-scale=", 17) == 0 ||
           std::strncmp(arg, "--glass-translation=", 20) == 0 ||
           std::strncmp(arg, "--rotor-translation=", 20) == 0 || std::strncmp(arg, "--rotation-speed=", 17) == 0 ||
           std::strncmp(arg, "--rotation-axis=", 16) == 0 || std::strncmp(arg, "--rotation-center=", 18) == 0 ||
           std::strncmp(arg, "--relax-steps=", 14) == 0 || std::strncmp(arg, "--relax-vtp-every=", 18) == 0 ||
           std::strncmp(arg, "--rho-rotor=", 12) == 0 || std::strcmp(arg, "--auto-center-rotor") == 0 ||
           std::strcmp(arg, "--no-auto-center-rotor") == 0 || std::strncmp(arg, "--bodies=", 9) == 0 ||
           std::strncmp(arg, "--rotor-level-set-refinement=", 29) == 0 ||
           std::strncmp(arg, "--glass-level-set-refinement=", 29) == 0 ||
           std::strncmp(arg, "--rotation-rpm=", 15) == 0 ||
           std::strncmp(arg, "--wall-thickness-factor=", 24) == 0 ||
           std::strncmp(arg, "--wall-mesh-resolution=", 23) == 0 ||
           std::strncmp(arg, "--coil-stack-height=", 20) == 0;
}

inline void applyFrenchStirringCommandLineOption(const char *arg, OphelieFrenchStirringCaseParams &stirring)
{
    if (std::strncmp(arg, "--glass-stl=", 12) == 0)
    {
        stirring.glass_stl_path = std::string(arg + 12);
    }
    else if (std::strncmp(arg, "--rotor-stl=", 12) == 0)
    {
        stirring.rotor_stl_path = std::string(arg + 12);
        stirring.rotor_bbox_stl_path.clear();
    }
    else if (std::strncmp(arg, "--rotor-bbox-stl=", 17) == 0)
    {
        stirring.rotor_bbox_stl_path = std::string(arg + 17);
    }
    else if (std::strncmp(arg, "--geometry-scale=", 17) == 0)
    {
        stirring.geometry_scale = static_cast<Real>(std::atof(arg + 17));
    }
    else if (std::strncmp(arg, "--glass-translation=", 20) == 0)
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (std::sscanf(arg + 20, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            stirring.glass_translation = Vecd(static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z));
        }
    }
    else if (std::strncmp(arg, "--rotor-translation=", 20) == 0)
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (std::sscanf(arg + 20, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            stirring.rotor_translation = Vecd(static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z));
            stirring.rotor_translation_user_set = true;
            stirring.auto_center_rotor_in_glass = false;
        }
    }
    else if (std::strcmp(arg, "--auto-center-rotor") == 0)
    {
        stirring.auto_center_rotor_in_glass = true;
        stirring.rotor_translation_user_set = false;
    }
    else if (std::strcmp(arg, "--no-auto-center-rotor") == 0)
    {
        stirring.auto_center_rotor_in_glass = false;
        stirring.rotor_translation_user_set = true;
    }
    else if (std::strncmp(arg, "--rotation-speed=", 17) == 0)
    {
        stirring.rotation_speed_rad_s = static_cast<Real>(std::atof(arg + 17));
    }
    else if (std::strncmp(arg, "--rotation-axis=", 16) == 0)
    {
        double x = 0.0;
        double y = 0.0;
        double z = 1.0;
        if (std::sscanf(arg + 16, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            stirring.rotation_axis = Vecd(static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z));
        }
    }
    else if (std::strncmp(arg, "--rotation-center=", 18) == 0)
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (std::sscanf(arg + 18, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            stirring.rotation_center = Vecd(static_cast<Real>(x), static_cast<Real>(y), static_cast<Real>(z));
        }
    }
    else if (std::strncmp(arg, "--relax-steps=", 14) == 0)
    {
        stirring.relax_steps = std::atoi(arg + 14);
    }
    else if (std::strncmp(arg, "--relax-vtp-every=", 18) == 0)
    {
        stirring.relax_vtp_every = std::atoi(arg + 18);
    }
    else if (std::strncmp(arg, "--rho-rotor=", 12) == 0)
    {
        stirring.rho_rotor = static_cast<Real>(std::atof(arg + 12));
    }
    else if (std::strncmp(arg, "--bodies=", 9) == 0)
    {
        stirring.relax_bodies = std::string(arg + 9);
    }
    else if (std::strncmp(arg, "--rotor-level-set-refinement=", 29) == 0)
    {
        stirring.rotor_level_set_refinement = static_cast<Real>(std::atof(arg + 29));
    }
    else if (std::strncmp(arg, "--glass-level-set-refinement=", 29) == 0)
    {
        stirring.glass_level_set_refinement = static_cast<Real>(std::atof(arg + 29));
    }
    else if (std::strncmp(arg, "--rotation-rpm=", 15) == 0)
    {
        stirring.rotation_speed_rad_s =
            static_cast<Real>(std::atof(arg + 15)) * Real(2.0 * M_PI) / Real(60.0);
    }
    else if (std::strncmp(arg, "--wall-thickness-factor=", 24) == 0)
    {
        stirring.wall_thickness_factor = static_cast<Real>(std::atof(arg + 24));
    }
    else if (std::strncmp(arg, "--wall-mesh-resolution=", 23) == 0)
    {
        stirring.wall_mesh_resolution = std::atoi(arg + 23);
    }
    else if (std::strncmp(arg, "--coil-stack-height=", 20) == 0)
    {
        stirring.coil_stack_height = static_cast<Real>(std::atof(arg + 20));
    }
    else if (isFrenchReducedCommandLineOption(arg))
    {
        applyFrenchReducedCommandLineOption(arg, stirring.french);
    }
}

inline StdVec<std::string> filterFrenchStirringCommandLine(int ac, char *av[], OphelieFrenchStirringCaseParams &stirring)
{
    StdVec<std::string> filtered_arguments;
    filtered_arguments.emplace_back(av[0]);
    for (int arg_index = 1; arg_index < ac; ++arg_index)
    {
        const char *arg = av[arg_index];
        if (isFrenchStirringCommandLineOption(arg) || isFrenchReducedCommandLineOption(arg))
        {
            applyFrenchStirringCommandLineOption(arg, stirring);
            continue;
        }
        filtered_arguments.emplace_back(arg);
    }
    resolveFrenchStirringGeometryPlacement(stirring);
    return filtered_arguments;
}

inline void printFrenchStirringCaseSummary(const OphelieFrenchStirringCaseParams &stirring)
{
    const BoundingBoxd glass_bounds =
        frenchStirringScaledStlBounds(stirring.glass_stl_path, stirring.glass_translation, stirring.geometry_scale);
    const BoundingBoxd rotor_bounds = frenchStirringScaledStlBounds(
        stirring.rotor_stl_path, stirring.rotor_translation, stirring.geometry_scale);
    std::cout << "[ophelie][stirring] CAD glass + paddle STLs\n"
              << "  glass_stl=" << stirring.glass_stl_path << " scale=" << stirring.geometry_scale
              << " translation=" << stirring.glass_translation.transpose() << " dp=" << stirring.french.dp << "\n"
              << "  glass_bbox min=(" << glass_bounds.lower_.transpose() << ") max=(" << glass_bounds.upper_.transpose()
              << ")\n"
              << "  rotor_stl=" << stirring.rotor_stl_path << " translation=" << stirring.rotor_translation.transpose()
              << " auto_center=" << (stirring.auto_center_rotor_in_glass ? 1 : 0) << "\n"
              << "  rotor_bbox min=(" << rotor_bounds.lower_.transpose() << ") max=(" << rotor_bounds.upper_.transpose()
              << ")\n"
              << "  rotation: center=" << stirring.rotation_center.transpose()
              << " axis=" << stirring.rotation_axis.transpose() << " omega=" << stirring.rotation_speed_rad_s
              << " rad/s (" << stirring.rotation_speed_rad_s * Real(60.0) / Real(2.0 * M_PI)
              << " rpm) U_ref~" << frenchStirringReferenceSpeed(stirring) << " m/s\n"
              << "  melt as z-cylinder: center=" << stirring.french.glass_center.transpose()
              << " R=" << stirring.french.glass_radius << " half_H=" << stirring.french.glass_half_height
              << " wall_thickness=" << frenchStirringWallThickness(stirring) << "\n"
              << "  induction: f=" << stirring.french.frequency_hz << " Hz turns=" << stirring.french.coil.num_loops
              << " R_coil=" << stirring.french.coil.loop_radius << " z=[" << stirring.french.coil.z_min << ", "
              << stirring.french.coil.z_max << "] sigma=" << stirring.french.sigma_glass
              << " S/m P_target=" << stirring.french.target_joule_power << " W\n";
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_STIRRING_GEOMETRY_H
