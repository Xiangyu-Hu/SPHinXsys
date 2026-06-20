#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H

#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_parameters.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * French-paper-inspired reduced cold-crucible case (Jacoutot et al. 2008).
 *
 * Literature-based defaults: D=650 mm, f=300 kHz, sigma=16 S/m @ 1473 K, P_joule ~50 kW example.
 * Remaining dimensions (glass height, coil radius/z/loops, dp, wall) are reduced-case assumptions.
 */
struct OphelieFrenchReducedCaseParams
{
    Vecd glass_center = Vecd(0.0, 0.0, 0.25);
    Real glass_radius = 0.325;
    Real glass_half_height = 0.25;
    Real dp = 0.02;
    /** Simbody triangle mesh resolution for TriangleMeshShapeCylinder (French glass). */
    int glass_mesh_resolution = 20;

    OphelieMultiloopCoilSpec coil;
    Real ampere_turns = 8.0;
    bool auto_coil_z = true;
    Real coil_z_inset = 0.05;
    bool coil_z_user_set = false;

    Real frequency_hz = 300.0e3;
    Real sigma_glass = 16.0;
    Real target_joule_power = 50.0e3;

    bool enable_coil_visual = false;
    bool enable_crucible_visual = false;
    Real crucible_wall_thickness = 0.02;
    Real coil_visual_thickness = 0.02;
};

inline Real frenchReducedGlassHeight(const OphelieFrenchReducedCaseParams &french)
{
    return Real(2.0) * french.glass_half_height;
}

inline void applyFrenchReducedDefaults(OphelieParameters &params, OphelieFrenchReducedCaseParams &french)
{
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    french.coil.stack_center = Vecd(french.glass_center[0], french.glass_center[1], 0.0);
    french.coil.loop_radius = 0.40;
    french.coil.z_min = french.coil_z_inset;
    french.coil.z_max = frenchReducedGlassHeight(french) - french.coil_z_inset;
    french.coil.num_loops = 8;
    french.coil.segments_per_loop = 256;
    french.coil.use_cell_centered_loops = true;
    french.ampere_turns = static_cast<Real>(french.coil.num_loops);
    french.coil.current_per_loop = french.ampere_turns / static_cast<Real>(french.coil.num_loops);

    params.glass_center_ = french.glass_center;
    params.coil_center_ = french.glass_center;
    params.frequency_ = french.frequency_hz;
    params.sigma_glass_ = french.sigma_glass;
    params.target_joule_power_ = french.target_joule_power;
    params.current_amplitude_ = french.coil.current_per_loop;
    params.number_of_turns_ = static_cast<Real>(french.coil.num_loops);
    params.softening_length_ = 0.25 * french.dp;
}

inline void syncFrenchReducedCoilCurrentFromAmpereTurns(OphelieFrenchReducedCaseParams &french)
{
    if (french.coil.num_loops == 0)
    {
        french.coil.current_per_loop = 0.0;
        return;
    }
    french.coil.current_per_loop = french.ampere_turns / static_cast<Real>(french.coil.num_loops);
}

inline void refreshFrenchReducedCoilStack(OphelieFrenchReducedCaseParams &french)
{
    french.coil.stack_center = Vecd(french.glass_center[0], french.glass_center[1], 0.0);
    if (french.auto_coil_z && !french.coil_z_user_set)
    {
        const Real z_bottom = french.glass_center[2] - french.glass_half_height;
        const Real z_top = french.glass_center[2] + french.glass_half_height;
        french.coil.z_min = z_bottom + french.coil_z_inset;
        french.coil.z_max = z_top - french.coil_z_inset;
    }
    syncFrenchReducedCoilCurrentFromAmpereTurns(french);
}

/** Solid glass cylinder from Simbody triangle mesh (SphinxSys official relax style). */
class OphelieFrenchReducedGlassCylinderShape : public ComplexShape
{
  public:
    OphelieFrenchReducedGlassCylinderShape(const std::string &shape_name, const Vecd &center, Real radius,
                                           Real half_height, int mesh_resolution = 20)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), radius, half_height, mesh_resolution, center);
    }
};

/** Thin annular cylinder at coil radius; visualization only. */
class OphelieFrenchCoilVisualShape : public ComplexShape
{
  public:
    OphelieFrenchCoilVisualShape(const std::string &shape_name, const OphelieFrenchReducedCaseParams &french)
        : ComplexShape(shape_name)
    {
        const Real z_center = 0.5 * (french.coil.z_min + french.coil.z_max);
        const Real half_height = 0.5 * (french.coil.z_max - french.coil.z_min) + french.coil_visual_thickness;
        const Vecd center(french.glass_center[0], french.glass_center[1], z_center);
        add<GeometricShapeCylinder>(Transform(center), french.coil.loop_radius + french.coil_visual_thickness,
                                    half_height);
        subtract<GeometricShapeCylinder>(Transform(center),
                                         std::max(french.coil.loop_radius - french.coil_visual_thickness, Real(0.01)),
                                         half_height);
    }
};

/** Hollow cylindrical crucible wall; visualization / future thermal BC only. */
class OphelieFrenchCrucibleWallVisualShape : public ComplexShape
{
  public:
    OphelieFrenchCrucibleWallVisualShape(const std::string &shape_name, const OphelieFrenchReducedCaseParams &french)
        : ComplexShape(shape_name)
    {
        const Real outer_radius = french.glass_radius + french.crucible_wall_thickness;
        add<GeometricShapeCylinder>(Transform(french.glass_center), outer_radius, french.glass_half_height);
        subtract<GeometricShapeCylinder>(Transform(french.glass_center), french.glass_radius, french.glass_half_height);
    }
};

inline BoundingBoxd frenchReducedDomainBounds(const OphelieFrenchReducedCaseParams &french, Real boundary_width)
{
    const Real coil_extent = french.coil.loop_radius + french.coil_visual_thickness + boundary_width;
    const Vecd lower(french.glass_center[0] - coil_extent, french.glass_center[1] - coil_extent,
                     french.glass_center[2] - french.glass_half_height - boundary_width);
    const Vecd upper(french.glass_center[0] + coil_extent, french.glass_center[1] + coil_extent,
                     french.glass_center[2] + french.glass_half_height + boundary_width);
    return BoundingBoxd(lower, upper);
}

inline void printFrenchReducedCaseSummary(const OphelieFrenchReducedCaseParams &french)
{
    std::cout << "[ophelie] French-paper-inspired reduced case (NOT exact CAD reconstruction)\n";
    std::cout << "[ophelie] literature_based: D=" << (2.0 * french.glass_radius)
              << "m f=" << french.frequency_hz << "Hz sigma=" << french.sigma_glass
              << " S/m target_P=" << french.target_joule_power << "W\n";
    std::cout << "[ophelie] reduced_assumptions: glass_center=" << french.glass_center.transpose()
              << " glass_radius=" << french.glass_radius << " glass_height=" << frenchReducedGlassHeight(french)
              << " dp=" << french.dp << " glass_mesh_resolution=" << french.glass_mesh_resolution
              << " coil_radius=" << french.coil.loop_radius
              << " coil_loops=" << french.coil.num_loops << " coil_z=[" << french.coil.z_min << "," << french.coil.z_max
              << "] coil_segments=" << french.coil.segments_per_loop
              << " ampere_turns=" << french.ampere_turns << " I_per_loop=" << french.coil.current_per_loop
              << " cell_centered_loops=" << (french.coil.use_cell_centered_loops ? 1 : 0)
              << " coil_visual=" << (french.enable_coil_visual ? 1 : 0)
              << " crucible_visual=" << (french.enable_crucible_visual ? 1 : 0) << std::endl;
}

inline void syncFrenchReducedToParameters(const OphelieFrenchReducedCaseParams &french, OphelieParameters &params)
{
    params.glass_center_ = french.glass_center;
    params.coil_center_ = french.glass_center;
    params.frequency_ = french.frequency_hz;
    params.sigma_glass_ = french.sigma_glass;
    params.target_joule_power_ = french.target_joule_power;
    params.current_amplitude_ = french.coil.current_per_loop;
    params.number_of_turns_ = static_cast<Real>(french.coil.num_loops);
    params.softening_length_ = 0.25 * french.dp;
}

inline bool isFrenchReducedCommandLineOption(const char *arg)
{
    return std::strncmp(arg, "--glass-radius=", 15) == 0 || std::strncmp(arg, "--glass-height=", 15) == 0 ||
           std::strncmp(arg, "--glass-center=", 15) == 0 || std::strncmp(arg, "--dp=", 5) == 0 ||
           std::strncmp(arg, "--glass-mesh-resolution=", 24) == 0 ||
           std::strncmp(arg, "--coil-radius=", 14) == 0 || std::strncmp(arg, "--coil-turns=", 13) == 0 ||
           std::strncmp(arg, "--coil-num-loops=", 17) == 0 || std::strncmp(arg, "--coil-current=", 15) == 0 || std::strncmp(arg, "--coil-z-min=", 13) == 0 ||
           std::strncmp(arg, "--coil-z-max=", 13) == 0 ||
           std::strncmp(arg, "--coil-segments-per-loop=", 25) == 0 || std::strncmp(arg, "--ampere-turns=", 15) == 0 ||
           std::strncmp(arg, "--frequency=", 12) == 0 || std::strcmp(arg, "--coil-visual") == 0 ||
           std::strcmp(arg, "--crucible-visual") == 0 ||
           std::strncmp(arg, "--crucible-wall-thickness=", 26) == 0;
}

inline void applyFrenchReducedCommandLineOption(const char *arg, OphelieFrenchReducedCaseParams &french)
{
    if (std::strncmp(arg, "--glass-radius=", 15) == 0)
    {
        french.glass_radius = static_cast<Real>(std::atof(arg + 15));
    }
    else if (std::strncmp(arg, "--glass-height=", 15) == 0)
    {
        french.glass_half_height = static_cast<Real>(std::atof(arg + 15)) * Real(0.5);
        french.glass_center[2] = french.glass_half_height;
    }
    else if (std::strncmp(arg, "--glass-center=", 15) == 0)
    {
        Real x = 0.0;
        Real y = 0.0;
        Real z = 0.0;
        if (std::sscanf(arg + 15, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            french.glass_center = Vecd(x, y, z);
        }
    }
    else if (std::strncmp(arg, "--dp=", 5) == 0)
    {
        french.dp = static_cast<Real>(std::atof(arg + 5));
    }
    else if (std::strncmp(arg, "--glass-mesh-resolution=", 24) == 0)
    {
        french.glass_mesh_resolution = std::atoi(arg + 24);
        if (french.glass_mesh_resolution < 4)
        {
            french.glass_mesh_resolution = 4;
        }
    }
    else if (std::strncmp(arg, "--coil-radius=", 14) == 0)
    {
        french.coil.loop_radius = static_cast<Real>(std::atof(arg + 14));
    }
    else if (std::strncmp(arg, "--coil-turns=", 13) == 0 || std::strncmp(arg, "--coil-num-loops=", 17) == 0)
    {
        const char *value = std::strncmp(arg, "--coil-turns=", 13) == 0 ? arg + 13 : arg + 17;
        french.coil.num_loops = static_cast<size_t>(std::atoi(value));
    }
    else if (std::strncmp(arg, "--coil-current=", 15) == 0)
    {
        french.coil.current_per_loop = static_cast<Real>(std::atof(arg + 15));
        french.ampere_turns = french.coil.current_per_loop * static_cast<Real>(french.coil.num_loops);
    }
    else if (std::strncmp(arg, "--coil-z-min=", 13) == 0)
    {
        french.coil.z_min = static_cast<Real>(std::atof(arg + 13));
        french.coil_z_user_set = true;
    }
    else if (std::strncmp(arg, "--coil-z-max=", 13) == 0)
    {
        french.coil.z_max = static_cast<Real>(std::atof(arg + 13));
        french.coil_z_user_set = true;
    }
    else if (std::strncmp(arg, "--coil-segments-per-loop=", 25) == 0)
    {
        french.coil.segments_per_loop = static_cast<size_t>(std::atoi(arg + 25));
    }
    else if (std::strncmp(arg, "--ampere-turns=", 15) == 0)
    {
        french.ampere_turns = static_cast<Real>(std::atof(arg + 15));
    }
    else if (std::strncmp(arg, "--frequency=", 12) == 0)
    {
        french.frequency_hz = static_cast<Real>(std::atof(arg + 12));
    }
    else if (std::strcmp(arg, "--coil-visual") == 0)
    {
        french.enable_coil_visual = true;
    }
    else if (std::strcmp(arg, "--crucible-visual") == 0)
    {
        french.enable_crucible_visual = true;
    }
    else if (std::strncmp(arg, "--crucible-wall-thickness=", 26) == 0)
    {
        french.crucible_wall_thickness = static_cast<Real>(std::atof(arg + 26));
    }
}

inline StdVec<std::string> filterFrenchReducedCommandLine(int ac, char *av[], OphelieFrenchReducedCaseParams &french)
{
    StdVec<std::string> filtered_arguments;
    filtered_arguments.emplace_back(av[0]);
    for (int arg_index = 1; arg_index < ac; ++arg_index)
    {
        const char *arg = av[arg_index];
        if (isFrenchReducedCommandLineOption(arg))
        {
            applyFrenchReducedCommandLineOption(arg, french);
            continue;
        }
        filtered_arguments.emplace_back(arg);
    }
    return filtered_arguments;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H
