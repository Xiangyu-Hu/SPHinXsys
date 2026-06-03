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

struct OphelieFrenchReducedCaseParams
{
    Vecd glass_center = Vecd(0.325, 0.325, 0.25);
    Real glass_radius = 0.325;
    Real glass_half_height = 0.25;
    OphelieMultiloopCoilSpec coil;
    Real frequency_hz = 300.0e3;
    Real sigma_glass = 16.0;
    Real target_joule_power = 50.0e3;
};

inline void applyFrenchReducedDefaults(OphelieParameters &params, OphelieFrenchReducedCaseParams &french)
{
    french.coil.stack_center = french.glass_center;
    french.coil.loop_radius = 0.40;
    french.coil.z_min = french.glass_center[2] - french.glass_half_height + 0.05;
    french.coil.z_max = french.glass_center[2] + french.glass_half_height - 0.05;
    french.coil.num_loops = 8;
    french.coil.current_per_loop = 1.0;
    french.coil.segments_per_loop = 180;

    params.glass_center_ = french.glass_center;
    params.coil_center_ = french.glass_center;
    params.frequency_ = french.frequency_hz;
    params.sigma_glass_ = french.sigma_glass;
    params.target_joule_power_ = french.target_joule_power;
    params.current_amplitude_ = french.coil.current_per_loop;
    params.number_of_turns_ = static_cast<Real>(french.coil.num_loops);
    params.softening_length_ = 0.25 * 0.05;
}

class OphelieFrenchReducedGlassCylinderShape : public ComplexShape
{
  public:
    OphelieFrenchReducedGlassCylinderShape(const std::string &shape_name, const Vecd &center, Real radius,
                                           Real half_height)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeCylinder>(Transform(center), radius, half_height);
    }
};

inline BoundingBoxd frenchReducedDomainBounds(const OphelieFrenchReducedCaseParams &french, Real boundary_width)
{
    const Real coil_extent = french.coil.loop_radius + boundary_width;
    const Vecd lower(french.glass_center[0] - coil_extent, french.glass_center[1] - coil_extent,
                     french.glass_center[2] - french.glass_half_height - boundary_width);
    const Vecd upper(french.glass_center[0] + coil_extent, french.glass_center[1] + coil_extent,
                     french.glass_center[2] + french.glass_half_height + boundary_width);
    return BoundingBoxd(lower, upper);
}

inline void printFrenchReducedCaseSummary(const OphelieFrenchReducedCaseParams &french)
{
    std::cout << "[ophelie] French reduced case glass_center=" << french.glass_center.transpose()
              << " glass_radius=" << french.glass_radius << " glass_half_height=" << french.glass_half_height
              << " coil_radius=" << french.coil.loop_radius << " coil_loops=" << french.coil.num_loops
              << " coil_z=[" << french.coil.z_min << "," << french.coil.z_max << "]"
              << " I_per_loop=" << french.coil.current_per_loop << " frequency=" << french.frequency_hz
              << " sigma=" << french.sigma_glass << std::endl;
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
}

inline bool isFrenchReducedCommandLineOption(const char *arg)
{
    return std::strncmp(arg, "--glass-radius=", 15) == 0 || std::strncmp(arg, "--glass-height=", 15) == 0 ||
           std::strncmp(arg, "--glass-center=", 15) == 0 || std::strncmp(arg, "--coil-radius=", 14) == 0 ||
           std::strncmp(arg, "--coil-turns=", 13) == 0 || std::strncmp(arg, "--coil-current=", 15) == 0 ||
           std::strncmp(arg, "--frequency=", 12) == 0;
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
    }
    else if (std::strncmp(arg, "--glass-center=", 15) == 0)
    {
        Real x = 0.0;
        Real y = 0.0;
        Real z = 0.0;
        if (std::sscanf(arg + 15, "%lf,%lf,%lf", &x, &y, &z) == 3)
        {
            french.glass_center = Vecd(x, y, z);
            french.coil.stack_center = french.glass_center;
        }
    }
    else if (std::strncmp(arg, "--coil-radius=", 14) == 0)
    {
        french.coil.loop_radius = static_cast<Real>(std::atof(arg + 14));
    }
    else if (std::strncmp(arg, "--coil-turns=", 13) == 0)
    {
        french.coil.num_loops = static_cast<size_t>(std::atoi(arg + 13));
    }
    else if (std::strncmp(arg, "--coil-current=", 15) == 0)
    {
        french.coil.current_per_loop = static_cast<Real>(std::atof(arg + 15));
    }
    else if (std::strncmp(arg, "--frequency=", 12) == 0)
    {
        french.frequency_hz = static_cast<Real>(std::atof(arg + 12));
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

inline void refreshFrenchReducedCoilStack(OphelieFrenchReducedCaseParams &french)
{
    french.coil.stack_center = french.glass_center;
    french.coil.z_min = french.glass_center[2] - french.glass_half_height + 0.05;
    french.coil.z_max = french.glass_center[2] + french.glass_half_height - 0.05;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H
