#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_REDUCED_GEOMETRY_H

#include "electromagnetic_ophelie_finite_section_coil.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_parameters.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <iostream>
#include <limits>

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
    /** Doc145: L0 filament | L1 finite_section_prescribed_J | L2 not implemented. */
    OphelieFrenchCoilExcitationLevel coil_excitation = OphelieFrenchCoilExcitationLevel::Filament;
    /** Used when coil_excitation=L1. coil.loop_radius is treated as centreline radius. */
    Real coil_section_radial_width_m = 0.02;
    Real coil_section_axial_height_m = 0.01;
    size_t coil_section_n_radial = 4;
    size_t coil_section_n_axial = 2;
    /** Imaginary part of peak phasor current per turn (Doc151 complex-I). */
    Real coil_current_imag = 0.0;
    bool actual_CEP_turns_verified = false;

    Real frequency_hz = 300.0e3;
    Real sigma_glass = 16.0;
    Real target_joule_power = 50.0e3;

    bool enable_coil_visual = false;
    bool enable_crucible_visual = false;
    Real crucible_wall_thickness = 0.02;
    Real coil_visual_thickness = 0.02;

    /**
     * Doc 66: bottom-at-zero (glass_center.z ≈ half_height) is a canonical EM
     * reference/export contract, not a universal world-frame law. Stirring STL
     * paths clear this; applyFrenchReducedDefaults sets it true.
     */
    bool require_canonical_bottom_at_zero = false;
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
    // Doc151: z_min/z_max are first/last turn centres for French excitation (L0+L1 share).
    french.coil.use_cell_centered_loops = false;
    french.ampere_turns = static_cast<Real>(french.coil.num_loops);
    french.coil.current_per_loop = french.ampere_turns / static_cast<Real>(french.coil.num_loops);
    // Canonical reduced EM profile: melt bottom on z=0 (Doc 62/66).
    french.require_canonical_bottom_at_zero = true;

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

/** Doc 66: reject silent defaults / trailing garbage / non-finite CLI values. */
inline bool parseStrictFiniteDouble(const char *value, double &out)
{
    if (value == nullptr || *value == '\0')
    {
        return false;
    }
    char *end = nullptr;
    errno = 0;
    out = std::strtod(value, &end);
    if (end == value || errno == ERANGE || !std::isfinite(out))
    {
        return false;
    }
    while (*end == ' ' || *end == '\t')
    {
        ++end;
    }
    return *end == '\0';
}

inline bool parseStrictCommaTripleDoubles(const char *value, double out[3])
{
    if (value == nullptr)
    {
        return false;
    }
    const char *p = value;
    for (int i = 0; i < 3; ++i)
    {
        char *end = nullptr;
        errno = 0;
        out[i] = std::strtod(p, &end);
        if (end == p || errno == ERANGE || !std::isfinite(out[i]))
        {
            return false;
        }
        p = end;
        if (i < 2)
        {
            while (*p == ' ' || *p == '\t')
            {
                ++p;
            }
            if (*p != ',')
            {
                return false;
            }
            ++p;
            while (*p == ' ' || *p == '\t')
            {
                ++p;
            }
        }
    }
    while (*p == ' ' || *p == '\t')
    {
        ++p;
    }
    return *p == '\0';
}

[[noreturn]] inline void frenchReducedFatalCli(const char *option, const char *value, const char *why)
{
    std::cerr << "[ophelie][FATAL] invalid " << option << "='" << (value != nullptr ? value : "") << "': " << why
              << "\n";
    std::exit(1);
}

inline Real frenchReducedStrictPositiveReal(const char *option, const char *value)
{
    double parsed = 0.0;
    if (!parseStrictFiniteDouble(value, parsed) || !(parsed > 0.0))
    {
        frenchReducedFatalCli(option, value, "need a finite positive number with no trailing garbage");
    }
    const Real casted = static_cast<Real>(parsed);
    if (!std::isfinite(static_cast<double>(casted)) || !(casted > Real(0)))
    {
        frenchReducedFatalCli(option, value, "overflowed or became non-finite as Real");
    }
    return casted;
}

inline Real frenchReducedStrictFiniteReal(const char *option, const char *value)
{
    double parsed = 0.0;
    if (!parseStrictFiniteDouble(value, parsed))
    {
        frenchReducedFatalCli(option, value, "need a finite number with no trailing garbage");
    }
    const Real casted = static_cast<Real>(parsed);
    if (!std::isfinite(static_cast<double>(casted)))
    {
        frenchReducedFatalCli(option, value, "overflowed to non-finite Real");
    }
    return casted;
}

inline size_t frenchReducedStrictPositiveCount(const char *option, const char *value)
{
    if (value == nullptr || *value == '\0')
    {
        frenchReducedFatalCli(option, value, "need a positive integer");
    }
    char *end = nullptr;
    errno = 0;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || errno == ERANGE || parsed <= 0)
    {
        frenchReducedFatalCli(option, value, "need a positive integer (negatives rejected; no silent size_t wrap)");
    }
    while (*end == ' ' || *end == '\t')
    {
        ++end;
    }
    if (*end != '\0')
    {
        frenchReducedFatalCli(option, value, "trailing garbage after integer");
    }
    return static_cast<size_t>(parsed);
}

inline void validateFrenchReducedFiniteGeometry(const OphelieFrenchReducedCaseParams &french)
{
    const auto bad = [](Real x) { return !std::isfinite(static_cast<double>(x)); };
    if (bad(french.glass_center[0]) || bad(french.glass_center[1]) || bad(french.glass_center[2]) ||
        bad(french.glass_radius) || !(french.glass_radius > Real(0)) || bad(french.glass_half_height) ||
        !(french.glass_half_height > Real(0)) || bad(french.dp) || !(french.dp > Real(0)) ||
        bad(french.coil.loop_radius) || !(french.coil.loop_radius > Real(0)) || bad(french.coil.z_min) ||
        bad(french.coil.z_max) || !(french.coil.z_max > french.coil.z_min) || french.coil.num_loops == 0 ||
        french.coil.segments_per_loop == 0 || bad(french.coil.current_per_loop) || bad(french.ampere_turns) ||
        bad(french.frequency_hz) || !(french.frequency_hz > Real(0)))
    {
        std::cerr << "[ophelie][FATAL] French reduced geometry has non-finite or non-positive dimensions "
                     "(center/radius/half_height/dp/coil/frequency). Rejecting before mesh construction.\n";
        std::exit(1);
    }
}

/**
 * Doc 62/66: catch the historical melt/coil z misalignment on the canonical bottom-at-zero
 * profile. General geometry only checks relative melt/coil extent + finiteness.
 */
inline void validateFrenchReducedGlassCoilFrame(const OphelieFrenchReducedCaseParams &french)
{
    validateFrenchReducedFiniteGeometry(french);

    const Real z_lo = french.glass_center[2] - french.glass_half_height;
    const Real z_hi = french.glass_center[2] + french.glass_half_height;
    if (french.require_canonical_bottom_at_zero)
    {
        const Real expected_center_z = french.glass_half_height;
        const Real dz = std::abs(french.glass_center[2] - expected_center_z);
        const Real tol = std::max(Real(1.0e-4), Real(0.01) * french.glass_half_height);
        if (dz > tol)
        {
            std::cerr << "[ophelie][FATAL] canonical bottom-at-zero profile requires glass_center.z≈half_height; "
                         "got center.z="
                      << french.glass_center[2] << " half_height=" << expected_center_z << " dz=" << dz
                      << " tol=" << tol
                      << " (Doc 62/66). For a rigidly translated melt+coil, clear "
                         "require_canonical_bottom_at_zero and keep auto_coil_z (or translate coil with the melt).\n";
            std::exit(1);
        }
    }

    if (french.coil.z_min < z_lo - french.glass_half_height || french.coil.z_max > z_hi + french.glass_half_height)
    {
        std::cerr << "[ophelie][warning] coil_z=[" << french.coil.z_min << "," << french.coil.z_max
                  << "] far from glass_z=[" << z_lo << "," << z_hi << "] — verify intentional.\n";
    }
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
    validateFrenchReducedGlassCoilFrame(french);
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

/**
 * Stage 4.1 fluid wall (open-top crucible): SPHinXsys cylinder style
 * (same pattern as tests_sycl mixed_poiseuille / taylor_bar):
 *   TriangleMeshShapeCylinder add/subtract → defineBodyLevelSetShape → NormalFromBodyShapeCK.
 * Cup = taller/lower outer cylinder minus melt cavity (no overlapping second add).
 */
class OphelieFrenchNaturalCrucibleWallShape : public ComplexShape
{
  public:
    /**
     * Open-top cup: annular side wall + bottom. `rim_height` extends the side wall
     * above the melt (still no lid). Default 0 keeps the historical flush rim.
     */
    OphelieFrenchNaturalCrucibleWallShape(const std::string &shape_name, const OphelieFrenchReducedCaseParams &french,
                                          Real wall_thickness, int mesh_resolution = 20, Real rim_height = 0.0)
        : ComplexShape(shape_name)
    {
        const Vecd &c = french.glass_center;
        const Real r_in = french.glass_radius;
        const Real r_out = french.glass_radius + wall_thickness;
        const Real half_h = french.glass_half_height;
        const Real rim = std::max(rim_height, Real(0));
        const Vecd outer_center(c[0], c[1], c[2] + Real(0.5) * rim - Real(0.5) * wall_thickness);
        const Real outer_half_h = half_h + Real(0.5) * rim + Real(0.5) * wall_thickness;
        const Vecd inner_center(c[0], c[1], c[2] + Real(0.5) * rim);
        const Real inner_half_h = half_h + Real(0.5) * rim;
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_out, outer_half_h, mesh_resolution, outer_center,
                                       "OuterBoundary");
        subtract<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_in, inner_half_h, mesh_resolution, inner_center,
                                            "InnerBoundary");
    }
};

/**
 * Doc126: closed crucible — annular side + bottom + flat lid (no open free surface).
 * Outer cylinder minus melt cylinder; lid thickness = wall_thickness.
 */
class OphelieFrenchClosedCrucibleWallShape : public ComplexShape
{
  public:
    OphelieFrenchClosedCrucibleWallShape(const std::string &shape_name, const OphelieFrenchReducedCaseParams &french,
                                         Real wall_thickness, int mesh_resolution = 20)
        : ComplexShape(shape_name)
    {
        const Vecd &c = french.glass_center;
        const Real r_in = french.glass_radius;
        const Real r_out = french.glass_radius + wall_thickness;
        const Real half_h = french.glass_half_height;
        const Vecd outer_center(c[0], c[1], c[2]);
        const Real outer_half_h = half_h + wall_thickness;
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_out, outer_half_h, mesh_resolution, outer_center,
                                       "OuterBoundary");
        subtract<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_in, half_h, mesh_resolution, c, "InnerMelt");
    }
};

/** Flat lid disk only (free-slip viscous partner); pressure contact may reuse full closed wall. */
class OphelieFrenchFlatLidWallShape : public ComplexShape
{
  public:
    OphelieFrenchFlatLidWallShape(const std::string &shape_name, const OphelieFrenchReducedCaseParams &french,
                                  Real wall_thickness, int mesh_resolution = 20)
        : ComplexShape(shape_name)
    {
        const Vecd &c = french.glass_center;
        const Real half_h = french.glass_half_height;
        const Real lid_half = Real(0.5) * wall_thickness;
        const Vecd lid_center(c[0], c[1], c[2] + half_h + lid_half);
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), french.glass_radius + wall_thickness, lid_half,
                                       mesh_resolution, lid_center, "FlatLid");
    }
};

inline BoundingBoxd frenchReducedDomainBounds(const OphelieFrenchReducedCaseParams &french, Real boundary_width)
{
    const Real coil_extent = french.coil.loop_radius + french.coil_visual_thickness + boundary_width;
    const Real glass_z_bottom = french.glass_center[2] - french.glass_half_height;
    const Real glass_z_top = french.glass_center[2] + french.glass_half_height;
    const Real z_lower = std::min(glass_z_bottom, french.coil.z_min);
    const Real z_upper = std::max(glass_z_top, french.coil.z_max);
    const Vecd lower(french.glass_center[0] - coil_extent, french.glass_center[1] - coil_extent,
                     z_lower - boundary_width);
    const Vecd upper(french.glass_center[0] + coil_extent, french.glass_center[1] + coil_extent,
                     z_upper + boundary_width);
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
              << " coil_radius_centerline=" << french.coil.loop_radius
              << " coil_loops=" << french.coil.num_loops << " coil_z=[" << french.coil.z_min << "," << french.coil.z_max
              << "] coil_segments=" << french.coil.segments_per_loop
              << " ampere_turns=" << french.ampere_turns << " I_per_loop=" << french.coil.current_per_loop
              << " cell_centered_loops=" << (french.coil.use_cell_centered_loops ? 1 : 0)
              << " coil_excitation=" << ophelieFrenchCoilExcitationLevelName(french.coil_excitation)
              << " section_rxh=" << french.coil_section_radial_width_m << "x" << french.coil_section_axial_height_m
              << " actual_CEP_turns_verified=" << (french.actual_CEP_turns_verified ? 1 : 0)
              << " coil_visual=" << (french.enable_coil_visual ? 1 : 0)
              << " crucible_visual=" << (french.enable_crucible_visual ? 1 : 0) << std::endl;
}

inline void applyFrenchReducedCoilExcitationBiotToGlass(BaseParticles &glass_particles,
                                                       const OphelieGlassFieldNames &glass_names,
                                                       const OphelieFrenchReducedCaseParams &french, Real mu0,
                                                       Real softening_length)
{
    // Doc151: French excitation always uses endpoint turn centres from coil.z_min/z_max
    // (same list for L0 and L1). cell_centered_loops is ignored on this path.
    applyFrenchCoilExcitationBiotToGlass(glass_particles, glass_names, french.coil, french.coil_excitation,
                                         french.coil_section_radial_width_m, french.coil_section_axial_height_m,
                                         french.coil_section_n_radial, french.coil_section_n_axial, mu0,
                                         softening_length, french.coil_current_imag);
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
           std::strncmp(arg, "--coil-num-loops=", 17) == 0 || std::strncmp(arg, "--coil-current-imag=", 20) == 0 ||
           std::strncmp(arg, "--coil-current=", 15) == 0 || std::strncmp(arg, "--coil-z-min=", 13) == 0 ||
           std::strncmp(arg, "--coil-z-max=", 13) == 0 ||
           std::strncmp(arg, "--coil-segments-per-loop=", 25) == 0 || std::strncmp(arg, "--ampere-turns=", 15) == 0 ||
           std::strncmp(arg, "--coil-excitation=", 18) == 0 ||
           std::strncmp(arg, "--coil-section-radial=", 22) == 0 ||
           std::strncmp(arg, "--coil-section-axial=", 21) == 0 ||
           std::strncmp(arg, "--coil-section-nr=", 18) == 0 || std::strncmp(arg, "--coil-section-nz=", 18) == 0 ||
           std::strncmp(arg, "--frequency=", 12) == 0 || std::strcmp(arg, "--coil-visual") == 0 ||
           std::strcmp(arg, "--crucible-visual") == 0 ||
           std::strncmp(arg, "--crucible-wall-thickness=", 26) == 0 ||
           std::strncmp(arg, "--dump-rz-csv=", 14) == 0 ||
           std::strncmp(arg, "--sigma-mode=", 13) == 0 ||
           std::strncmp(arg, "--sigma-field-file=", 19) == 0 ||
           std::strncmp(arg, "--em-result-json=", 17) == 0 ||
           std::strncmp(arg, "--dump-em-particles-csv=", 24) == 0;
}

inline void applyFrenchReducedCommandLineOption(const char *arg, OphelieFrenchReducedCaseParams &french)
{
    if (std::strncmp(arg, "--glass-radius=", 15) == 0)
    {
        french.glass_radius = frenchReducedStrictPositiveReal("--glass-radius", arg + 15);
    }
    else if (std::strncmp(arg, "--glass-height=", 15) == 0)
    {
        const Real height = frenchReducedStrictPositiveReal("--glass-height", arg + 15);
        french.glass_half_height = Real(0.5) * height;
        french.glass_center[2] = french.glass_half_height;
    }
    else if (std::strncmp(arg, "--glass-center=", 15) == 0)
    {
        // Doc 62: parse into double then cast — Real is float under SYCL.
        // Doc 66: require exactly three finite values and full string consumption.
        double xyz[3] = {0.0, 0.0, 0.0};
        if (!parseStrictCommaTripleDoubles(arg + 15, xyz))
        {
            frenchReducedFatalCli("--glass-center", arg + 15,
                                  "need exactly x,y,z finite doubles (no NaN/Inf/trailing garbage); "
                                  "never silently keep the previous center");
        }
        const Real rx = static_cast<Real>(xyz[0]);
        const Real ry = static_cast<Real>(xyz[1]);
        const Real rz = static_cast<Real>(xyz[2]);
        if (!std::isfinite(static_cast<double>(rx)) || !std::isfinite(static_cast<double>(ry)) ||
            !std::isfinite(static_cast<double>(rz)))
        {
            frenchReducedFatalCli("--glass-center", arg + 15, "overflowed to non-finite Real");
        }
        french.glass_center = Vecd(rx, ry, rz);
    }
    else if (std::strncmp(arg, "--dp=", 5) == 0)
    {
        french.dp = frenchReducedStrictPositiveReal("--dp", arg + 5);
    }
    else if (std::strncmp(arg, "--glass-mesh-resolution=", 24) == 0)
    {
        const size_t res = frenchReducedStrictPositiveCount("--glass-mesh-resolution", arg + 24);
        french.glass_mesh_resolution = static_cast<int>(std::max<size_t>(res, 4));
    }
    else if (std::strncmp(arg, "--coil-radius=", 14) == 0)
    {
        french.coil.loop_radius = frenchReducedStrictPositiveReal("--coil-radius", arg + 14);
    }
    else if (std::strncmp(arg, "--coil-turns=", 13) == 0 || std::strncmp(arg, "--coil-num-loops=", 17) == 0)
    {
        const char *opt = std::strncmp(arg, "--coil-turns=", 13) == 0 ? "--coil-turns" : "--coil-num-loops";
        const char *value = std::strncmp(arg, "--coil-turns=", 13) == 0 ? arg + 13 : arg + 17;
        french.coil.num_loops = frenchReducedStrictPositiveCount(opt, value);
    }
    else if (std::strncmp(arg, "--coil-current=", 15) == 0)
    {
        french.coil.current_per_loop = frenchReducedStrictFiniteReal("--coil-current", arg + 15);
        french.ampere_turns = french.coil.current_per_loop * static_cast<Real>(french.coil.num_loops);
    }
    else if (std::strncmp(arg, "--coil-z-min=", 13) == 0)
    {
        french.coil.z_min = frenchReducedStrictFiniteReal("--coil-z-min", arg + 13);
        french.coil_z_user_set = true;
    }
    else if (std::strncmp(arg, "--coil-z-max=", 13) == 0)
    {
        french.coil.z_max = frenchReducedStrictFiniteReal("--coil-z-max", arg + 13);
        french.coil_z_user_set = true;
    }
    else if (std::strncmp(arg, "--coil-segments-per-loop=", 25) == 0)
    {
        french.coil.segments_per_loop = frenchReducedStrictPositiveCount("--coil-segments-per-loop", arg + 25);
    }
    else if (std::strncmp(arg, "--coil-excitation=", 18) == 0)
    {
        french.coil_excitation = parseOphelieFrenchCoilExcitationLevel(std::string(arg + 18));
        if (french.coil_excitation == OphelieFrenchCoilExcitationLevel::SolvedConductor)
        {
            std::cerr << "[ophelie][error] --coil-excitation=solved_conductor (L2) is not implemented. "
                         "Doc151 forbids fallback to L0.\n";
            std::exit(2);
        }
    }
    else if (std::strncmp(arg, "--coil-current-imag=", 20) == 0)
    {
        french.coil_current_imag = frenchReducedStrictFiniteReal("--coil-current-imag", arg + 20);
    }
    else if (std::strncmp(arg, "--coil-section-radial=", 22) == 0)
    {
        french.coil_section_radial_width_m =
            frenchReducedStrictPositiveReal("--coil-section-radial", arg + 22);
    }
    else if (std::strncmp(arg, "--coil-section-axial=", 21) == 0)
    {
        french.coil_section_axial_height_m =
            frenchReducedStrictPositiveReal("--coil-section-axial", arg + 21);
    }
    else if (std::strncmp(arg, "--coil-section-nr=", 18) == 0)
    {
        french.coil_section_n_radial = frenchReducedStrictPositiveCount("--coil-section-nr", arg + 18);
    }
    else if (std::strncmp(arg, "--coil-section-nz=", 18) == 0)
    {
        french.coil_section_n_axial = frenchReducedStrictPositiveCount("--coil-section-nz", arg + 18);
    }
    else if (std::strncmp(arg, "--ampere-turns=", 15) == 0)
    {
        french.ampere_turns = frenchReducedStrictFiniteReal("--ampere-turns", arg + 15);
    }
    else if (std::strncmp(arg, "--frequency=", 12) == 0)
    {
        french.frequency_hz = frenchReducedStrictPositiveReal("--frequency", arg + 12);
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
        french.crucible_wall_thickness = frenchReducedStrictPositiveReal("--crucible-wall-thickness", arg + 26);
    }
    else if (std::strncmp(arg, "--dump-rz-csv=", 14) == 0)
    {
        // Consumed by aind/picard diagnostic mains; ignore here so CLI filtering drops it.
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
