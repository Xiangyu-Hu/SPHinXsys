/**
 * @file electromagnetic_ophelie_finite_section_coil.h
 * @brief Doc145/151: L0 filament / L1 finite-section prescribed J coil excitation.
 *
 * Doc151: display, L0 and L1 share the same turn_center_z_m[]. L2 is not implemented
 * and must hard-fail (no silent fallback to L0).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_FINITE_SECTION_COIL_H
#define ELECTROMAGNETIC_OPHELIE_FINITE_SECTION_COIL_H

#include "electromagnetic_ophelie_current_moment_source.h"
#include "electromagnetic_ophelie_multiloop_source.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** French / Doc145 coil excitation level (distinct from TEAM7 OphelieCoilSourceModel). */
enum class OphelieFrenchCoilExcitationLevel
{
    Filament = 0,                 // L0
    FiniteSectionPrescribedJ = 1, // L1
    SolvedConductor = 2           // L2 — not implemented
};

inline const char *ophelieFrenchCoilExcitationLevelName(OphelieFrenchCoilExcitationLevel level)
{
    switch (level)
    {
    case OphelieFrenchCoilExcitationLevel::Filament:
        return "filament";
    case OphelieFrenchCoilExcitationLevel::FiniteSectionPrescribedJ:
        return "finite_section_prescribed_J";
    case OphelieFrenchCoilExcitationLevel::SolvedConductor:
        return "solved_conductor";
    }
    return "unknown";
}

inline OphelieFrenchCoilExcitationLevel parseOphelieFrenchCoilExcitationLevel(const std::string &name)
{
    if (name == "finite_section" || name == "finite_section_prescribed_J" || name == "L1" || name == "l1")
    {
        return OphelieFrenchCoilExcitationLevel::FiniteSectionPrescribedJ;
    }
    if (name == "solved_conductor" || name == "L2" || name == "l2")
    {
        return OphelieFrenchCoilExcitationLevel::SolvedConductor;
    }
    if (name == "filament" || name == "L0" || name == "l0")
    {
        return OphelieFrenchCoilExcitationLevel::Filament;
    }
    std::cerr << "[ophelie][error] unknown --coil-excitation='" << name
              << "' (expected filament|finite_section|solved_conductor)\n";
    std::exit(2);
}

struct OphelieFiniteSectionCoilSpec
{
    Vecd stack_axis_xy = Vecd::Zero();
    size_t num_turns = 5;
    Real coil_inner_radius_m = 0.285;
    Real conductor_radial_width_m = 0.02;
    Real conductor_axial_height_m = 0.01;
    Real z_envelope_low_m = -0.0225;
    Real z_envelope_high_m = 0.2075;
    /** Doc151 single source of turn centres. If empty, derived from envelope endpoints. */
    StdVec<Real> turn_center_z_m;
    Real current_per_turn_real = 100.0;
    Real current_per_turn_imag = 0.0;
    size_t n_radial_quad = 4;
    size_t n_axial_quad = 2;
    size_t segments_per_subring = 96;
    std::string section_model = "solid_rectangular_effective";
    std::string geometry_id = "A5_assumed";
    bool actual_CEP_turns_verified = false;
};

inline Real ophelieCoilCenterlineRadius(const OphelieFiniteSectionCoilSpec &spec)
{
    return spec.coil_inner_radius_m + Real(0.5) * spec.conductor_radial_width_m;
}

inline Real ophelieCoilTurnSectionArea(const OphelieFiniteSectionCoilSpec &spec)
{
    return spec.conductor_radial_width_m * spec.conductor_axial_height_m;
}

inline void ophelieCoilPrescribedJTheta(const OphelieFiniteSectionCoilSpec &spec, Real &j_real, Real &j_imag)
{
    const Real area = ophelieCoilTurnSectionArea(spec);
    j_real = spec.current_per_turn_real / (area + TinyReal);
    j_imag = spec.current_per_turn_imag / (area + TinyReal);
}

/** Endpoint-uniform centres between z_first and z_last (inclusive). */
inline void ophelieCoilFillEndpointTurnCenters(size_t num_turns, Real z_first, Real z_last, StdVec<Real> &z_centers)
{
    z_centers.clear();
    if (num_turns == 0)
    {
        return;
    }
    if (num_turns == 1)
    {
        z_centers.push_back(Real(0.5) * (z_first + z_last));
        return;
    }
    const Real dz = (z_last - z_first) / static_cast<Real>(num_turns - 1);
    for (size_t k = 0; k < num_turns; ++k)
    {
        z_centers.push_back(z_first + static_cast<Real>(k) * dz);
    }
}

inline void ophelieCoilUniformTurnCentersFromEnvelope(const OphelieFiniteSectionCoilSpec &spec, StdVec<Real> &z_centers)
{
    const Real z_first = spec.z_envelope_low_m + Real(0.5) * spec.conductor_axial_height_m;
    const Real z_last = spec.z_envelope_high_m - Real(0.5) * spec.conductor_axial_height_m;
    ophelieCoilFillEndpointTurnCenters(spec.num_turns, z_first, z_last, z_centers);
}

/** Resolve turn centres: prefer explicit list; else envelope endpoints. */
inline void ophelieCoilResolveTurnCenters(const OphelieFiniteSectionCoilSpec &spec, StdVec<Real> &z_centers)
{
    if (!spec.turn_center_z_m.empty())
    {
        if (spec.turn_center_z_m.size() != spec.num_turns)
        {
            std::cerr << "[ophelie][error] turn_center_z_m size=" << spec.turn_center_z_m.size()
                      << " != num_turns=" << spec.num_turns << "\n";
            std::exit(2);
        }
        z_centers = spec.turn_center_z_m;
        return;
    }
    ophelieCoilUniformTurnCentersFromEnvelope(spec, z_centers);
}

/** Backward-compatible name used by older call sites. */
inline void ophelieCoilUniformTurnCenters(const OphelieFiniteSectionCoilSpec &spec, StdVec<Real> &z_centers)
{
    ophelieCoilResolveTurnCenters(spec, z_centers);
}

inline bool ophelieCoilTurnsOverlapAxially(const OphelieFiniteSectionCoilSpec &spec)
{
    if (spec.num_turns < 2)
    {
        return false;
    }
    StdVec<Real> z_centers;
    ophelieCoilResolveTurnCenters(spec, z_centers);
    const Real pitch = z_centers[1] - z_centers[0];
    return pitch < spec.conductor_axial_height_m - TinyReal;
}

/**
 * Build unit-current geometric filament moments for a circular loop, then scale by I_re/I_im.
 * Avoids pure-imag failure of "scale imag from real" when I_re=0.
 */
inline void buildCircularLoopComplexCurrentMoments(const Vecd &center, Real radius, size_t segments, Real i_real,
                                                   Real i_imag, StdVec<OphelieCurrentMomentSample> &moments)
{
    moments.clear();
    moments.reserve(std::max<size_t>(1, segments));
    const Real dtheta = 2.0 * Pi / static_cast<Real>(std::max<size_t>(1, segments));
    for (size_t i = 0; i < segments; ++i)
    {
        const Real theta = static_cast<Real>(i) * dtheta;
        const Real theta_next = static_cast<Real>(i + 1) * dtheta;
        const Vecd p0(center[0] + radius * std::cos(theta), center[1] + radius * std::sin(theta), center[2]);
        const Vecd p1(center[0] + radius * std::cos(theta_next), center[1] + radius * std::sin(theta_next), center[2]);
        const Vecd midpoint = 0.5 * (p0 + p1);
        const Vecd dl = p1 - p0;
        moments.push_back(makeFilamentCurrentMoment(midpoint, dl, i_real, i_imag));
    }
}

/**
 * L1: Nr×Nz subrings; Σ I_sub = I_turn. Geometry moments independent of current phase.
 */
inline void buildFiniteSectionPrescribedJMoments(const OphelieFiniteSectionCoilSpec &spec,
                                                 StdVec<OphelieCurrentMomentSample> &moments)
{
    moments.clear();
    const size_t nr = std::max<size_t>(1, spec.n_radial_quad);
    const size_t nz = std::max<size_t>(1, spec.n_axial_quad);
    const Real i_cell_real = spec.current_per_turn_real / static_cast<Real>(nr * nz);
    const Real i_cell_imag = spec.current_per_turn_imag / static_cast<Real>(nr * nz);
    const Real r_inner = spec.coil_inner_radius_m;
    const Real dr = spec.conductor_radial_width_m / static_cast<Real>(nr);
    const Real dh = spec.conductor_axial_height_m / static_cast<Real>(nz);
    const size_t segs = std::max<size_t>(8, spec.segments_per_subring);

    StdVec<Real> z_centers;
    ophelieCoilResolveTurnCenters(spec, z_centers);

    for (size_t turn = 0; turn < z_centers.size(); ++turn)
    {
        const Real z0 = z_centers[turn] - Real(0.5) * spec.conductor_axial_height_m;
        for (size_t ir = 0; ir < nr; ++ir)
        {
            const Real r_sub = r_inner + (static_cast<Real>(ir) + Real(0.5)) * dr;
            for (size_t iz = 0; iz < nz; ++iz)
            {
                const Real z_sub = z0 + (static_cast<Real>(iz) + Real(0.5)) * dh;
                const Vecd center(spec.stack_axis_xy[0], spec.stack_axis_xy[1], z_sub);
                StdVec<OphelieCurrentMomentSample> loop_moments;
                buildCircularLoopComplexCurrentMoments(center, r_sub, segs, i_cell_real, i_cell_imag, loop_moments);
                moments.insert(moments.end(), loop_moments.begin(), loop_moments.end());
            }
        }
    }
}

/** L0 centreline filaments on the same turn centres / complex current. */
inline void buildFiniteSectionL0FilamentMoments(const OphelieFiniteSectionCoilSpec &spec,
                                               StdVec<OphelieCurrentMomentSample> &moments,
                                               size_t segments_per_loop = 128)
{
    moments.clear();
    StdVec<Real> z_centers;
    ophelieCoilResolveTurnCenters(spec, z_centers);
    const Real r_c = ophelieCoilCenterlineRadius(spec);
    const size_t segs = std::max<size_t>(8, segments_per_loop);
    for (Real z : z_centers)
    {
        const Vecd center(spec.stack_axis_xy[0], spec.stack_axis_xy[1], z);
        StdVec<OphelieCurrentMomentSample> loop_moments;
        buildCircularLoopComplexCurrentMoments(center, r_c, segs, spec.current_per_turn_real,
                                               spec.current_per_turn_imag, loop_moments);
        moments.insert(moments.end(), loop_moments.begin(), loop_moments.end());
    }
}

inline Real ophelieFilamentLoopCurrentEstimate(const StdVec<OphelieCurrentMomentSample> &moments, size_t begin,
                                              size_t end)
{
    Real sum_moment = 0.0;
    Real sum_dl = 0.0;
    if (end <= begin || end > moments.size())
    {
        return 0.0;
    }
    for (size_t i = begin; i < end; ++i)
    {
        sum_moment += moments[i].moment_real.norm();
        const size_t j = (i + 1 == end) ? begin : (i + 1);
        sum_dl += (moments[j].position - moments[i].position).norm();
    }
    return sum_moment / (sum_dl + TinyReal);
}

inline OphelieMultiloopCoilSpec ophelieFiniteSectionToL0Multiloop(const OphelieFiniteSectionCoilSpec &spec)
{
    OphelieMultiloopCoilSpec m;
    m.stack_center = Vecd(spec.stack_axis_xy[0], spec.stack_axis_xy[1], 0.0);
    m.loop_radius = ophelieCoilCenterlineRadius(spec);
    m.num_loops = spec.num_turns;
    m.current_per_loop = spec.current_per_turn_real;
    m.segments_per_loop = 128;
    m.use_cell_centered_loops = false; // Doc151: endpoint centres
    StdVec<Real> z_centers;
    ophelieCoilResolveTurnCenters(spec, z_centers);
    if (!z_centers.empty())
    {
        m.z_min = z_centers.front();
        m.z_max = z_centers.back();
    }
    return m;
}

/**
 * Build French L0/L1 excitation from a shared finite-section geometry.
 * coil.z_min / coil.z_max are interpreted as first/last turn centres (not cell envelope).
 * L2 hard-fails. Invalid geometry hard-fails (no L0 fallback).
 */
inline OphelieFiniteSectionCoilSpec makeFrenchSharedCoilSpec(const OphelieMultiloopCoilSpec &coil,
                                                            Real section_radial_width_m, Real section_axial_height_m,
                                                            size_t n_radial_quad, size_t n_axial_quad,
                                                            Real current_imag = 0.0)
{
    if (!(section_radial_width_m > TinyReal) || !(section_axial_height_m > TinyReal))
    {
        std::cerr << "[ophelie][error] coil section sizes must be > 0 (radial=" << section_radial_width_m
                  << " axial=" << section_axial_height_m << ")\n";
        std::exit(2);
    }
    if (coil.num_loops == 0)
    {
        std::cerr << "[ophelie][error] coil num_loops == 0\n";
        std::exit(2);
    }
    OphelieFiniteSectionCoilSpec spec;
    spec.stack_axis_xy = coil.stack_center;
    spec.num_turns = coil.num_loops;
    spec.conductor_radial_width_m = section_radial_width_m;
    spec.conductor_axial_height_m = section_axial_height_m;
    spec.coil_inner_radius_m = coil.loop_radius - Real(0.5) * section_radial_width_m;
    if (!(spec.coil_inner_radius_m > TinyReal))
    {
        std::cerr << "[ophelie][error] L1 coil_inner_radius<=0 from centreline=" << coil.loop_radius
                  << " radial_width=" << section_radial_width_m << "\n";
        std::exit(2);
    }
    ophelieCoilFillEndpointTurnCenters(spec.num_turns, coil.z_min, coil.z_max, spec.turn_center_z_m);
    spec.z_envelope_low_m = spec.turn_center_z_m.front() - Real(0.5) * section_axial_height_m;
    spec.z_envelope_high_m = spec.turn_center_z_m.back() + Real(0.5) * section_axial_height_m;
    spec.current_per_turn_real = coil.current_per_loop;
    spec.current_per_turn_imag = current_imag;
    spec.n_radial_quad = std::max<size_t>(1, n_radial_quad);
    spec.n_axial_quad = std::max<size_t>(1, n_axial_quad);
    spec.segments_per_subring = std::max<size_t>(8, coil.segments_per_loop);
    if (ophelieCoilTurnsOverlapAxially(spec))
    {
        std::cerr << "[ophelie][error] L1 turns overlap axially (pitch < axial height). Refuse to continue.\n";
        std::exit(2);
    }
    return spec;
}

inline void buildFrenchCoilExcitationMoments(const OphelieMultiloopCoilSpec &coil,
                                            OphelieFrenchCoilExcitationLevel level, Real section_radial_width_m,
                                            Real section_axial_height_m, size_t n_radial_quad, size_t n_axial_quad,
                                            StdVec<OphelieCurrentMomentSample> &moments, bool quiet = false,
                                            Real current_imag = 0.0)
{
    if (level == OphelieFrenchCoilExcitationLevel::SolvedConductor)
    {
        std::cerr << "[ophelie][error] coil_excitation=solved_conductor (L2) is not implemented. "
                     "Doc151 forbids silent fallback to L0.\n";
        std::exit(2);
    }

    const OphelieFiniteSectionCoilSpec spec = makeFrenchSharedCoilSpec(
        coil, section_radial_width_m, section_axial_height_m, n_radial_quad, n_axial_quad, current_imag);

    if (level == OphelieFrenchCoilExcitationLevel::FiniteSectionPrescribedJ)
    {
        buildFiniteSectionPrescribedJMoments(spec, moments);
        if (!quiet)
        {
            std::cout << "[ophelie] coil_excitation=finite_section_prescribed_J (L1)"
                      << " geometry_id=" << spec.geometry_id << " turns=" << spec.num_turns
                      << " inner_R=" << spec.coil_inner_radius_m << " centerline_R=" << ophelieCoilCenterlineRadius(spec)
                      << " section=" << spec.conductor_radial_width_m << "x" << spec.conductor_axial_height_m
                      << " Nr×Nz=" << spec.n_radial_quad << "x" << spec.n_axial_quad
                      << " z_centers=[" << spec.turn_center_z_m.front() << "," << spec.turn_center_z_m.back() << "]"
                      << " n_moments=" << moments.size() << std::endl;
        }
        return;
    }

    // L0: same turn centres / centreline radius / complex current as L1 parent geometry.
    buildFiniteSectionL0FilamentMoments(spec, moments, coil.segments_per_loop);
    if (!quiet)
    {
        std::cout << "[ophelie] coil_excitation=filament (L0) same_turn_centers_as_L1"
                  << " turns=" << spec.num_turns << " centerline_R=" << ophelieCoilCenterlineRadius(spec)
                  << " z_centers=[" << spec.turn_center_z_m.front() << "," << spec.turn_center_z_m.back() << "]"
                  << " n_moments=" << moments.size() << std::endl;
    }
}

inline void applyFrenchCoilExcitationBiotToGlass(BaseParticles &glass_particles, const OphelieGlassFieldNames &glass_names,
                                                const OphelieMultiloopCoilSpec &coil,
                                                OphelieFrenchCoilExcitationLevel level, Real section_radial_width_m,
                                                Real section_axial_height_m, size_t n_radial_quad, size_t n_axial_quad,
                                                Real mu0, Real softening_length, Real current_imag = 0.0)
{
    StdVec<OphelieCurrentMomentSample> moments;
    buildFrenchCoilExcitationMoments(coil, level, section_radial_width_m, section_axial_height_m, n_radial_quad,
                                     n_axial_quad, moments, /*quiet=*/false, current_imag);
    applyFilamentBiotSavartToGlassHost(glass_particles, glass_names, moments, mu0, softening_length);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FINITE_SECTION_COIL_H
