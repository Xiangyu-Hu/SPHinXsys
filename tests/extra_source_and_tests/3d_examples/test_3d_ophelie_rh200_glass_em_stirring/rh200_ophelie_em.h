/**
 * @file rh200_ophelie_em.h
 * @brief Step 2: RH200 glass bbox → multiloop coil → OPHELIE complex edge-flux EM solve once.
 */
#ifndef RH200_OPHELIE_EM_H
#define RH200_OPHELIE_EM_H

#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "rh200_excitation_joule_audit.h"
#include "rh200_fake_joule_heat.h"
#include "rh200_joule_heat_grid.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace rh200
{

inline constexpr const char *kRh200CoilGeometryCsv = "./output/rh200_coil_geometry_log.csv";
inline constexpr const char *kRh200EmSolveCsv = "./output/rh200_em_solve_summary.csv";

struct Rh200EmCli
{
    Real target_power = 50000.0;
    Real sigma0 = 16.0;
    Real frequency_hz = 300000.0;
    Real coil_radius_factor = 1.15;
    Real coil_z_margin_low = 0.05;
    Real coil_z_margin_high = 0.05;
    size_t coil_turns = 6;
};

struct Rh200GlassBboxAudit
{
    Vecd bmin = Vecd::Zero();
    Vecd bmax = Vecd::Zero();
    Vecd size = Vecd::Zero();
    Vecd center = Vecd::Zero();
    Real radial_extent_xy = Real(0);
    size_t particle_count = 0;
};

struct Rh200EmSolveResult
{
    Rh200GlassBboxAudit glass_bbox;
    electromagnetics::ophelie::OphelieFrenchEmSolveResult em;
    Real joule_power_w = Real(0);
    Real phi_eq_res_vol = Real(0);
    Real joule_heat_min = Real(0);
    Real joule_heat_max = Real(0);
    Real base_current_per_loop = Real(0);
    Real equivalent_current_scale = Real(1);
    Real em_power_scale_factor = Real(1);
    Rh200GlassEmFieldAuditSnapshot raw_fields;
    Rh200GlassEmFieldAuditSnapshot scaled_fields;
};

inline Rh200GlassBboxAudit measureGlassParticleBbox(BaseParticles &particles)
{
    Rh200GlassBboxAudit audit;
    audit.particle_count = particles.TotalRealParticles();
    if (audit.particle_count == 0)
    {
        return audit;
    }
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    audit.bmin = pos[0];
    audit.bmax = pos[0];
    for (size_t i = 1; i < audit.particle_count; ++i)
    {
        audit.bmin = audit.bmin.cwiseMin(pos[i]);
        audit.bmax = audit.bmax.cwiseMax(pos[i]);
    }
    audit.size = audit.bmax - audit.bmin;
    audit.center = Real(0.5) * (audit.bmin + audit.bmax);
    for (size_t i = 0; i < audit.particle_count; ++i)
    {
        const Vecd rel = pos[i] - audit.center;
        audit.radial_extent_xy = std::max(audit.radial_extent_xy, std::sqrt(rel[0] * rel[0] + rel[1] * rel[1]));
    }
    return audit;
}

inline Rh200GlassBboxAudit measureStlGlassBbox(const std::string &stl_path, Real geometry_scale)
{
    TriangleMeshShapeSTL shape(stl_path, Vec3d::Zero(), geometry_scale, "glass_fluid");
    const BoundingBoxd bounds = shape.findBounds();
    Rh200GlassBboxAudit audit;
    audit.bmin = Vecd(bounds.lower_[0], bounds.lower_[1], bounds.lower_[2]);
    audit.bmax = Vecd(bounds.upper_[0], bounds.upper_[1], bounds.upper_[2]);
    audit.size = audit.bmax - audit.bmin;
    audit.center = Real(0.5) * (audit.bmin + audit.bmax);
    audit.radial_extent_xy = Real(0.5) * std::max(audit.size[0], audit.size[1]);
    return audit;
}

inline void fillFrenchParamsFromRh200Bbox(const Rh200GlassBboxAudit &bbox, Real dp, const Rh200EmCli &em_cli,
                                          electromagnetics::ophelie::OphelieFrenchReducedCaseParams &french,
                                          electromagnetics::ophelie::OphelieParameters &params)
{
    using namespace electromagnetics::ophelie;
    french.dp = dp;
    french.glass_center = bbox.center;
    french.glass_radius = std::max(bbox.radial_extent_xy, Real(1.0e-3));
    french.glass_half_height = std::max(Real(0.5) * bbox.size[2], Real(1.0e-3));
    french.sigma_glass = em_cli.sigma0;
    french.frequency_hz = em_cli.frequency_hz;
    french.target_joule_power = em_cli.target_power;
    french.auto_coil_z = false;
    french.coil_z_user_set = true;

    const Real h_fluid = std::max(bbox.size[2], Real(1.0e-3));
    french.coil.stack_center = Vecd(bbox.center[0], bbox.center[1], Real(0));
    french.coil.loop_radius = em_cli.coil_radius_factor * french.glass_radius;
    french.coil.z_min = bbox.bmin[2] + em_cli.coil_z_margin_low * h_fluid;
    french.coil.z_max = bbox.bmax[2] - em_cli.coil_z_margin_high * h_fluid;
    if (french.coil.z_max <= french.coil.z_min)
    {
        french.coil.z_min = bbox.bmin[2] + Real(0.02) * h_fluid;
        french.coil.z_max = bbox.bmax[2] - Real(0.02) * h_fluid;
    }
    french.coil.num_loops = std::max<size_t>(em_cli.coil_turns, 1);
    french.coil.segments_per_loop = 256;
    french.coil.use_cell_centered_loops = true;
    french.coil.current_per_loop = Real(1);
    french.ampere_turns = static_cast<Real>(french.coil.num_loops);
    syncFrenchReducedCoilCurrentFromAmpereTurns(french);
    syncFrenchReducedToParameters(french, params);
    params.softening_length_ = Real(0.25) * dp;
}

inline void writeRh200CoilGeometryLogCsv(const std::string &path, Real geometry_scale,
                                         const Rh200GlassBboxAudit &stl_bbox, const Rh200GlassBboxAudit &particle_bbox,
                                         const electromagnetics::ophelie::OphelieFrenchReducedCaseParams &french)
{
    std::ofstream out(path);
    out << "label,geometry_scale,cx,cy,cz,Lx,Ly,Lz,radial_extent_xy,particle_count,"
        << "coil_R,coil_z_min,coil_z_max,coil_turns,coil_segments,frequency_hz,sigma,target_power_W\n";
    out << "stl_glass," << geometry_scale << "," << stl_bbox.center[0] << "," << stl_bbox.center[1] << ","
        << stl_bbox.center[2] << "," << stl_bbox.size[0] << "," << stl_bbox.size[1] << "," << stl_bbox.size[2] << ","
        << stl_bbox.radial_extent_xy << ",0," << french.coil.loop_radius << "," << french.coil.z_min << ","
        << french.coil.z_max << "," << french.coil.num_loops << "," << french.coil.segments_per_loop << ","
        << french.frequency_hz << "," << french.sigma_glass << "," << french.target_joule_power << "\n";
    out << "glass_particles," << geometry_scale << "," << particle_bbox.center[0] << "," << particle_bbox.center[1]
        << "," << particle_bbox.center[2] << "," << particle_bbox.size[0] << "," << particle_bbox.size[1] << ","
        << particle_bbox.size[2] << "," << particle_bbox.radial_extent_xy << "," << particle_bbox.particle_count << ","
        << french.coil.loop_radius << "," << french.coil.z_min << "," << french.coil.z_max << ","
        << french.coil.num_loops << "," << french.coil.segments_per_loop << "," << french.frequency_hz << ","
        << french.sigma_glass << "," << french.target_joule_power << "\n";
    std::cout << "[rh200] wrote coil geometry log: " << path << std::endl;
}

inline void appendRh200GlassParticleBboxToAuditCsv(const std::string &path, const Rh200GlassBboxAudit &particle_bbox)
{
    std::ofstream out(path, std::ios::app);
    out << "glass_particles_reload,1," << particle_bbox.bmin[0] << "," << particle_bbox.bmin[1] << ","
        << particle_bbox.bmin[2] << "," << particle_bbox.bmax[0] << "," << particle_bbox.bmax[1] << ","
        << particle_bbox.bmax[2] << "," << particle_bbox.size[0] << "," << particle_bbox.size[1] << ","
        << particle_bbox.size[2] << "," << particle_bbox.center[0] << "," << particle_bbox.center[1] << ","
        << particle_bbox.center[2] << "\n";
}

inline void writeRh200EmSolveSummaryCsv(const std::string &path, const Rh200EmSolveResult &result,
                                        const electromagnetics::ophelie::OphelieFrenchReducedCaseParams &french)
{
    std::ofstream out(path);
    out << "P_joule_W,phi_eq_res_vol,joule_heat_min,joule_heat_max,coil_R,coil_z_min,coil_z_max,coil_turns,"
           "frequency_hz,sigma,target_power_W,n_glass\n";
    out << result.joule_power_w << "," << result.phi_eq_res_vol << "," << result.joule_heat_min << ","
        << result.joule_heat_max << "," << french.coil.loop_radius << "," << french.coil.z_min << ","
        << french.coil.z_max << "," << french.coil.num_loops << "," << french.frequency_hz << "," << french.sigma_glass
        << "," << french.target_joule_power << "," << result.glass_bbox.particle_count << "\n";
    std::cout << "[rh200] wrote EM solve summary: " << path << std::endl;
}

inline void configureRh200OphelieEmParameters(electromagnetics::ophelie::OphelieParameters &params,
                                              electromagnetics::ophelie::OphelieTestCliOptions &cli_options)
{
    using namespace electromagnetics::ophelie;
    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_power_scaling_ = false;
    cli_options.no_power_scaling = true;
    finalizeOphelieCurrentFormConfiguration(params, cli_options);
}

inline void hostExtractEmDepositionFields(
    BaseParticles &particles, const electromagnetics::ophelie::OphelieGlassFieldNames &glass_names,
    Rh200EmParticleDepositionHost &dep)
{
    using namespace electromagnetics::ophelie;
    const size_t n = particles.TotalRealParticles();
    dep.position.resize(n);
    dep.volumetric_measure.resize(n);
    dep.joule_heat.resize(n);
    dep.phi_real.resize(n);
    dep.phi_imag.resize(n);
    dep.e_real.resize(n);
    dep.e_imag.resize(n);
    dep.j_real.resize(n);
    dep.j_imag.resize(n);

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, glass_names.joule_heat);
    syncVariableToHost<Real>(particles, glass_names.phi_real);
    syncVariableToHost<Real>(particles, glass_names.phi_imag);
    syncVariableToHost<Vecd>(particles, glass_names.e_real);
    syncVariableToHost<Vecd>(particles, glass_names.e_imag);
    syncVariableToHost<Vecd>(particles, glass_names.j_real);
    syncVariableToHost<Vecd>(particles, glass_names.j_imag);

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *joule = particles.getVariableDataByName<Real>(glass_names.joule_heat);
    const Real *phi_r = particles.getVariableDataByName<Real>(glass_names.phi_real);
    const Real *phi_i = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    const Vecd *e_r = particles.getVariableDataByName<Vecd>(glass_names.e_real);
    const Vecd *e_i = particles.getVariableDataByName<Vecd>(glass_names.e_imag);
    const Vecd *j_r = particles.getVariableDataByName<Vecd>(glass_names.j_real);
    const Vecd *j_i = particles.getVariableDataByName<Vecd>(glass_names.j_imag);
    for (size_t i = 0; i < n; ++i)
    {
        dep.position[i] = pos[i];
        dep.volumetric_measure[i] = vol[i];
        dep.joule_heat[i] = joule[i];
        dep.phi_real[i] = phi_r[i];
        dep.phi_imag[i] = phi_i[i];
        dep.e_real[i] = e_r[i];
        dep.e_imag[i] = e_i[i];
        dep.j_real[i] = j_r[i];
        dep.j_imag[i] = j_i[i];
    }
}

template <class ExecutionPolicy>
inline Rh200EmSolveResult runRh200OphelieEmSolveOnce(
    SolidBody &glass, Inner<> &glass_inner, const electromagnetics::ophelie::OphelieGlassFieldNames &glass_names,
    electromagnetics::ophelie::OphelieParameters &params,
    electromagnetics::ophelie::OphelieFrenchReducedCaseParams &french,
    Rh200EmParticleDepositionHost *raw_deposition_out = nullptr)
{
    using namespace electromagnetics::ophelie;
    BaseParticles &particles = glass.getBaseParticles();
    Rh200EmSolveResult result;
    result.glass_bbox = measureGlassParticleBbox(particles);
    result.base_current_per_loop = french.coil.current_per_loop;

    StateDynamics<ExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();

    result.em = runFrenchReducedEmPipeline<ExecutionPolicy>(glass, glass_inner, glass_names, params, french);
    const size_t n = particles.TotalRealParticles();
    result.raw_fields = collectRh200GlassEmFieldAuditSnapshot(particles, glass_names, n);
    result.raw_fields.p_joule_particle = result.em.joule_power_raw;
    if (raw_deposition_out != nullptr)
    {
        hostExtractEmDepositionFields(particles, glass_names, *raw_deposition_out);
    }

    result.equivalent_current_scale =
        calibrateFrenchCoilCurrentToTargetPower(french, params, result.em.joule_power_raw);
    result.em_power_scale_factor = result.equivalent_current_scale * result.equivalent_current_scale;

    result.em = runFrenchReducedEmPipeline<ExecutionPolicy>(glass, glass_inner, glass_names, params, french);

    if (ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        syncOphelieJouleHeatPrimaryForThermalOneWay<ExecutionPolicy>(glass, glass_names, params);
    }

    result.scaled_fields = collectRh200GlassEmFieldAuditSnapshot(particles, glass_names, n);
    syncVariableToHost<Real>(particles, glass_names.joule_heat);
    const Real *joule = particles.getVariableDataByName<Real>(glass_names.joule_heat);
    result.joule_heat_min = joule[0];
    result.joule_heat_max = joule[0];
    for (size_t i = 1; i < n; ++i)
    {
        result.joule_heat_min = std::min(result.joule_heat_min, joule[i]);
        result.joule_heat_max = std::max(result.joule_heat_max, joule[i]);
    }

    result.joule_power_w = hostVolWeightedSum(particles, glass_names.joule_heat, n);
    result.scaled_fields.p_joule_particle = result.joule_power_w;
    result.phi_eq_res_vol = result.em.phi_eq_res_vol;

    std::cout << "[rh200] EM solve once: n_glass=" << n << " P_raw=" << result.raw_fields.p_joule_particle
              << " W em_power_scale=" << result.em_power_scale_factor << " P_scaled=" << result.joule_power_w << " W"
              << " phi_eq_res_vol=" << result.phi_eq_res_vol << " joule_heat=[" << result.joule_heat_min << ","
              << result.joule_heat_max << "] W/m^3" << std::endl;
    return result;
}

inline void hostExtractJouleHeatField(BaseParticles &particles, const std::string &joule_field, StdVec<Real> &joule_out)
{
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, joule_field);
    const size_t n = particles.TotalRealParticles();
    const Real *joule = particles.getVariableDataByName<Real>(joule_field);
    joule_out.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        joule_out[i] = joule[i];
    }
}

inline void hostInstallJouleHeatOnBody(SPHBody &body, const StdVec<Real> &joule_host, const std::string &joule_field = kJouleHeatField)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    if (joule_host.size() != n)
    {
        std::cerr << "[rh200] ERROR: JouleHeat copy size mismatch: host=" << joule_host.size() << " particles=" << n
                  << std::endl;
        std::exit(1);
    }
    particles.registerStateVariable<Real>(joule_field, Real(0));
    Real *joule = particles.getVariableDataByName<Real>(joule_field);
    for (size_t i = 0; i < n; ++i)
    {
        joule[i] = joule_host[i];
    }
#if SPHINXSYS_USE_SYCL
    electromagnetics::ophelie::syncVariableToDevice<Real>(particles, joule_field);
#endif
}

inline void hostInstallEmCurrentOnBody(SPHBody &body, const Rh200EmParticleDepositionHost &dep,
                                       const electromagnetics::ophelie::OphelieGlassFieldNames &glass_names)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    if (dep.j_real.size() != n || dep.j_imag.size() != n)
    {
        std::cerr << "[rh200] ERROR: J current copy size mismatch: j_real=" << dep.j_real.size()
                  << " j_imag=" << dep.j_imag.size() << " particles=" << n << std::endl;
        std::exit(1);
    }
    particles.registerStateVariable<Vecd>(glass_names.j_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(glass_names.j_imag, ZeroData<Vecd>::value);
    Vecd *j_r = particles.getVariableDataByName<Vecd>(glass_names.j_real);
    Vecd *j_i = particles.getVariableDataByName<Vecd>(glass_names.j_imag);
    for (size_t i = 0; i < n; ++i)
    {
        j_r[i] = dep.j_real[i];
        j_i[i] = dep.j_imag[i];
    }
#if SPHINXSYS_USE_SYCL
    electromagnetics::ophelie::syncVariableToDevice<Vecd>(particles, glass_names.j_real);
    electromagnetics::ophelie::syncVariableToDevice<Vecd>(particles, glass_names.j_imag);
#endif
}

} // namespace rh200
} // namespace SPH

#endif // RH200_OPHELIE_EM_H
