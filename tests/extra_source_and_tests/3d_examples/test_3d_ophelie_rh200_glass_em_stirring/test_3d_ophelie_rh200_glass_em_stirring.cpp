/**
 * @file test_3d_ophelie_rh200_glass_em_stirring.cpp
 * @brief Step 0: RH200 single-phase glass stirring + thermal (glass-rotor only; wall adiabatic).
 *
 * Phases (single executable):
 *   --relax=1                     particle relaxation -> ./reload
 *   --reload=1 --run=1            flow + thermal stirring (no EM yet)
 *
 * STL units: meters (--geometry-scale=1.0 default).
 */
#include "io_environment.h"
#include "sphinxsys.h"
#include "rh200_fake_joule_heat.h"
#include "rh200_excitation_joule_audit.h"
#include "rh200_joule_heat_grid.h"
#include "rh200_material_preset.h"
#include "rh200_ophelie_em.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::rh200;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

std::string g_path_rotor = "./input/rotor.stl";
std::string g_path_full_glass = "./input/full_oil.stl";
std::string g_path_glass = "./input/oil.stl";

constexpr Real kRcd = 0.178 * 2.0;
constexpr Real kRotationSpeed = 10.472; /**< rad/s, ~100 RPM (2*pi*100/60). */
constexpr Real kRotorRevolutionPeriod = Real(2) * Pi / kRotationSpeed;
/** Default phase slip when sampling ~once per revolution (deg advance per output frame). */
constexpr Real kDefaultStateRecordPhaseSlipDeg = Real(24);

/** ~1 frame/rev plus phase slip: interval = T*(1 + slip/360); rotor advances slip deg per frame. */
inline Real rh200StaggeredOncePerRevRecordInterval(Real phase_slip_deg)
{
    return kRotorRevolutionPeriod * (Real(1) + phase_slip_deg / Real(360));
}

inline Real rh200DefaultStateRecordInterval()
{
    return rh200StaggeredOncePerRevRecordInterval(kDefaultStateRecordPhaseSlipDeg);
}

inline Real rh200StateRecordIntervalDegreesPerFrame(Real interval_s)
{
    return kRotationSpeed * interval_s * Real(180) / Pi;
}

inline size_t rh200EstimatedStateRecordFrames(Real end_time, Real interval_s)
{
    if (!(interval_s > TinyReal) || !(end_time >= Real(0)))
    {
        return 1;
    }
    return static_cast<size_t>(std::floor(end_time / interval_s)) + 1;
}

struct Rh200GlassStirringCli
{
    bool relax = false;
    bool run_flow = false;
    bool use_reload = false;
    Real dp = 0.008;
    Real geometry_scale = 1.0;
    Real end_time = 1.0;
    size_t relax_steps = 5000;
    std::string reload_dir = "./reload";
    bool state_recording = false;
    bool em_solve = false;
    FakeJouleMode joule_mode = FakeJouleMode::Off;
    Real joule_grid_spacing_factor = 1.0;
    Real joule_grid_bbox_margin = 2.0;
    Real state_record_interval = rh200DefaultStateRecordInterval();
    Real state_record_phase_slip_deg = kDefaultStateRecordPhaseSlipDeg;
    size_t state_record_frames_per_rev = 0;
    bool state_record_interval_user_set = false;
    bool state_record_frames_per_rev_user_set = false;
    Rh200MaterialPreset preset = Rh200MaterialPreset::DemoCurrent;
    Rh200FlowMaterialParams material;
    Rh200EmCli em;
};

inline void parseRh200GlassStirringCli(int ac, char *av[], Rh200GlassStirringCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--relax=1") == 0 || std::strcmp(av[i], "--relax") == 0 ||
            std::strcmp(av[i], "--relax=true") == 0)
        {
            cli.relax = true;
        }
        else if (std::strcmp(av[i], "--run=1") == 0 || std::strcmp(av[i], "--run") == 0)
        {
            cli.run_flow = true;
        }
        else if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0 ||
                 std::strcmp(av[i], "--reload=true") == 0)
        {
            cli.use_reload = true;
        }
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--geometry-scale=", 17) == 0)
        {
            cli.geometry_scale = static_cast<Real>(std::atof(av[i] + 17));
        }
        else if (std::strncmp(av[i], "--end-time=", 11) == 0)
        {
            cli.end_time = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--relax-steps=", 14) == 0)
        {
            cli.relax_steps = static_cast<size_t>(std::atoi(av[i] + 14));
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
        }
        else if (std::strcmp(av[i], "--state-recording=1") == 0 || std::strcmp(av[i], "--state-recording") == 0)
        {
            cli.state_recording = true;
        }
        else if (std::strncmp(av[i], "--joule-mode=", 13) == 0)
        {
            cli.joule_mode = parseFakeJouleMode(std::string(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--joule-grid-spacing-factor=", 28) == 0)
        {
            cli.joule_grid_spacing_factor = static_cast<Real>(std::atof(av[i] + 28));
        }
        else if (std::strncmp(av[i], "--joule-grid-bbox-margin=", 25) == 0)
        {
            cli.joule_grid_bbox_margin = static_cast<Real>(std::atof(av[i] + 25));
        }
        else if (std::strncmp(av[i], "--state-record-interval=", 24) == 0)
        {
            cli.state_record_interval = static_cast<Real>(std::atof(av[i] + 24));
            cli.state_record_interval_user_set = true;
        }
        else if (std::strncmp(av[i], "--state-record-frames-per-rev=", 30) == 0)
        {
            cli.state_record_frames_per_rev = static_cast<size_t>(std::atoi(av[i] + 30));
            cli.state_record_frames_per_rev_user_set = true;
        }
        else if (std::strncmp(av[i], "--state-record-phase-slip-deg=", 32) == 0)
        {
            cli.state_record_phase_slip_deg = static_cast<Real>(std::atof(av[i] + 32));
        }
        else if (std::strcmp(av[i], "--em-solve=1") == 0 || std::strcmp(av[i], "--em-solve") == 0)
        {
            cli.em_solve = true;
        }
        else if (std::strncmp(av[i], "--target-power=", 15) == 0)
        {
            cli.em.target_power = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strncmp(av[i], "--sigma0=", 9) == 0)
        {
            cli.em.sigma0 = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--frequency=", 12) == 0)
        {
            cli.em.frequency_hz = static_cast<Real>(std::atof(av[i] + 12));
        }
        else if (std::strncmp(av[i], "--coil-radius-factor=", 21) == 0)
        {
            cli.em.coil_radius_factor = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strncmp(av[i], "--coil-turns=", 13) == 0)
        {
            cli.em.coil_turns = static_cast<size_t>(std::atoi(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--coil-z-margin-low=", 20) == 0)
        {
            cli.em.coil_z_margin_low = static_cast<Real>(std::atof(av[i] + 20));
        }
        else if (std::strncmp(av[i], "--coil-z-margin-high=", 21) == 0)
        {
            cli.em.coil_z_margin_high = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strncmp(av[i], "--preset=", 9) == 0)
        {
            cli.preset = parseRh200MaterialPreset(std::string(av[i] + 9));
        }
    }
    cli.material = makeRh200FlowMaterialParams(cli.preset);
    applyRh200MaterialPresetToEmDefaults(cli.preset, cli.em.sigma0, cli.em.frequency_hz, cli.em.target_power);

    if (!cli.state_record_interval_user_set)
    {
        if (cli.state_record_frames_per_rev_user_set && cli.state_record_frames_per_rev > 0)
        {
            cli.state_record_interval =
                kRotorRevolutionPeriod / static_cast<Real>(cli.state_record_frames_per_rev);
        }
        else
        {
            cli.state_record_interval = rh200StaggeredOncePerRevRecordInterval(cli.state_record_phase_slip_deg);
        }
    }
}

inline bool isRh200CustomCommandLineOption(const char *arg)
{
    return std::strncmp(arg, "--dp=", 5) == 0 || std::strncmp(arg, "--geometry-scale=", 17) == 0 ||
           std::strncmp(arg, "--end-time=", 11) == 0 || std::strncmp(arg, "--relax-steps=", 14) == 0 ||
           std::strncmp(arg, "--reload-dir=", 13) == 0 || std::strcmp(arg, "--relax=1") == 0 ||
           std::strcmp(arg, "--relax") == 0 || std::strcmp(arg, "--relax=true") == 0 ||
           std::strcmp(arg, "--run=1") == 0 || std::strcmp(arg, "--run") == 0 || std::strcmp(arg, "--reload=1") == 0 ||
           std::strcmp(arg, "--reload") == 0 || std::strcmp(arg, "--reload=true") == 0 ||
           std::strcmp(arg, "--state-recording=1") == 0 || std::strcmp(arg, "--state-recording") == 0 ||
           std::strncmp(arg, "--joule-mode=", 13) == 0 ||
           std::strncmp(arg, "--joule-grid-spacing-factor=", 28) == 0 ||
           std::strncmp(arg, "--joule-grid-bbox-margin=", 25) == 0 ||
           std::strncmp(arg, "--state-record-interval=", 24) == 0 ||
           std::strncmp(arg, "--state-record-frames-per-rev=", 30) == 0 ||
           std::strncmp(arg, "--state-record-phase-slip-deg=", 32) == 0 ||
           std::strcmp(arg, "--em-solve=1") == 0 || std::strcmp(arg, "--em-solve") == 0 ||
           std::strncmp(arg, "--target-power=", 15) == 0 || std::strncmp(arg, "--sigma0=", 9) == 0 ||
           std::strncmp(arg, "--frequency=", 12) == 0 || std::strncmp(arg, "--coil-radius-factor=", 21) == 0 ||
           std::strncmp(arg, "--coil-turns=", 13) == 0 || std::strncmp(arg, "--coil-z-margin-low=", 20) == 0 ||
           std::strncmp(arg, "--coil-z-margin-high=", 21) == 0 || std::strncmp(arg, "--preset=", 9) == 0;
}

inline StdVec<std::string> filterRh200SphCommandLine(int ac, char *av[])
{
    StdVec<std::string> filtered;
    filtered.emplace_back(av[0]);
    for (int i = 1; i < ac; ++i)
    {
        if (!isRh200CustomCommandLineOption(av[i]))
        {
            filtered.emplace_back(av[i]);
        }
    }
    return filtered;
}

inline void callRh200SphCommandLineOptions(SPHSystem &sph_system, int ac, char *av[])
{
    const StdVec<std::string> filtered_arguments = filterRh200SphCommandLine(ac, av);
    StdVec<char *> filtered_argv;
    filtered_argv.reserve(filtered_arguments.size());
    for (auto &argument : filtered_arguments)
    {
        filtered_argv.push_back(const_cast<char *>(argument.c_str()));
    }
    sph_system.handleCommandlineOptions(static_cast<int>(filtered_argv.size()), filtered_argv.data());
}

inline void logStlBoundingBox(const std::string &label, const std::string &stl_path, Real geometry_scale)
{
    TriangleMeshShapeSTL shape(stl_path, Vec3d::Zero(), geometry_scale, label);
    const BoundingBoxd bounds = shape.findBounds();
    const Vecd bmin(bounds.lower_[0], bounds.lower_[1], bounds.lower_[2]);
    const Vecd bmax(bounds.upper_[0], bounds.upper_[1], bounds.upper_[2]);
    const Vecd size = bmax - bmin;
    const Vecd center = Real(0.5) * (bmin + bmax);
    std::cout << "[rh200] STL bbox " << label << " scale=" << geometry_scale << "\n"
              << "  min=(" << bmin.transpose() << ") max=(" << bmax.transpose() << ")\n"
              << "  size=(" << size.transpose() << ") center=(" << center.transpose() << ")\n";
}

inline void writeRh200GeometryAuditCsv(const std::string &path, Real geometry_scale)
{
    std::ofstream out(path);
    if (!out)
    {
        std::cerr << "[rh200] warning: could not write " << path << std::endl;
        return;
    }
    out << "label,geometry_scale,xmin,ymin,zmin,xmax,ymax,zmax,Lx,Ly,Lz,cx,cy,cz\n";
    const StdVec<std::pair<std::string, std::string>> meshes = {
        {"glass_fluid", g_path_glass},
        {"wall_container", g_path_full_glass},
        {"rotor", g_path_rotor},
    };
    for (const auto &entry : meshes)
    {
        TriangleMeshShapeSTL shape(entry.second, Vec3d::Zero(), geometry_scale, entry.first);
        const BoundingBoxd bounds = shape.findBounds();
        const Vecd bmin(bounds.lower_[0], bounds.lower_[1], bounds.lower_[2]);
        const Vecd bmax(bounds.upper_[0], bounds.upper_[1], bounds.upper_[2]);
        const Vecd size = bmax - bmin;
        const Vecd center = Real(0.5) * (bmin + bmax);
        out << entry.first << "," << geometry_scale << "," << bmin[0] << "," << bmin[1] << "," << bmin[2] << ","
            << bmax[0] << "," << bmax[1] << "," << bmax[2] << "," << size[0] << "," << size[1] << "," << size[2]
            << "," << center[0] << "," << center[1] << "," << center[2] << "\n";
    }
    std::cout << "[rh200] wrote geometry audit: " << path << std::endl;
}

class Rotor : public ComplexShape
{
  public:
    explicit Rotor(const std::string &shape_name, Real geometry_scale)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(g_path_rotor, Vec3d::Zero(), geometry_scale);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name, Real geometry_scale, Real resolution_ref)
        : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(4.0 * resolution_ref, g_path_full_glass, Vec3d::Zero(), geometry_scale);
        subtract<TriangleMeshShapeSTL>(g_path_full_glass, Vec3d::Zero(), geometry_scale);
    }
};

class GlassGeometry : public ComplexShape
{
  public:
    explicit GlassGeometry(const std::string &shape_name, Real geometry_scale)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(g_path_glass, Vec3d::Zero(), geometry_scale);
    }
};

class BoundingBoxCalculator : public TriangleMeshShapeSTL
{
  public:
    BoundingBoxCalculator(const std::string &stl_path, const Vec3d &translation, Real scale, Real resolution_ref)
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

class MergedBoundingBox
{
  public:
    MergedBoundingBox(const StdVec<std::string> &files, const Vec3d &translation, Real scale, Real resolution_ref)
        : files_(files), translation_(translation), scale_(scale), resolution_ref_(resolution_ref)
    {
    }

    BoundingBoxd calculate()
    {
        Vec3d global_min = Vec3d::Constant(std::numeric_limits<Real>::max());
        Vec3d global_max = Vec3d::Constant(-std::numeric_limits<Real>::max());
        for (const std::string &fpath : files_)
        {
            BoundingBoxCalculator calc(fpath, translation_, scale_, resolution_ref_);
            const BoundingBoxd b = calc.calculate();
            for (int d = 0; d < 3; ++d)
            {
                global_min[d] = std::min(global_min[d], b.lower_[d]);
                global_max[d] = std::max(global_max[d], b.upper_[d]);
            }
        }
        return BoundingBoxd(global_min, global_max);
    }

  private:
    StdVec<std::string> files_;
    Vec3d translation_;
    Real scale_;
    Real resolution_ref_;
};

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

inline std::string rh200AnimationFrameTag(size_t frame_index)
{
    std::ostringstream oss;
    oss << "ite_" << std::setw(10) << std::setfill('0') << frame_index;
    return oss.str();
}

struct Rh200SynchronizedAnimationWriter
{
    size_t frame_index = 0;
    std::string glass_name_;
    std::string rotor_proxy_name_;
    std::string rotor_sph_name_;
    bool write_rotor_sph_ = false;
    std::ofstream manifest_;

    void open(const std::string &manifest_path, const std::string &glass_name, const std::string &rotor_proxy_name,
              const std::string &rotor_sph_name = "")
    {
        glass_name_ = glass_name;
        rotor_proxy_name_ = rotor_proxy_name;
        rotor_sph_name_ = rotor_sph_name;
        write_rotor_sph_ = !rotor_sph_name_.empty();
        manifest_.open(manifest_path);
        manifest_ << "frame_index,frame_tag,physical_time_s,glass_vtp,rotor_proxy_vtp";
        if (write_rotor_sph_)
        {
            manifest_ << ",rotor_sph_vtp";
        }
        manifest_ << "\n";
    }

    template <typename GlassRecorder, typename UpdateRotorDynamics, typename RotorWriter, typename RotorSphRecorder>
    void write(Real physical_time, GlassRecorder &glass_recorder, UpdateRotorDynamics &update_rotor,
               RotorWriter &rotor_writer, RotorSphRecorder *rotor_sph_recorder = nullptr)
    {
        update_rotor.exec(physical_time);
        glass_recorder.writeToFile(frame_index);
        rotor_writer.writeToFile(frame_index);
        if (write_rotor_sph_ && rotor_sph_recorder != nullptr)
        {
            rotor_sph_recorder->writeToFile(frame_index);
        }
        const std::string tag = rh200AnimationFrameTag(frame_index);
        manifest_ << frame_index << "," << tag << "," << std::fixed << std::setprecision(9) << physical_time << ","
                  << glass_name_ << "_" << tag << ".vtp," << rotor_proxy_name_ << "_" << tag << ".vtp";
        if (write_rotor_sph_)
        {
            manifest_ << "," << rotor_sph_name_ << "_" << tag << ".vtp";
        }
        manifest_ << "\n";
        ++frame_index;
    }

    void close()
    {
        if (manifest_.is_open())
        {
            manifest_.close();
        }
    }
};

inline Real measureStlGlassVolume(const std::string &stl_path, Real geometry_scale)
{
    const Rh200GlassBboxAudit stl_bbox = measureStlGlassBbox(stl_path, geometry_scale);
    return stl_bbox.size[0] * stl_bbox.size[1] * stl_bbox.size[2];
}

inline void writeRh200ExcitationToJouleAuditFromEmSolve(
    const Rh200GlassStirringCli &cli, const Rh200EmSolveResult &em_result,
    const electromagnetics::ophelie::OphelieFrenchReducedCaseParams &french, BaseParticles &glass_particles,
    Real rho_glass, Real cp_glass, const Rh200EmParticleDepositionHost *raw_dep,
    const Rh200EmParticleDepositionHost *scaled_dep, const Rh200JouleHeatGridSpec *grid_spec, Real grid_rescale_factor,
    Real p_scaled_grid_sample)
{
    Rh200ExcitationToJouleAuditRecord rec;
    rec.case_name = std::string("rh200_glass_em_stirring|") + rh200MaterialPresetName(cli.preset);
    rec.joule_mode = fakeJouleModeName(cli.joule_mode);
    rec.em_mode = "edge-flux-complex";
    rec.frequency_hz = cli.em.frequency_hz;
    rec.omega = Real(2) * Pi * cli.em.frequency_hz;
    rec.coil_radius = french.coil.loop_radius;
    rec.coil_z_min = french.coil.z_min;
    rec.coil_z_max = french.coil.z_max;
    rec.coil_turns = french.coil.num_loops;
    rec.coil_center_x = french.coil.stack_center[0];
    rec.coil_center_y = french.coil.stack_center[1];
    rec.coil_center_z = french.coil.stack_center[2];
    rec.base_current_per_loop = em_result.base_current_per_loop;
    rec.equivalent_current_scale = em_result.equivalent_current_scale;
    rec.em_power_scale_factor = em_result.em_power_scale_factor;
    rec.grid_rescale_factor = grid_rescale_factor;
    rec.target_power = cli.em.target_power;
    rec.raw_fields = em_result.raw_fields;
    rec.scaled_fields = em_result.scaled_fields;
    rec.thermal_mass = measureRh200GlassParticleThermalMass(glass_particles, rho_glass, cp_glass);
    rec.stl_glass_volume = measureStlGlassVolume(g_path_glass, cli.geometry_scale);
    if (raw_dep != nullptr && grid_spec != nullptr)
    {
        rec.p_joule_raw_grid_sample_initial = hostSampleRawGridPowerFromDeposition(*grid_spec, *raw_dep);
    }
    rec.p_joule_scaled_grid_sample_initial = p_scaled_grid_sample;
    if (rec.p_joule_scaled_grid_sample_initial <= TinyReal && scaled_dep != nullptr && grid_spec != nullptr)
    {
        Real unused = Real(1);
        rec.p_joule_scaled_grid_sample_initial =
            hostSampleGridPowerFromDeposition(*grid_spec, *scaled_dep, true, grid_rescale_factor, unused);
    }
    writeRh200ExcitationToJouleAuditCsv(kExcitationToJouleAuditCsv, rec);
}

inline Rh200EmSolveResult prepareRh200EmJouleHeat(const Rh200GlassStirringCli &cli, int ac, char *av[],
                                                  StdVec<Real> &joule_heat_out,
                                                  Rh200EmParticleDepositionHost *deposition_out = nullptr,
                                                  Rh200EmParticleDepositionHost *raw_deposition_out = nullptr,
                                                  electromagnetics::ophelie::OphelieFrenchReducedCaseParams *french_out =
                                                      nullptr)
{
    const StdVec<std::string> stl_files = {g_path_glass};
    MergedBoundingBox merger(stl_files, Vec3d::Zero(), cli.geometry_scale, cli.dp);
    const BoundingBoxd system_domain_bounds = merger.calculate();

    SPHSystem em_system(system_domain_bounds, cli.dp);
    em_system.setReloadParticles(true);
    callRh200SphCommandLineOptions(em_system, ac, av);
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    const Real rho0_solid = cli.material.rho_glass;
    SolidBody glass_em(em_system, makeShared<GlassGeometry>("GlassGeometry", cli.geometry_scale));
    glass_em.defineAdaptationRatios(1.3, 1.0);
    glass_em.defineMatterMaterial<Solid>(rho0_solid);
    glass_em.generateParticles<BaseParticles, Reload>(glass_em.Name());

    em_system.initializeSystemCellLinkedLists();
    em_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_em, glass_names);
    (void)register_glass;

    Inner<> glass_inner_em(glass_em);

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    params.target_joule_power_ = cli.em.target_power;
    params.sigma_glass_ = cli.em.sigma0;
    params.frequency_ = cli.em.frequency_hz;

    OphelieTestCliOptions ophelie_cli;
    configureRh200OphelieEmParameters(params, ophelie_cli);

    const Rh200GlassBboxAudit particle_bbox = measureGlassParticleBbox(glass_em.getBaseParticles());
    fillFrenchParamsFromRh200Bbox(particle_bbox, cli.dp, cli.em, french, params);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    const Rh200GlassBboxAudit stl_bbox = measureStlGlassBbox(g_path_glass, cli.geometry_scale);
    writeRh200CoilGeometryLogCsv(kRh200CoilGeometryCsv, cli.geometry_scale, stl_bbox, particle_bbox, french);

    std::cout << "[rh200] OPHELIE solve on isolated glass reload (dp=" << cli.dp
              << " n_glass=" << particle_bbox.particle_count << ")" << std::endl;

    const Rh200EmSolveResult em_result = runRh200OphelieEmSolveOnce<MainExecutionPolicy>(
        glass_em, glass_inner_em, glass_names, params, french, raw_deposition_out);
    hostExtractJouleHeatField(glass_em.getBaseParticles(), glass_names.joule_heat, joule_heat_out);
    if (deposition_out != nullptr)
    {
        hostExtractEmDepositionFields(glass_em.getBaseParticles(), glass_names, *deposition_out);
    }
    if (french_out != nullptr)
    {
        *french_out = french;
    }
    writeRh200EmSolveSummaryCsv(kRh200EmSolveCsv, em_result, french);
    return em_result;
}

inline Rh200JouleHeatGridBuildReport buildRh200EulerianJouleHeatGridFromEm(
    const Rh200GlassStirringCli &cli, Rh200JouleHeatEmGridBundle &grid_bundle, const Rh200EmParticleDepositionHost &em_dep,
    Rh200JouleHeatGridSpec *grid_spec_out = nullptr)
{
    const Real grid_spacing = cli.joule_grid_spacing_factor * cli.dp;
    const Rh200JouleHeatGridSpec spec = makeRh200JouleHeatGridSpecFromStl(g_path_full_glass, cli.geometry_scale,
                                                                          grid_spacing, cli.joule_grid_bbox_margin);
    if (grid_spec_out != nullptr)
    {
        *grid_spec_out = spec;
    }
    logStlBoundingBox("full_oil_joule_grid", g_path_full_glass, cli.geometry_scale);
    std::cout << "[rh200] em-grid: domain=full_oil.stl spacing=" << grid_spacing << " (" << cli.joule_grid_spacing_factor
              << "*dp) bbox_margin=" << cli.joule_grid_bbox_margin << "*spacing"
              << " grid_bbox_min=(" << spec.lower_bound_.transpose() << ") max=(" << spec.upper_bound_.transpose()
              << ") nodes=" << spec.nx_ << "x" << spec.ny_ << "x" << spec.nz_ << std::endl;
    const Rh200JouleHeatGridBuildReport report =
        buildRh200JouleHeatEmGridBundle(grid_bundle, spec, em_dep, cli.em.target_power);
    writeRh200JouleHeatGridSummaryCsv(kJouleHeatGridSummaryCsv, spec, report);
    writeRh200JouleHeatEmGridVti(kJouleHeatGridVti, grid_bundle);
    std::cout << "[rh200] em-grid build: P_em_particle=" << report.p_em_particle
              << " P_grid_sample_after_scale=" << report.p_grid_sample_after_scale
              << " rel_l2=" << report.sample_vs_em_rel_l2 << " rel_max=" << report.sample_vs_em_rel_max << std::endl;
    return report;
}

int runParticleRelaxation(const Rh200GlassStirringCli &cli, int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    std::cerr << "[rh200] ERROR: particle relaxation requires SYCL/GPU build (SPHINXSYS_USE_SYCL)." << std::endl;
    return 1;
#else
    const StdVec<std::string> stl_files = {g_path_rotor, g_path_full_glass};
    MergedBoundingBox merger(stl_files, Vec3d::Zero(), cli.geometry_scale, cli.dp);
    const BoundingBoxd system_domain_bounds = merger.calculate();

    SPHSystem sph_system(system_domain_bounds, cli.dp);
    sph_system.setRunParticleRelaxation(true);
    callRh200SphCommandLineOptions(sph_system, ac, av);

    FluidBody glass(sph_system, makeShared<GlassGeometry>("GlassGeometry", cli.geometry_scale));
    glass.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary", cli.geometry_scale, cli.dp));
    LevelSetShape &wall_level_set_shape =
        wall_boundary.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody rotor(sph_system, makeShared<Rotor>("Rotor", cli.geometry_scale));
    LevelSetShape &rotor_level_set_shape = rotor.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    rotor.generateParticles<BaseParticles, Lattice>();

    NearShapeSurface wall_near_surface(wall_boundary);
    NearShapeSurface rotor_near_surface(rotor);

    Inner<> wall_inner(wall_boundary);
    Inner<> rotor_inner(rotor);
    Inner<> glass_inner(glass);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    StdVec<RealBody *> relax_bodies = {&wall_boundary, &rotor};
    StdVec<RealBody *> all_bodies = {&wall_boundary, &rotor, &glass};

    host_methods.addStateDynamics<RandomizeParticlePositionCK>(wall_boundary).exec();
    host_methods.addStateDynamics<RandomizeParticlePositionCK>(rotor).exec();

    ParticleDynamicsGroup update_cell_linked_list = main_methods.addCellLinkedListDynamics(all_bodies);
    ParticleDynamicsGroup update_relation;
    update_relation.add(&main_methods.addRelationDynamics(wall_inner));
    update_relation.add(&main_methods.addRelationDynamics(rotor_inner));
    update_relation.add(&main_methods.addRelationDynamics(glass_inner));
    ParticleDynamicsGroup update_configuration = update_cell_linked_list + update_relation;

    ParticleDynamicsGroup relaxation_residual;
    relaxation_residual.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(wall_inner)
             .addPostStateDynamics<LevelsetKernelGradientIntegral>(wall_boundary, wall_level_set_shape));
    relaxation_residual.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(rotor_inner)
             .addPostStateDynamics<LevelsetKernelGradientIntegral>(rotor, rotor_level_set_shape));

    ReduceDynamicsGroup relaxation_scaling =
        main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(relax_bodies);

    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(relax_bodies);
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(wall_near_surface));
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(rotor_near_surface));

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(sph_system);

    body_state_recorder.writeToFile(0);

    const size_t vtp_interval = cli.relax_steps >= 1000 ? 1000 : std::max<size_t>(100, cli.relax_steps / 5);
    std::cout << "[rh200] SYCL-CK relax: dp=" << cli.dp << " geometry_scale=" << cli.geometry_scale
              << " steps=" << cli.relax_steps << " n_glass=" << glass.getBaseParticles().TotalRealParticles()
              << " n_wall=" << wall_boundary.getBaseParticles().TotalRealParticles()
              << " n_rotor=" << rotor.getBaseParticles().TotalRealParticles() << std::endl;

    for (size_t ite = 0; ite < cli.relax_steps; ++ite)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        const Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);

        if (ite % vtp_interval == 0 && ite != 0)
        {
            std::cout << "[rh200] relaxation step " << ite << std::endl;
            body_state_recorder.writeToFile(static_cast<int>(ite));
        }
    }

    body_state_recorder.writeToFile(static_cast<int>(cli.relax_steps));
    write_particle_reload_files.writeToFile(0);
    std::cout << "[rh200] SYCL-CK relaxation done; reload written (GlassGeometry, WallBoundary, Rotor)" << std::endl;
    return 0;
#endif
}

int runGlassStirringFlow(const Rh200GlassStirringCli &cli, int ac, char *av[])
{
    if (cli.joule_mode != FakeJouleMode::Off && cli.joule_mode != FakeJouleMode::EmFixed &&
        cli.joule_mode != FakeJouleMode::EmGrid)
    {
        std::cerr << "[rh200] ERROR: joule-mode supports only off|em-fixed|em-grid in current production flow."
                  << std::endl;
        return 1;
    }

    StdVec<Real> em_fixed_joule_host;
    Rh200EmParticleDepositionHost em_dep;
    Rh200EmParticleDepositionHost em_dep_raw;
    Rh200EmSolveResult em_prepare_result;
    OphelieFrenchReducedCaseParams em_french_for_audit;
    Rh200JouleHeatGridSpec joule_grid_spec;
    Rh200JouleHeatGridBuildReport joule_grid_report;
    const Rh200JouleHeatGridSpec *joule_grid_spec_ptr = nullptr;
    Real joule_grid_rescale_factor = Real(1);
    Real joule_grid_scaled_sample_power = Real(0);
    Rh200JouleHeatEmGridBundle em_grid_bundle;
#if SPHINXSYS_USE_SYCL
    Rh200JouleHeatGridDevice joule_grid_device;
    Rh200JouleHeatGridDevice j_real_x_grid_device;
    Rh200JouleHeatGridDevice j_real_y_grid_device;
    Rh200JouleHeatGridDevice j_real_z_grid_device;
    Rh200JouleHeatGridDevice j_imag_x_grid_device;
    Rh200JouleHeatGridDevice j_imag_y_grid_device;
    Rh200JouleHeatGridDevice j_imag_z_grid_device;
#endif
    if (cli.joule_mode == FakeJouleMode::EmFixed)
    {
        em_prepare_result = prepareRh200EmJouleHeat(cli, ac, av, em_fixed_joule_host, &em_dep, &em_dep_raw,
                                                      &em_french_for_audit);
    }
    else if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
        em_prepare_result = prepareRh200EmJouleHeat(cli, ac, av, em_fixed_joule_host, &em_dep, &em_dep_raw,
                                                      &em_french_for_audit);
        joule_grid_report = buildRh200EulerianJouleHeatGridFromEm(cli, em_grid_bundle, em_dep, &joule_grid_spec);
        joule_grid_spec_ptr = &joule_grid_spec;
        joule_grid_rescale_factor = joule_grid_report.grid_scale_factor;
        joule_grid_scaled_sample_power = joule_grid_report.p_grid_sample_after_scale;
#if SPHINXSYS_USE_SYCL
        joule_grid_device.upload(em_grid_bundle.joule_heat);
        j_real_x_grid_device.upload(em_grid_bundle.j_real_x);
        j_real_y_grid_device.upload(em_grid_bundle.j_real_y);
        j_real_z_grid_device.upload(em_grid_bundle.j_real_z);
        j_imag_x_grid_device.upload(em_grid_bundle.j_imag_x);
        j_imag_y_grid_device.upload(em_grid_bundle.j_imag_y);
        j_imag_z_grid_device.upload(em_grid_bundle.j_imag_z);
#else
        std::cerr << "[rh200] ERROR: em-grid requires SYCL/GPU build." << std::endl;
        return 1;
#endif
    }

    const StdVec<std::string> stl_files = {g_path_rotor, g_path_full_glass};
    MergedBoundingBox merger(stl_files, Vec3d::Zero(), cli.geometry_scale, cli.dp);
    const BoundingBoxd system_domain_bounds = merger.calculate();

    logStlBoundingBox("glass_fluid", g_path_glass, cli.geometry_scale);
    writeRh200GeometryAuditCsv("./output/rh200_geometry_audit.csv", cli.geometry_scale);

    SPHSystem sph_system(system_domain_bounds, cli.dp);
    sph_system.setReloadParticles(true);
    callRh200SphCommandLineOptions(sph_system, ac, av);
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    const Real rho0_glass = cli.material.rho_glass;
    const Real gravity_g = 9.81;
    const Real U_glass = 1.0 * kRotationSpeed * 0.5 * kRcd;
    const Real c_glass = 10.0 * U_glass;
    const Real mu_glass = cli.material.mu_glass;
    const Real rho0_solid = cli.material.rho_rotor;

    const Real cp_glass = cli.material.cp_glass;
    const Real cv_glass = cp_glass * rho0_glass;
    const Real k_glass = cli.material.k_glass;
    const Real k_rotor_solid = cli.material.k_rotor_solid;
    const Real k_glass_rotor_contact = 2.0 * k_glass * k_rotor_solid / (k_glass + k_rotor_solid);
    const Real initial_temperature_glass = cli.material.t_initial_glass;
    const Real cp_rotor = cli.material.cp_rotor;
    const Real cv_rotor = cp_rotor * rho0_solid;

    logRh200MaterialPreset(cli.preset, cli.material, cli.em.sigma0, cli.em.frequency_hz, cli.em.target_power);

    SolidBody rotor(sph_system, makeShared<Rotor>("Rotor", cli.geometry_scale));
    rotor.defineAdaptationRatios(1.3, 1.0);
    rotor.defineMatterMaterial<Solid>(rho0_solid);
    rotor.generateParticles<BaseParticles, Reload>(rotor.Name());

    FluidBody glass(sph_system, makeShared<GlassGeometry>("GlassGeometry", cli.geometry_scale));
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_glass, c_glass);
    glass.addMaterialProperty<Viscosity>(mu_glass);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());

    if (cli.joule_mode == FakeJouleMode::EmFixed || cli.joule_mode == FakeJouleMode::EmGrid)
    {
        const Rh200EmParticleDepositionHost *scaled_dep_ptr = &em_dep;
        writeRh200ExcitationToJouleAuditFromEmSolve(
            cli, em_prepare_result, em_french_for_audit, glass.getBaseParticles(), rho0_glass, cp_glass, &em_dep_raw,
            scaled_dep_ptr, joule_grid_spec_ptr, joule_grid_rescale_factor, joule_grid_scaled_sample_power);
    }

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary", cli.geometry_scale, cli.dp));
    wall_boundary.defineAdaptationRatios(1.3, 1.0);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.Name());

    TriangleMeshShapeSTL rotor_shape(g_path_rotor, Vec3d::Zero(), cli.geometry_scale);
    ObserverBody rotor_proxy(sph_system, "RotorProxy");
    rotor_proxy.defineAdaptationRatios(2.0);
    rotor_proxy.generateParticles<ObserverParticles>(rotor_shape);

    Inner<> glass_inner(glass);
    Inner<> rotor_inner(rotor);
    Contact<> glass_wall_contact(glass, {&wall_boundary, &rotor});
    Contact<> glass_rotor_thermal_contact(glass, {&rotor});
    Contact<> rotor_fluid_contact(rotor, {&glass});

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary).exec();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(rotor).exec();

    auto &update_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(glass);
    auto &update_rotor_cell_linked_list = main_methods.addCellLinkedListDynamics(rotor);
    auto &update_wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &update_glass_relations = main_methods.addRelationDynamics(glass_inner, glass_wall_contact);
    auto &update_rotor_fluid_contact = main_methods.addRelationDynamics(rotor_fluid_contact);
    auto &fluid_particle_sorting = main_methods.addSortDynamics(glass);

    Gravity gravity(Vecd(0.0, 0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(glass, gravity);
    auto &fluid_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(glass);
    auto &fluid_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(glass);

    auto &acoustic_step_1st_half =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(glass_wall_contact);
    auto &acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK,
                                                    NoKernelCorrectionCK>(glass_inner);
    auto &acoustic_step_2nd_half_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_wall_contact);
    acoustic_step_2nd_half.addPostContactInteraction(acoustic_step_2nd_half_wall);
    auto &density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
            .addPostContactInteraction(glass_wall_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(glass);

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(
            glass_inner);
    auto &fluid_viscous_force_from_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(
            glass_wall_contact);
    fluid_viscous_force.addPostContactInteraction(fluid_viscous_force_from_wall);
    auto &viscous_force_from_fluid_on_rotor =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::ViscousForceFromFluid, std::remove_reference_t<decltype(fluid_viscous_force_from_wall)>>(
            rotor_fluid_contact);
    auto &pressure_force_from_fluid_on_rotor =
        main_methods.addInteractionDynamicsWithUpdate<
            FSI::PressureForceFromFluid, std::remove_reference_t<decltype(acoustic_step_2nd_half_wall)>>(
            rotor_fluid_contact);

    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(glass, U_glass);
    auto &acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(glass, 0.55);

    using DiffusiveInitialCondition =
        StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, ConstantValue<Real>>>;
    const std::string temperature_name = "Temperature";
    DiffusiveInitialCondition initialize_glass_temperature(glass, temperature_name, initial_temperature_glass);
    DiffusiveInitialCondition initialize_rotor_temperature(rotor, temperature_name, initial_temperature_glass);

    IsotropicDiffusion glass_heat_diffusion(temperature_name, k_glass, cv_glass);
    IsotropicDiffusion heat_diffusion_glass_rotor(temperature_name, k_glass_rotor_contact, cv_glass);
    IsotropicDiffusion rotor_heat_diffusion(temperature_name, k_rotor_solid, cv_rotor);
    IsotropicDiffusion heat_diffusion_rotor_glass(temperature_name, k_glass_rotor_contact, cv_rotor);

    using ThermalRelaxationComplex = RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                              Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>,
                              Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>>>>;

    ThermalRelaxationComplex glass_thermal_relaxation(DynamicsArgs(glass_inner, &glass_heat_diffusion),
                                                      DynamicsArgs(glass_rotor_thermal_contact,
                                                                   &heat_diffusion_glass_rotor));
    Contact<> rotor_glass_thermal_contact(rotor, {&glass});
    ThermalRelaxationComplex rotor_thermal_relaxation(DynamicsArgs(rotor_inner, &rotor_heat_diffusion),
                                                      DynamicsArgs(rotor_glass_thermal_contact,
                                                                   &heat_diffusion_rotor_glass));

    const bool enable_joule_heat = jouleHeatEnabled(cli.joule_mode);
    const Real joule_target_power = cli.em.target_power;
    const Rh200GlassThermalMassAudit glass_thermal_mass =
        enable_joule_heat ? measureRh200GlassParticleThermalMass(glass.getBaseParticles(), rho0_glass, cp_glass)
                          : Rh200GlassThermalMassAudit{};
    if (cli.joule_mode == FakeJouleMode::EmFixed)
    {
        hostInstallJouleHeatOnBody(glass, em_fixed_joule_host, kJouleHeatField);
        std::cout << "[rh200] em-fixed: installed frozen JouleHeat on fluid glass (" << em_fixed_joule_host.size()
                  << " particles)" << std::endl;
    }
    else if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
        glass.getBaseParticles().registerStateVariable<Real>(kJouleHeatField, Real(0));
        std::cout << "[rh200] em-grid: JouleHeat sampled each step from Eulerian grid (full_oil.stl domain)"
                  << std::endl;
    }

    OphelieGlassFieldNames flow_glass_em_names;
    const bool write_em_current_to_vtp =
        cli.joule_mode == FakeJouleMode::EmFixed || cli.joule_mode == FakeJouleMode::EmGrid;
    if (cli.joule_mode == FakeJouleMode::EmFixed)
    {
        hostInstallEmCurrentOnBody(glass, em_dep, flow_glass_em_names);
        std::cout << "[rh200] installed frozen JReal/JImag on fluid glass (" << em_dep.j_real.size()
                  << " particles) for VTP output" << std::endl;
    }
    else if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
        BaseParticles &glass_particles = glass.getBaseParticles();
        glass_particles.registerStateVariable<Vecd>(flow_glass_em_names.j_real, ZeroData<Vecd>::value);
        glass_particles.registerStateVariable<Vecd>(flow_glass_em_names.j_imag, ZeroData<Vecd>::value);
        std::cout << "[rh200] em-grid: JReal/JImag sampled each step from Eulerian current grids" << std::endl;
    }

    StateDynamics<MainExecutionPolicy, ApplyFakeJouleHeatToTemperatureCK> apply_fake_joule_heat(glass, cp_glass,
                                                                                                  rho0_glass);
#if SPHINXSYS_USE_SYCL
    std::unique_ptr<StateDynamics<MainExecutionPolicy, SampleJouleHeatFromGridCK>> sample_joule_heat_from_grid;
    std::unique_ptr<StateDynamics<MainExecutionPolicy, SampleCurrentFromGridCK>> sample_current_from_grid;
    std::unique_ptr<ReduceDynamicsCK<MainExecutionPolicy, OutOfGridJouleSampleReduceCK>> out_of_grid_joule_reduce;
    if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
        sample_joule_heat_from_grid =
            std::make_unique<StateDynamics<MainExecutionPolicy, SampleJouleHeatFromGridCK>>(glass, joule_grid_device);
        sample_current_from_grid = std::make_unique<StateDynamics<MainExecutionPolicy, SampleCurrentFromGridCK>>(
            glass, j_real_x_grid_device, j_real_y_grid_device, j_real_z_grid_device, j_imag_x_grid_device,
            j_imag_y_grid_device, j_imag_z_grid_device, flow_glass_em_names.j_real, flow_glass_em_names.j_imag);
        out_of_grid_joule_reduce = std::make_unique<ReduceDynamicsCK<MainExecutionPolicy, OutOfGridJouleSampleReduceCK>>(
            glass, joule_grid_device);
    }
#endif
    ReduceDynamicsCK<MainExecutionPolicy, FakeJoulePowerReduceCK> fake_joule_power_reduce(glass);
    ReduceDynamicsCK<MainExecutionPolicy, ThermalEnergyAboveInitialReduceCK> thermal_energy_reduce(
        glass, initial_temperature_glass, cp_glass, rho0_glass);
    ReduceDynamicsCK<MainExecutionPolicy, TemperatureMaxReduceCK> temperature_max_reduce(glass);
    ReduceDynamicsCK<MainExecutionPolicy, TemperatureVolWeightedMeanReduceCK> temperature_mean_reduce(glass);

    auto &glass_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass);
    glass_state_recorder.addToWrite<Real>(glass, temperature_name);
    glass_state_recorder.addToWrite<Real>(glass, "Pressure");
    if (enable_joule_heat)
    {
        glass_state_recorder.addToWrite<Real>(glass, kJouleHeatField);
    }
    if (write_em_current_to_vtp)
    {
        glass_state_recorder.addToWrite<Vecd>(glass, flow_glass_em_names.j_real);
        glass_state_recorder.addToWrite<Vecd>(glass, flow_glass_em_names.j_imag);
    }

    auto &rotor_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(rotor);
    rotor_state_recorder.addToWrite<Real>(rotor, temperature_name);

    StateDynamics<MainExecutionPolicy, UpdateSpinningParticlePosition> update_rotor_proxy_positions(
        rotor_proxy, Vec3d(0.0, 0.0, 0.0), Vec3d(0.0, 0.0, 1.0), kRotationSpeed);
    BodyStatesRecordingToTriangleMeshVtpCK<MainExecutionPolicy> write_rotor_surface(rotor_proxy, rotor_shape);
    Rh200SynchronizedAnimationWriter animation_writer;
    if (cli.state_recording)
    {
        animation_writer.open("./output/rh200_animation_manifest.csv", glass.Name(), rotor_proxy.Name(), rotor.Name());
    }

    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);

    SolidBodyPartForSimbody rotor_constraint_area(rotor, makeShared<Rotor>("Rotor", cli.geometry_scale));
    SimTK::Body::Rigid info(*rotor_constraint_area.body_part_mass_properties_);

    const Vecd rotor_axis(0.0, 0.0, 1.0);
    const SimTK::Vec3 normalized_rotor_axis = EigenToSimTK(rotor_axis.normalized());
    const SimTK::UnitVec3 unit_rotor_axis(static_cast<SimTK::Real>(normalized_rotor_axis[0]),
                                          static_cast<SimTK::Real>(normalized_rotor_axis[1]),
                                          static_cast<SimTK::Real>(normalized_rotor_axis[2]));
    SimTK::Rotation R_Rotor;
    R_Rotor.setRotationFromOneAxis(unit_rotor_axis, SimTK::ZAxis);
    const SimTK::Vec3 rotor_screw_origin(0.0, 0.0, 0.0);
    const SimTK::Transform offset_transform(R_Rotor, rotor_screw_origin);
    const SimTK::Transform child_transform(R_Rotor, SimTK::Vec3(0.0));
    SimTK::MobilizedBody::Pin mob_body_rotor(matter.Ground(), offset_transform, info, child_transform);
    MBsystem.realizeTopology();

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    StateDynamics<MainExecutionPolicy, solid_dynamics::ConstraintBodyPartBySimBodyCK> constraint_rotation_rotor(
        rotor_constraint_area, MBsystem, mob_body_rotor, integ);

    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    update_glass_cell_linked_list.exec();
    update_rotor_cell_linked_list.exec();
    update_wall_cell_linked_list.exec();
    constant_gravity.exec();
    update_glass_relations.exec();
    update_rotor_fluid_contact.exec();
    density_regularization.exec();
    fluid_advection_step_setup.exec();
    initialize_glass_temperature.exec();
    initialize_rotor_temperature.exec();
    if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
#if SPHINXSYS_USE_SYCL
        if (sample_joule_heat_from_grid)
        {
            sample_joule_heat_from_grid->exec();
        }
        if (sample_current_from_grid)
        {
            sample_current_from_grid->exec();
        }
#endif
    }
    if (cli.state_recording)
    {
        animation_writer.write(Real(0), glass_state_recorder, update_rotor_proxy_positions, write_rotor_surface,
                               &rotor_state_recorder);
    }

    SimTK::State state = MBsystem.getDefaultState();
    mob_body_rotor.setOneU(state, 0, kRotationSpeed);
    integ.setAccuracy(1e-5);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    size_t advection_steps = 0;
    const int screening_interval = 10;
    const Real state_record_interval = std::max(cli.state_record_interval, Real(1.0e-6));
    auto &state_recording = time_stepper.addTriggerByInterval(state_record_interval);
    const Real energy_log_interval = state_record_interval;
    auto &energy_logging = time_stepper.addTriggerByInterval(energy_log_interval);
    Real integrated_joule_energy = Real(0);
    Real t_mean_initial_record = initial_temperature_glass;

    if (enable_joule_heat)
    {
        if (cli.joule_mode == FakeJouleMode::EmGrid)
        {
            writeEnergyBudgetGridCsvHeader(kEnergyBudgetCsv);
        }
        else
        {
            writeEnergyBudgetCsvHeader(kEnergyBudgetCsv);
        }
        writeRh200HeatingRateAuditCsvHeader(kHeatingRateAuditCsv);
        if (cli.joule_mode == FakeJouleMode::EmGrid)
        {
#if SPHINXSYS_USE_SYCL
            if (sample_joule_heat_from_grid)
            {
                sample_joule_heat_from_grid->exec();
            }
            if (sample_current_from_grid)
            {
                sample_current_from_grid->exec();
            }
#endif
        }
        const Real p0 = fake_joule_power_reduce.exec();
        const Real e0 = thermal_energy_reduce.exec();
        const Real t_mean_initial = temperature_mean_reduce.exec();
        t_mean_initial_record = t_mean_initial;
        const Real t_min_initial = hostGlassTemperatureMin(glass.getBaseParticles(), temperature_name);
        const Real t_max_initial = temperature_max_reduce.exec();
        if (cli.joule_mode == FakeJouleMode::EmGrid)
        {
#if SPHINXSYS_USE_SYCL
            const Real n_glass = static_cast<Real>(glass.getBaseParticles().TotalRealParticles());
            const Real out_count = out_of_grid_joule_reduce ? out_of_grid_joule_reduce->exec() : Real(0);
            appendEnergyBudgetGridCsv(kEnergyBudgetCsv, Real(0), p0, e0, integrated_joule_energy,
                                      t_mean_initial, t_max_initial, out_count, out_count / (n_glass + TinyReal));
            appendRh200HeatingRateAuditCsv(
                kHeatingRateAuditCsv,
                makeRh200HeatingRateAuditRecord(Real(0), glass_thermal_mass, joule_target_power, p0,
                                                integrated_joule_energy, e0, t_mean_initial, t_min_initial,
                                                t_max_initial, t_mean_initial, out_count,
                                                out_count / (n_glass + TinyReal)));
#endif
        }
        else
        {
            appendEnergyBudgetCsv(kEnergyBudgetCsv, Real(0), p0, e0, integrated_joule_energy, t_mean_initial,
                                  t_max_initial);
            appendRh200HeatingRateAuditCsv(
                kHeatingRateAuditCsv,
                makeRh200HeatingRateAuditRecord(Real(0), glass_thermal_mass, joule_target_power, p0,
                                                integrated_joule_energy, e0, t_mean_initial, t_min_initial,
                                                t_max_initial, t_mean_initial, Real(0), Real(0)));
        }
    }

    std::cout << "[rh200] flow run: dp=" << cli.dp << " geometry_scale=" << cli.geometry_scale
              << " end_time=" << cli.end_time << " n_glass=" << glass.getBaseParticles().TotalRealParticles()
              << " n_wall=" << wall_boundary.getBaseParticles().TotalRealParticles()
              << " n_rotor=" << rotor.getBaseParticles().TotalRealParticles()
              << " (reload must match --dp; wall extrude + lattice use same dp)"
              << " thermal: Joule + glass/rotor diffusion (T0=" << initial_temperature_glass
              << " K glass&rotor; wall adiabatic; no stress/strain output)"
              << " joule_mode=" << fakeJouleModeName(cli.joule_mode)
              << " state_record_interval=" << state_record_interval << " s"
              << " (~" << rh200EstimatedStateRecordFrames(cli.end_time, state_record_interval)
              << " VTP frames for end_time=" << cli.end_time << " s";
    if (cli.state_record_frames_per_rev_user_set && cli.state_record_frames_per_rev > 0)
    {
        std::cout << ", dense " << cli.state_record_frames_per_rev << " frames/rev";
    }
    else
    {
        std::cout << ", staggered ~1/rev +" << cli.state_record_phase_slip_deg
                  << " deg rotor slip/frame (~" << rh200StateRecordIntervalDegreesPerFrame(state_record_interval)
                  << " deg total/frame)";
    }
    std::cout << ", rotor ~100 RPM)" << std::endl;

    while (!time_stepper.isEndTime(cli.end_time))
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        pressure_force_from_fluid_on_rotor.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);

        integ.stepBy(acoustic_dt);
        SimTK::State &state_for_update = integ.updAdvancedState();
        force_on_bodies.clearAllBodyForces(state_for_update);
        constraint_rotation_rotor.exec();

        if (enable_joule_heat)
        {
            if (cli.joule_mode == FakeJouleMode::EmGrid)
            {
#if SPHINXSYS_USE_SYCL
                if (sample_joule_heat_from_grid)
                {
                    sample_joule_heat_from_grid->exec();
                }
                if (sample_current_from_grid)
                {
                    sample_current_from_grid->exec();
                }
#endif
            }
            apply_fake_joule_heat.exec(acoustic_dt);
            integrated_joule_energy += fake_joule_power_reduce.exec() * acoustic_dt;
        }

        glass_thermal_relaxation.exec(acoustic_dt);
        rotor_thermal_relaxation.exec(acoustic_dt);

        if (advection_step(advection_time_step))
        {
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(6) << "[rh200] N=" << advection_steps
                          << " t=" << time_stepper.getPhysicalTime()
                          << " advection_dt=" << advection_step.getInterval()
                          << " acoustic_dt=" << time_stepper.getGlobalTimeStepSize() << std::endl;
            }
            ++advection_steps;
            fluid_update_particle_position.exec();

            const bool should_record = state_recording();
            if (should_record && cli.state_recording)
            {
                animation_writer.write(time_stepper.getPhysicalTime(), glass_state_recorder,
                                       update_rotor_proxy_positions, write_rotor_surface, &rotor_state_recorder);
            }

            if (enable_joule_heat && energy_logging())
            {
                const Real p_joule = fake_joule_power_reduce.exec();
                const Real t_mean = temperature_mean_reduce.exec();
                const Real t_max = temperature_max_reduce.exec();
                const Real t_min = hostGlassTemperatureMin(glass.getBaseParticles(), temperature_name);
                const Real thermal_energy = thermal_energy_reduce.exec();
                if (cli.joule_mode == FakeJouleMode::EmGrid)
                {
#if SPHINXSYS_USE_SYCL
                    const Real n_glass = static_cast<Real>(glass.getBaseParticles().TotalRealParticles());
                    const Real out_count = out_of_grid_joule_reduce ? out_of_grid_joule_reduce->exec() : Real(0);
                    appendEnergyBudgetGridCsv(kEnergyBudgetCsv, time_stepper.getPhysicalTime(), p_joule,
                                              thermal_energy, integrated_joule_energy, t_mean, t_max, out_count,
                                              out_count / (n_glass + TinyReal));
                    appendRh200HeatingRateAuditCsv(
                        kHeatingRateAuditCsv,
                        makeRh200HeatingRateAuditRecord(time_stepper.getPhysicalTime(), glass_thermal_mass,
                                                        joule_target_power, p_joule, integrated_joule_energy,
                                                        thermal_energy, t_mean, t_min, t_max, t_mean_initial_record,
                                                        out_count, out_count / (n_glass + TinyReal)));
#endif
                }
                else
                {
                    appendEnergyBudgetCsv(kEnergyBudgetCsv, time_stepper.getPhysicalTime(), p_joule, thermal_energy,
                                          integrated_joule_energy, t_mean, t_max);
                    appendRh200HeatingRateAuditCsv(
                        kHeatingRateAuditCsv,
                        makeRh200HeatingRateAuditRecord(time_stepper.getPhysicalTime(), glass_thermal_mass,
                                                        joule_target_power, p_joule, integrated_joule_energy,
                                                        thermal_energy, t_mean, t_min, t_max, t_mean_initial_record,
                                                        Real(0), Real(0)));
                }
            }

            if (advection_steps % 100 == 0 && advection_steps != 1)
            {
                fluid_particle_sorting.exec();
            }
            update_glass_cell_linked_list.exec();
            update_rotor_cell_linked_list.exec();
            update_glass_relations.exec();
            update_rotor_fluid_contact.exec();

            fluid_advection_step_setup.exec();
            density_regularization.exec();
            viscous_force_from_fluid_on_rotor.exec();
        }
    }

    animation_writer.close();
    std::cout << "[rh200] flow run finished at t=" << time_stepper.getPhysicalTime()
              << " animation_frames=" << animation_writer.frame_index << std::endl;
    return 0;
}

int runRh200EmSolveOnly(const Rh200GlassStirringCli &cli, int ac, char *av[])
{
    const StdVec<std::string> stl_files = {g_path_glass};
    MergedBoundingBox merger(stl_files, Vec3d::Zero(), cli.geometry_scale, cli.dp);
    const BoundingBoxd system_domain_bounds = merger.calculate();

    logStlBoundingBox("glass_fluid", g_path_glass, cli.geometry_scale);
    writeRh200GeometryAuditCsv("./output/rh200_geometry_audit.csv", cli.geometry_scale);

    SPHSystem sph_system(system_domain_bounds, cli.dp);
    sph_system.setReloadParticles(true);
    callRh200SphCommandLineOptions(sph_system, ac, av);
    IO::getEnvironment().resetReloadFolder(cli.reload_dir, true);

    logRh200MaterialPreset(cli.preset, cli.material, cli.em.sigma0, cli.em.frequency_hz, cli.em.target_power);
    const Real rho0_solid = cli.material.rho_glass;

    SolidBody glass(sph_system, makeShared<GlassGeometry>("GlassGeometry", cli.geometry_scale));
    glass.defineAdaptationRatios(1.3, 1.0);
    glass.defineMatterMaterial<Solid>(rho0_solid);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass, glass_names);
    (void)register_glass;

    Inner<> glass_inner(glass);

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    params.target_joule_power_ = cli.em.target_power;
    params.sigma_glass_ = cli.em.sigma0;
    params.frequency_ = cli.em.frequency_hz;

    OphelieTestCliOptions ophelie_cli;
    configureRh200OphelieEmParameters(params, ophelie_cli);
    logOphelieFinalParams(params, ophelie_cli);

    const Rh200GlassBboxAudit stl_bbox = measureStlGlassBbox(g_path_glass, cli.geometry_scale);
    const Rh200GlassBboxAudit particle_bbox = measureGlassParticleBbox(glass.getBaseParticles());
    fillFrenchParamsFromRh200Bbox(particle_bbox, cli.dp, cli.em, french, params);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    writeRh200CoilGeometryLogCsv(kRh200CoilGeometryCsv, cli.geometry_scale, stl_bbox, particle_bbox, french);
    appendRh200GlassParticleBboxToAuditCsv("./output/rh200_geometry_audit.csv", particle_bbox);

    std::cout << "[rh200] EM solve setup: dp=" << cli.dp << " n_glass=" << particle_bbox.particle_count
              << " glass_center=(" << particle_bbox.center.transpose() << ") R_xy=" << particle_bbox.radial_extent_xy
              << " coil_R=" << french.coil.loop_radius << " coil_z=[" << french.coil.z_min << "," << french.coil.z_max
              << "] turns=" << french.coil.num_loops << " target_power=" << cli.em.target_power << " W sigma="
              << cli.em.sigma0 << " f=" << cli.em.frequency_hz << " Hz" << std::endl;

    Rh200EmParticleDepositionHost em_dep_raw;
    Rh200EmParticleDepositionHost em_dep;
    const Rh200EmSolveResult em_result = runRh200OphelieEmSolveOnce<MainExecutionPolicy>(
        glass, glass_inner, glass_names, params, french, &em_dep_raw);
    hostExtractEmDepositionFields(glass.getBaseParticles(), glass_names, em_dep);

    Rh200JouleHeatGridSpec grid_spec;
    const Rh200JouleHeatGridSpec *grid_spec_ptr = nullptr;
    Real grid_rescale_factor = Real(1);
    Real p_scaled_grid_sample = em_result.joule_power_w;
    if (cli.joule_mode == FakeJouleMode::EmGrid)
    {
        Rh200JouleHeatEmGridBundle grid_bundle;
        const Rh200JouleHeatGridBuildReport grid_report =
            buildRh200EulerianJouleHeatGridFromEm(cli, grid_bundle, em_dep, &grid_spec);
        grid_spec_ptr = &grid_spec;
        grid_rescale_factor = grid_report.grid_scale_factor;
        p_scaled_grid_sample = grid_report.p_grid_sample_after_scale;
    }

    writeRh200ExcitationToJouleAuditFromEmSolve(cli, em_result, french, glass.getBaseParticles(), rho0_solid,
                                                  cli.material.cp_glass, &em_dep_raw, &em_dep, grid_spec_ptr,
                                                  grid_rescale_factor, p_scaled_grid_sample);
    writeRh200EmSolveSummaryCsv(kRh200EmSolveCsv, em_result, french);

    BodyStatesRecordingToVtpCK<MainExecutionPolicy> em_vtp_recorder(sph_system);
    em_vtp_recorder.addToWrite<Real>(glass, glass_names.sigma);
    em_vtp_recorder.addToWrite<Real>(glass, glass_names.joule_heat);
    em_vtp_recorder.addToWrite<Vecd>(glass, glass_names.e_real);
    em_vtp_recorder.addToWrite<Vecd>(glass, glass_names.e_imag);
    em_vtp_recorder.addToWrite<Vecd>(glass, glass_names.j_real);
    em_vtp_recorder.addToWrite<Vecd>(glass, glass_names.j_imag);
    em_vtp_recorder.addToWrite<Real>(glass, glass_names.phi_real);
    em_vtp_recorder.addToWrite<Real>(glass, glass_names.phi_imag);
    em_vtp_recorder.writeToFile();

    const Real power_rel_err =
        std::abs(em_result.joule_power_w - cli.em.target_power) / (cli.em.target_power + TinyReal);
    const bool em_ok = particle_bbox.particle_count > 0 && std::isfinite(em_result.joule_power_w) &&
                       em_result.joule_power_w > TinyReal && em_result.joule_heat_min >= Real(0) &&
                       power_rel_err < Real(0.08);
    std::cout << "[rh200] EM solve acceptance: P_joule=" << em_result.joule_power_w
              << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << em_result.phi_eq_res_vol
              << " em_ok=" << (em_ok ? 1 : 0) << std::endl;
    return em_ok ? 0 : 1;
}

} // namespace

int main(int ac, char *av[])
{
    Rh200GlassStirringCli cli;
    parseRh200GlassStirringCli(ac, av, cli);

    if (cli.relax)
    {
        return runParticleRelaxation(cli, ac, av);
    }

    if (cli.run_flow && cli.use_reload)
    {
        return runGlassStirringFlow(cli, ac, av);
    }

    if (cli.em_solve && cli.use_reload)
    {
        return runRh200EmSolveOnly(cli, ac, av);
    }

    std::cout
        << "test_3d_ophelie_rh200_glass_em_stirring (RH200 glass stirring + thermal + OPHELIE EM)\n"
        << "  --relax=1 [--dp=0.008] [--geometry-scale=1.0] [--relax-steps=5000]\n"
        << "  --reload=1 --em-solve=1 [--target-power=50000] [--sigma0=16] [--frequency=300000]\n"
        << "      [--coil-radius-factor=1.15] [--coil-turns=6] [--coil-z-margin-low=0.05] "
           "[--coil-z-margin-high=0.05]\n"
        << "  --reload=1 --run=1 [--end-time=1] [--joule-mode=off|em-fixed|em-grid]\n"
        << "      em-grid (recommended): auto EM solve + Eulerian JouleHeat grid on full_oil.stl, sample Q(x_i(t))\n"
        << "      em-fixed (debug): OPHELIE solve once then frozen particle JouleHeat + stirring\n"
        << "  [--joule-grid-spacing-factor=1.0]  grid spacing = factor * dp (1.0=dp, 1.5=1.5*dp)\n"
        << "  [--joule-grid-bbox-margin=2.0]  expand full_oil bbox by margin*spacing each side\n"
        << "  [--state-record-phase-slip-deg=24]  default: ~1 VTP/rev + 24 deg slip (smooth video, ~100 frames/60s)\n"
        << "  [--state-record-frames-per-rev=N]  optional dense mode (many VTPs; e.g. 15 => ~900 frames/60s)\n"
        << "  [--state-record-interval=SEC]  override interval in seconds\n"
        << "  [--target-power=50000] [--sigma0=16] [--frequency=300000]\n"
        << "Resolution: --dp applies to SPHSystem + wall shell thickness; relax/reload/em/run must use the SAME dp.\n"
        << "  Default dp=0.008. Re-relax if you change dp: --relax=1 --dp=0.008\n"
        << "EM solve outputs: ./output/rh200_coil_geometry_log.csv, rh200_em_solve_summary.csv, glass VTP with "
           "JouleHeat/E/J.\n"
        << "em-grid outputs: ./output/JouleHeatGrid.vti (JouleHeat, PhiReal, PhiImag, EMagnitude), "
           "rh200_joule_heat_grid_summary.csv\n"
        << "Animation (--state-recording=1): glass + RotorProxy written with SAME ite_XXXXXXXX tag; "
           "see ./output/rh200_animation_manifest.csv\n"
        << "Thermal: glass inner diffusion + glass-rotor contact; wall adiabatic (no wall thermal contact).\n"
        << "STL units: meters (geometry-scale=1.0).\n";
    return 1;
}
