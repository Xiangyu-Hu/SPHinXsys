#ifndef APHI_TEAM7_NATIVE_RELOAD_GEOMETRY_HELPERS_H
#define APHI_TEAM7_NATIVE_RELOAD_GEOMETRY_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_team7_shapes.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_coil_source_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include "io_vtk.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** TEAM7 native geometry / reload (SI [m]; must match particle_generation_em after STL scale=0.001). */
struct AphiTeam7NativeGeometryConstants
{
    /** mm -> m at STL import (disk STL vertices are in mm). */
    static constexpr Real stl_scale = 0.001;
    static constexpr Real dp_reference_default = 0.006;
    static constexpr Real dp_0_default = dp_reference_default;
    static constexpr Real sigma_aluminum = 3.526e7;
    static constexpr Real sigma_air = 0.0;
    static constexpr Real sigma_coil_source = 0.0;

    static const Vec3d air_box_lower_small;
    static const Vec3d air_box_upper_small;
    static const std::string path_coil_stl;
    static const std::string path_plate_stl;
};

/** SI magnetic reluctivity nu = 1/mu0 for native reference path (not nu=1 nondimensional). */
inline Real team7NativeSiMagneticReluctivity()
{
    static constexpr Real mu0 = 4.0 * Pi * 1.0e-7;
    return 1.0 / mu0;
}

inline const Vec3d AphiTeam7NativeGeometryConstants::air_box_lower_small(-0.050, -0.050, -0.050);
inline const Vec3d AphiTeam7NativeGeometryConstants::air_box_upper_small(0.350, 0.350, 0.200);
inline const std::string AphiTeam7NativeGeometryConstants::path_coil_stl = "./input/coil.stl";
inline const std::string AphiTeam7NativeGeometryConstants::path_plate_stl = "./input/plate.stl";

class AphiTeam7NativeCoilShape : public ComplexShape
{
  public:
    explicit AphiTeam7NativeCoilShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(AphiTeam7NativeGeometryConstants::path_coil_stl, Vec3d::Zero(),
                                  AphiTeam7NativeGeometryConstants::stl_scale);
    }
};

class AphiTeam7NativePlateShape : public ComplexShape
{
  public:
    explicit AphiTeam7NativePlateShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(AphiTeam7NativeGeometryConstants::path_plate_stl, Vec3d::Zero(),
                                  AphiTeam7NativeGeometryConstants::stl_scale);
    }
};

class AphiTeam7NativeAirShape : public ComplexShape
{
  public:
    explicit AphiTeam7NativeAirShape(const std::string &shape_name, const Vec3d &box_lower, const Vec3d &box_upper)
        : ComplexShape(shape_name)
    {
        Vecd halfsize = 0.5 * (box_upper - box_lower);
        Vecd center = 0.5 * (box_lower + box_upper);
        Transform translation(center);
        add<GeometricShapeBox>(translation, halfsize, "OuterBoundary");
        subtract<TriangleMeshShapeSTL>(AphiTeam7NativeGeometryConstants::path_coil_stl, Vec3d::Zero(),
                                        AphiTeam7NativeGeometryConstants::stl_scale, "Coil");
        subtract<TriangleMeshShapeSTL>(AphiTeam7NativeGeometryConstants::path_plate_stl, Vec3d::Zero(),
                                        AphiTeam7NativeGeometryConstants::stl_scale, "Plate");
    }
};

/** Uniform impressed-current RHS on every particle of the source (coil) body. */
class AssignUniformImpressedCurrentRhsCK : public LocalDynamics
{
  public:
    AssignUniformImpressedCurrentRhsCK(SPHBody &sph_body, const AphiBlockNames &rhs_block, const Vecd &current_real,
                                       const Vecd &current_imag, Real amplitude)
        : LocalDynamics(sph_body),
          dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
          dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)),
          current_real_(current_real), current_imag_(current_imag), amplitude_(amplitude)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
              rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)), current_real_(encloser.current_real_),
              current_imag_(encloser.current_imag_), amplitude_(encloser.amplitude_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            rhs_a_real_[index_i] = amplitude_ * current_real_;
            rhs_a_imag_[index_i] = amplitude_ * current_imag_;
        }

      protected:
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        Vecd current_real_;
        Vecd current_imag_;
        Real amplitude_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    Vecd current_real_;
    Vecd current_imag_;
    Real amplitude_;
};

struct AphiTeam7NativeReloadBoundingBox
{
    Vecd lower = Vecd::Constant(std::numeric_limits<Real>::max());
    Vecd upper = Vecd::Constant(std::numeric_limits<Real>::lowest());
    bool valid = false;
};

inline AphiTeam7NativeReloadBoundingBox hostParticleBoundingBox(SPHBody &body)
{
    BaseParticles &particles = body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const size_t count = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    AphiTeam7NativeReloadBoundingBox box;
    if (count == 0)
    {
        return box;
    }
    box.valid = true;
    for (size_t i = 0; i != count; ++i)
    {
        box.lower = box.lower.cwiseMin(positions[i]);
        box.upper = box.upper.cwiseMax(positions[i]);
    }
    return box;
}

inline bool team7NativeReloadXmlHasSmoothingLengthRatio(const std::string &body_name,
                                                        const std::string &reload_dir = team7NativeWorkingReloadDirectory())
{
    std::ifstream input(reload_dir + "/" + body_name + "_rld.xml");
    if (!input)
    {
        return false;
    }
    std::string line;
    while (std::getline(input, line))
    {
        if (line.find("SmoothingLengthRatio") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

inline bool nativeReloadXmlExists(const std::string &body_name,
                                  const std::string &reload_dir = team7NativeWorkingReloadDirectory())
{
    const std::string path = reload_dir + "/" + body_name + "_rld.xml";
    return std::filesystem::exists(path);
}

inline bool allNativeReloadXmlExist()
{
    return nativeReloadXmlExists("Coil") && nativeReloadXmlExists("Plate") && nativeReloadXmlExists("Air");
}

struct AphiTeam7NativeGeometryAuditSummary
{
    bool reload_files_present = false;
    size_t air_particles = 0;
    size_t plate_particles = 0;
    size_t coil_particles = 0;
    AphiTeam7NativeReloadBoundingBox air_bbox{};
    AphiTeam7NativeReloadBoundingBox plate_bbox{};
    AphiTeam7NativeReloadBoundingBox coil_bbox{};
    std::string vtp_output_dir = "output";
};

inline void printVec3dMm(const char *label, const Vecd &value)
{
    std::cout << label << "_mm=" << (1.0e3 * value[0]) << "," << (1.0e3 * value[1]) << "," << (1.0e3 * value[2]);
}

inline void printTeam7NativeGeometryAudit(const char *test_name, const AphiTeam7NativeGeometryAuditSummary &summary)
{
    std::cout << test_name << " reload_files_present=" << (summary.reload_files_present ? 1 : 0)
              << " air_particles=" << summary.air_particles << " plate_particles=" << summary.plate_particles
              << " coil_particles=" << summary.coil_particles;
    if (summary.air_bbox.valid)
    {
        std::cout << " ";
        printVec3dMm("air_bbox_lower", summary.air_bbox.lower);
        std::cout << " ";
        printVec3dMm("air_bbox_upper", summary.air_bbox.upper);
    }
    if (summary.plate_bbox.valid)
    {
        std::cout << " ";
        printVec3dMm("plate_bbox_lower", summary.plate_bbox.lower);
        std::cout << " ";
        printVec3dMm("plate_bbox_upper", summary.plate_bbox.upper);
    }
    if (summary.coil_bbox.valid)
    {
        std::cout << " ";
        printVec3dMm("coil_bbox_lower", summary.coil_bbox.lower);
        std::cout << " ";
        printVec3dMm("coil_bbox_upper", summary.coil_bbox.upper);
    }
    std::cout << " probe_A1B1_y_mm=72 probe_A1B1_z_mm=34 probe_A2B2_y_mm=144 probe_A2B2_z_mm=34"
              << " vtp_output_dir=" << summary.vtp_output_dir << std::endl;
}

/** Fail if reload still looks like legacy mm coordinates (e.g. plate extent ~294 instead of ~0.294 m). */
inline bool team7NativeReloadCoordinatesLookSiMeters(const AphiTeam7NativeReloadBoundingBox &plate_bbox)
{
    if (!plate_bbox.valid)
    {
        return false;
    }
    constexpr Real max_extent_si_m = 2.0;
    const Vecd extent = plate_bbox.upper - plate_bbox.lower;
    return extent.maxCoeff() < max_extent_si_m;
}

inline bool team7NativeGeometryAuditPassed(const AphiTeam7NativeGeometryAuditSummary &summary)
{
    return summary.reload_files_present && summary.air_particles > 0 && summary.plate_particles > 0 &&
           summary.coil_particles > 0 && summary.air_bbox.valid && summary.plate_bbox.valid && summary.coil_bbox.valid &&
           team7NativeReloadCoordinatesLookSiMeters(summary.plate_bbox);
}

/** Framework adaptive bodies (same pattern as particle_generation_em / particle_relaxation). */
using Team7NativeAdaptiveCoilBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearSurface, SolidBody>;
using Team7NativeAdaptivePlateBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearSurface, SolidBody>;
using Team7NativeAdaptiveAirBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearInnerSurface, RealBody>;
using Team7NativeAdaptiveCoilInner = Inner<Relation<Team7NativeAdaptiveCoilBody>>;
using Team7NativeAdaptivePlateInner = Inner<Relation<Team7NativeAdaptivePlateBody>>;
using Team7NativeAdaptiveAirInner = Inner<Relation<Team7NativeAdaptiveAirBody>>;
using Team7AirToCoilContact = Contact<Relation<Team7NativeAdaptiveAirBody, Team7NativeAdaptiveCoilBody>>;
using Team7AirToPlateContact = Contact<Relation<Team7NativeAdaptiveAirBody, Team7NativeAdaptivePlateBody>>;
using Team7CoilToAirContact = Contact<Relation<Team7NativeAdaptiveCoilBody, Team7NativeAdaptiveAirBody>>;
using Team7PlateToAirContact = Contact<Relation<Team7NativeAdaptivePlateBody, Team7NativeAdaptiveAirBody>>;

/** Three-body contact case with body names matching particle_generation_em reload XML. */
struct AphiTeam7NativeReloadContactCase
{
    AphiTeam7NativeGeometryConfig geometry_config;
    SPHSystem sph_system;
    Team7NativeAdaptiveCoilBody *coil_body_ptr_ = nullptr;
    Team7NativeAdaptivePlateBody *plate_body_ptr_ = nullptr;
    Team7NativeAdaptiveAirBody *air_body_ptr_ = nullptr;
    Team7NativeAdaptiveCoilInner *coil_inner_ck_ = nullptr;
    Team7NativeAdaptivePlateInner *plate_inner_ck_ = nullptr;
    Team7NativeAdaptiveAirInner *air_inner_ck_ = nullptr;
    Team7AirToCoilContact *air_to_coil_contact_ck_ = nullptr;
    Team7AirToPlateContact *air_to_plate_contact_ck_ = nullptr;
    Team7CoilToAirContact *coil_to_air_ck_ = nullptr;
    Team7PlateToAirContact *plate_to_air_ck_ = nullptr;

    RealBody &air_body() { return *air_body_ptr_; }
    const RealBody &air_body() const { return *air_body_ptr_; }
    SolidBody &coil_body() { return *coil_body_ptr_; }
    const SolidBody &coil_body() const { return *coil_body_ptr_; }
    SolidBody &plate_body() { return *plate_body_ptr_; }
    const SolidBody &plate_body() const { return *plate_body_ptr_; }

    AphiTeam7NativeReloadContactCase(int ac, char *av[])
        : geometry_config(resolveTeam7NativeGeometryConfig(ac, av)),
          sph_system(BoundingBoxd(geometry_config.air_box_lower, geometry_config.air_box_upper),
                     team7NativeSphReferenceDp(geometry_config))
    {
        printTeam7NativeGeometryConfig(geometry_config);
        sph_system.setRunParticleRelaxation(false);
        sph_system.setReloadParticles(true);
        if (ac > 0)
        {
            StdVec<std::string> sph_arg_storage;
            StdVec<char *> sph_argv;
            filterSphinxsysCommandlineArgs(ac, av, sph_arg_storage, sph_argv);
            sph_system.handleCommandlineOptions(static_cast<int>(sph_argv.size()), sph_argv.data());
        }

        const Real solid_global_resolution = team7NativeSolidAdaptiveGlobalResolution(geometry_config);
        const int air_refinement_levels = team7NativeAirRefinementLevels(geometry_config);
        const int solid_refinement_levels = team7NativeSolidRefinementLevels(geometry_config);
        auto &inner_boundary_shape =
            sph_system.addShape<AphiTeam7NativeInnerBoundaryShape>("InnerBoundary", geometry_config.stl_scale);

        AphiTeam7NativeAdaptiveNearSurface solid_adaptation(solid_global_resolution, 1.15, 1.0, solid_refinement_levels);
        AphiTeam7NativeAdaptiveNearInnerSurface air_adaptation(
            team7NativeAirAdaptiveGlobalResolution(geometry_config), 1.15, 1.0, air_refinement_levels,
            &inner_boundary_shape);

        Team7NativeAdaptiveCoilBody &coil = sph_system.addAdaptiveBody<SolidBody, AphiTeam7NativeAdaptiveNearSurface>(
            solid_adaptation, makeShared<AphiTeam7NativeCoilShape>("Coil"));
        Team7NativeAdaptivePlateBody &plate = sph_system.addAdaptiveBody<SolidBody, AphiTeam7NativeAdaptiveNearSurface>(
            solid_adaptation, makeShared<AphiTeam7NativePlateShape>("Plate"));
        Team7NativeAdaptiveAirBody &air = sph_system.addAdaptiveBody<RealBody, AphiTeam7NativeAdaptiveNearInnerSurface>(
            air_adaptation,
            makeShared<AphiTeam7NativeAirShape>("Air", geometry_config.air_box_lower, geometry_config.air_box_upper));

        coil_body_ptr_ = &coil;
        plate_body_ptr_ = &plate;
        air_body_ptr_ = &air;

        auto setup_adaptive_body = [&](SPHBody &body, const std::string &reload_body_name) {
            body.defineMaterial<Solid>();
            if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
            {
                BaseParticles *particles = body.generateParticles<BaseParticles, Reload>(reload_body_name);
                if (team7NativeReloadXmlHasSmoothingLengthRatio(reload_body_name))
                {
                    particles->reloadExtraVariable<Real>("SmoothingLengthRatio");
                }
                else if (team7NativeRefinementActive(geometry_config))
                {
                    std::cerr << "error: reload \"" << reload_body_name
                              << "\" missing SmoothingLengthRatio; regenerate with particle_generation_em\n";
                    std::exit(1);
                }
                else
                {
                    particles->registerStateVariable<Real>("SmoothingLengthRatio");
                    particles->addEvolvingVariable<Real>("SmoothingLengthRatio");
                    Real *h_ratio = particles->getVariableDataByName<Real>("SmoothingLengthRatio");
                    for (size_t i = 0; i != particles->TotalRealParticles(); ++i)
                    {
                        h_ratio[i] = 1.0;
                    }
                }
            }
            else
            {
                body.defineBodyLevelSetShape(par_ck);
                body.generateParticles<BaseParticles, Lattice>();
            }
        };

        setup_adaptive_body(coil, "Coil");
        setup_adaptive_body(plate, "Plate");
        setup_adaptive_body(air, "Air");

        std::cout << "team7_native_reload adaptive_bodies=1 coil_plate_levels=" << solid_refinement_levels
                  << " air_levels=" << air_refinement_levels
                  << " air_particles=" << air.getBaseParticles().TotalRealParticles()
                  << " coil_particles=" << coil.getBaseParticles().TotalRealParticles()
                  << " plate_particles=" << plate.getBaseParticles().TotalRealParticles() << std::endl;

        coil_inner_ck_ = &sph_system.addInnerRelation(coil);
        plate_inner_ck_ = &sph_system.addInnerRelation(plate);
        air_inner_ck_ = &sph_system.addInnerRelation(air);
        air_to_coil_contact_ck_ = &sph_system.addContactRelation(air, StdVec<Team7NativeAdaptiveCoilBody *>{&coil});
        air_to_plate_contact_ck_ = &sph_system.addContactRelation(air, StdVec<Team7NativeAdaptivePlateBody *>{&plate});
        coil_to_air_ck_ = &sph_system.addContactRelation(coil, StdVec<Team7NativeAdaptiveAirBody *>{&air});
        plate_to_air_ck_ = &sph_system.addContactRelation(plate, StdVec<Team7NativeAdaptiveAirBody *>{&air});

        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
    }

    Team7NativeAdaptiveAirInner &air_inner() { return *air_inner_ck_; }
    Team7NativeAdaptiveCoilInner &coil_inner() { return *coil_inner_ck_; }
    Team7NativeAdaptivePlateInner &plate_inner() { return *plate_inner_ck_; }
    Team7AirToCoilContact &air_to_coil_contact() { return *air_to_coil_contact_ck_; }
    Team7AirToPlateContact &air_to_plate_contact() { return *air_to_plate_contact_ck_; }
    Team7CoilToAirContact &coil_to_air() { return *coil_to_air_ck_; }
    Team7PlateToAirContact &plate_to_air() { return *plate_to_air_ck_; }

    void updateRelations()
    {
        UpdateCellLinkedList<MainExecutionPolicy, Team7NativeAdaptiveAirBody> update_air_cell_linked_list(*air_body_ptr_);
        UpdateCellLinkedList<MainExecutionPolicy, Team7NativeAdaptiveCoilBody> update_coil_cell_linked_list(*coil_body_ptr_);
        UpdateCellLinkedList<MainExecutionPolicy, Team7NativeAdaptivePlateBody> update_plate_cell_linked_list(
            *plate_body_ptr_);
        UpdateRelation<MainExecutionPolicy, Team7NativeAdaptiveAirInner> update_air_inner(air_inner());
        UpdateRelation<MainExecutionPolicy, Team7NativeAdaptiveCoilInner> update_coil_inner(coil_inner());
        UpdateRelation<MainExecutionPolicy, Team7NativeAdaptivePlateInner> update_plate_inner(plate_inner());
        UpdateRelation<MainExecutionPolicy, Team7AirToCoilContact> update_air_to_coil(air_to_coil_contact());
        UpdateRelation<MainExecutionPolicy, Team7AirToPlateContact> update_air_to_plate(air_to_plate_contact());
        UpdateRelation<MainExecutionPolicy, Team7CoilToAirContact> update_coil_to_air(coil_to_air());
        UpdateRelation<MainExecutionPolicy, Team7PlateToAirContact> update_plate_to_air(plate_to_air());
        update_air_cell_linked_list.exec();
        update_coil_cell_linked_list.exec();
        update_plate_cell_linked_list.exec();
        update_air_inner.exec();
        update_coil_inner.exec();
        update_plate_inner.exec();
        update_air_to_coil.exec();
        update_air_to_plate.exec();
        update_coil_to_air.exec();
        update_plate_to_air.exec();
    }
};

inline void initializeTeam7NativeReloadAphiFields(AphiTeam7NativeReloadContactCase &case_setup,
                                                    const AphiVariableNames &names,
                                                    Real plate_sigma = AphiTeam7NativeGeometryConstants::sigma_aluminum)
{
    const Real nu_si = team7NativeSiMagneticReluctivity();
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_air(
        case_setup.air_body(), AphiTeam7NativeGeometryConstants::sigma_air, nu_si, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_coil(
        case_setup.coil_body(), AphiTeam7NativeGeometryConstants::sigma_coil_source, nu_si, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_plate(case_setup.plate_body(), plate_sigma,
                                                                                 nu_si, names);

    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_air_material(
        case_setup.air_body(), AphiTeam7NativeGeometryConstants::sigma_air, nu_si, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_coil_material(
        case_setup.coil_body(), AphiTeam7NativeGeometryConstants::sigma_coil_source, nu_si, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_plate_material(case_setup.plate_body(), plate_sigma,
                                                                                       nu_si, names.material);

    initialize_air.exec();
    initialize_coil.exec();
    initialize_plate.exec();
    set_air_material.exec();
    set_coil_material.exec();
    set_plate_material.exec();
}

inline AphiTeam7NativeGeometryAuditSummary runTeam7NativeGeometryAudit(int ac, char *av[])
{
    AphiTeam7NativeGeometryAuditSummary summary;
    summary.reload_files_present = allNativeReloadXmlExist();
    if (!summary.reload_files_present)
    {
        std::cout << "team7_native_geometry_audit missing reload XML under ./reload/ "
                     "(run particle_generation_em from bin/, then copy reload/ and input/ here)"
                  << std::endl;
        return summary;
    }

    AphiTeam7NativeReloadContactCase case_setup(ac, av);
    summary.air_particles = case_setup.air_body().getBaseParticles().TotalRealParticles();
    summary.plate_particles = case_setup.plate_body().getBaseParticles().TotalRealParticles();
    summary.coil_particles = case_setup.coil_body().getBaseParticles().TotalRealParticles();
    summary.air_bbox = hostParticleBoundingBox(case_setup.air_body());
    summary.plate_bbox = hostParticleBoundingBox(case_setup.plate_body());
    summary.coil_bbox = hostParticleBoundingBox(case_setup.coil_body());

    IOEnvironment io_environment(case_setup.sph_system);
    (void)io_environment;
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_vtp(case_setup.sph_system);
    write_vtp.writeToFile(0);
    summary.vtp_output_dir = "output";
    return summary;
}

struct AphiTeam7NativeReloadSourceDrivenSmokeSpec
{
    Real omega = 2.0 * Pi * 50.0;
    Real phi_gauge_penalty = 100.0;
    Real impressed_current_amplitude = 8.0;
    Real tolerance = 5.0e-3;
    /** Stage 10.17A: allow loose diagnostic convergence on first native reload mesh. */
    Real diagnostic_max_true_rel = 0.1;
    UnsignedInt restart_dimension = 30;
    UnsignedInt max_outer_iterations = 80;
    Real max_air_to_conductor_joule_ratio = 0.05;
    bool write_vtp = true;
    bool use_team7_coil_path_rhs = true;
    Team7CoilPathSourceSpec coil_path_source{};
    Real plate_sigma = AphiTeam7NativeGeometryConstants::sigma_aluminum;
    /** If true, skip plate Joule gate (vacuum / source-only B sanity). */
    bool vacuum_plate_nonconductive = false;
};

struct AphiTeam7NativeReloadSourceDrivenSmokeMetrics
{
    bool reload_files_present = false;
    bool converged = false;
    Real final_true_rel = 0.0;
    Real coil_rhs_norm = 0.0;
    Real plate_joule_integral = 0.0;
    Real air_joule_integral = 0.0;
    Real coil_joule_integral = 0.0;
    bool finite_fields = false;
    size_t air_particles = 0;
    size_t plate_particles = 0;
    size_t coil_particles = 0;
};

inline Real hostBodyJouleIntegral(SPHBody &body, const std::string &joule_variable_name)
{
    BaseParticles &particles = body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const size_t count = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto all_particles = [](const Vecd &) { return true; };
    return hostParticleRegionVolWeightedJoulePower(particles, positions, count, all_particles, joule_variable_name);
}

inline void runTeam7NativeContactJoulePipeline(AphiTeam7NativeReloadContactCase &case_setup,
                                               const AphiVariableNames &names,
                                               const AphiJouleHeatingFieldNames &joule_names, Real omega)
{
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Team7NativeAdaptiveAirInner>> air_grad_phi_inner(
        DynamicsArgs(case_setup.air_inner(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Team7AirToCoilContact>> air_grad_phi_coil(
        DynamicsArgs(case_setup.air_to_coil_contact(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Team7AirToPlateContact>> air_grad_phi_plate(
        DynamicsArgs(case_setup.air_to_plate_contact(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy,
                          AphiComputeScalarPhiGradientCK<Team7NativeAdaptiveCoilInner, Team7CoilToAirContact>>
        coil_grad_phi(DynamicsArgs(case_setup.coil_inner(), names.solution, joule_names),
                      DynamicsArgs(case_setup.coil_to_air(), names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy,
                          AphiComputeScalarPhiGradientCK<Team7NativeAdaptivePlateInner, Team7PlateToAirContact>>
        plate_grad_phi(DynamicsArgs(case_setup.plate_inner(), names.solution, joule_names),
                       DynamicsArgs(case_setup.plate_to_air(), names.solution, joule_names));
    air_grad_phi_inner.exec();
    air_grad_phi_coil.exec();
    air_grad_phi_plate.exec();
    coil_grad_phi.exec();
    plate_grad_phi.exec();

    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> air_electric_field(
        case_setup.air_body(), omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> coil_electric_field(
        case_setup.coil_body(), omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> plate_electric_field(
        case_setup.plate_body(), omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> air_joule(case_setup.air_body(), names.material,
                                                                               joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> coil_joule(case_setup.coil_body(), names.material,
                                                                                joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> plate_joule(case_setup.plate_body(), names.material,
                                                                                 joule_names);

    air_electric_field.exec();
    coil_electric_field.exec();
    plate_electric_field.exec();
    air_joule.exec();
    coil_joule.exec();
    plate_joule.exec();
}

inline Team7CoilBoundingBox team7CoilBoundingBoxFromReload(const AphiTeam7NativeReloadBoundingBox &bbox)
{
    Team7CoilBoundingBox out;
    out.lower = bbox.lower;
    out.upper = bbox.upper;
    out.valid = bbox.valid;
    return out;
}

inline Real prepareTeam7NativeCoilPathRhs(AphiTeam7NativeReloadContactCase &case_setup,
                                          const AphiTeam7NativeReloadBoundingBox &coil_bbox,
                                          const Team7CoilPathSourceSpec &spec)
{
    const Team7CoilBoundingBox bbox_si = team7CoilBoundingBoxFromReload(coil_bbox);
    const StdVec<Vecd> centerline = buildTeam7CoilCenterlinePolyline(bbox_si, spec);
    hostStoreTeam7CoilSourceTangents(case_setup.coil_body(), centerline, spec);
    return team7CoilPathCurrentDensityJ0(bbox_si, spec, hostCoilVolume(case_setup.coil_body()));
}

inline Real hostBodyRhsBlockNorm(SPHBody &body, const AphiBlockNames &rhs_block)
{
    BaseParticles &particles = body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const size_t count = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto all_particles = [](const Vecd &) { return true; };
    return hostParticleRegionBlockNorm(particles, rhs_block, positions, count, all_particles);
}

template <typename InnerRelationType, typename ContactRelationType>
inline void team7NativeBindGmresApplyOverride(AphiMultiBodyContactEntry &entry, SPHBody &body,
                                              InnerRelationType &inner_relation, ContactRelationType &contact_relation,
                                              const AphiVariableNames &names, const AphiLhsAssemblyOptions &options)
{
    entry.apply_z_override = [&, inner_ptr = &inner_relation, contact_ptr = &contact_relation](
                                 const AphiBlockNames &input_block, const AphiBlockNames &output_block) {
        StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_output(body, output_block);
        InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<InnerRelationType>> apply_inner(
            DynamicsArgs(*inner_ptr, input_block, output_block, names.material, options.omega, options, Real(0.01)));
        InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<ContactRelationType>> apply_contact(DynamicsArgs(
            *contact_ptr, input_block, output_block, names.material, options.omega, options, Real(0.01)));
        zero_output.exec();
        apply_inner.exec();
        apply_contact.exec();
    };
}

template <typename InnerRelationType, typename ContactCoilType, typename ContactPlateType>
inline void team7NativeBindGmresAirApplyOverride(AphiMultiBodyContactEntry &entry, SPHBody &body,
                                                 InnerRelationType &inner_relation, ContactCoilType &air_to_coil,
                                                 ContactPlateType &air_to_plate, const AphiVariableNames &names,
                                                 const AphiLhsAssemblyOptions &options)
{
    entry.apply_z_override = [&, inner_ptr = &inner_relation, coil_ptr = &air_to_coil, plate_ptr = &air_to_plate](
                                 const AphiBlockNames &input_block, const AphiBlockNames &output_block) {
        StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_output(body, output_block);
        InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<InnerRelationType>> apply_inner(
            DynamicsArgs(*inner_ptr, input_block, output_block, names.material, options.omega, options, Real(0.01)));
        InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<ContactCoilType>> apply_coil(DynamicsArgs(
            *coil_ptr, input_block, output_block, names.material, options.omega, options, Real(0.01)));
        InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<ContactPlateType>> apply_plate(DynamicsArgs(
            *plate_ptr, input_block, output_block, names.material, options.omega, options, Real(0.01)));
        zero_output.exec();
        apply_inner.exec();
        apply_coil.exec();
        apply_plate.exec();
    };
}

template <typename InnerRelationType, typename ContactRelationType>
inline void team7NativeBindGmresJacobiOverride(AphiMultiBodyContactEntry &entry, SPHBody &body,
                                               InnerRelationType &inner_relation, ContactRelationType &contact_relation,
                                               const AphiVariableNames &names, const AphiLhsAssemblyOptions &options)
{
    entry.compute_jacobi_override = [&, inner_ptr = &inner_relation, contact_ptr = &contact_relation]() {
        InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<InnerRelationType>> compute_inner(
            DynamicsArgs(*inner_ptr, names.material, options.omega, options));
        InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<ContactRelationType>>
            compute_contact(DynamicsArgs(*contact_ptr, names.material, options.omega, options));
        compute_inner.exec();
        compute_contact.exec();
    };
}

template <typename InnerRelationType, typename ContactCoilType, typename ContactPlateType>
inline void team7NativeBindGmresAirJacobiOverride(AphiMultiBodyContactEntry &entry, SPHBody &body,
                                                  InnerRelationType &inner_relation, ContactCoilType &air_to_coil,
                                                  ContactPlateType &air_to_plate, const AphiVariableNames &names,
                                                  const AphiLhsAssemblyOptions &options)
{
    entry.compute_jacobi_override = [&, inner_ptr = &inner_relation, coil_ptr = &air_to_coil, plate_ptr = &air_to_plate]() {
        InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<InnerRelationType>> compute_inner(
            DynamicsArgs(*inner_ptr, names.material, options.omega, options));
        InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<ContactCoilType>> compute_coil(
            DynamicsArgs(*coil_ptr, names.material, options.omega, options));
        InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<ContactPlateType>> compute_plate(
            DynamicsArgs(*plate_ptr, names.material, options.omega, options));
        compute_inner.exec();
        compute_coil.exec();
        compute_plate.exec();
    };
}

inline StdVec<AphiMultiBodyContactEntry> buildTeam7NativeGmresEntries(AphiTeam7NativeReloadContactCase &case_setup,
                                                                      const AphiVariableNames &names,
                                                                      const AphiLhsAssemblyOptions &options)
{
    AphiLhsAssemblyOptions operator_options = options;
    operator_options.use_a_divergence_penalty = false;
    StdVec<AphiMultiBodyContactEntry> entries;
    AphiMultiBodyContactEntry air_entry{case_setup.air_body()};
    team7NativeBindGmresAirJacobiOverride(air_entry, case_setup.air_body(), case_setup.air_inner(),
                                          case_setup.air_to_coil_contact(), case_setup.air_to_plate_contact(), names,
                                          operator_options);
    team7NativeBindGmresAirApplyOverride(air_entry, case_setup.air_body(), case_setup.air_inner(),
                                         case_setup.air_to_coil_contact(), case_setup.air_to_plate_contact(), names,
                                         operator_options);
    entries.push_back(air_entry);

    AphiMultiBodyContactEntry coil_entry{case_setup.coil_body()};
    team7NativeBindGmresJacobiOverride(coil_entry, case_setup.coil_body(), case_setup.coil_inner(), case_setup.coil_to_air(),
                                       names, operator_options);
    team7NativeBindGmresApplyOverride(coil_entry, case_setup.coil_body(), case_setup.coil_inner(), case_setup.coil_to_air(),
                                      names, operator_options);
    entries.push_back(coil_entry);

    AphiMultiBodyContactEntry plate_entry{case_setup.plate_body()};
    team7NativeBindGmresJacobiOverride(plate_entry, case_setup.plate_body(), case_setup.plate_inner(),
                                       case_setup.plate_to_air(), names, operator_options);
    team7NativeBindGmresApplyOverride(plate_entry, case_setup.plate_body(), case_setup.plate_inner(),
                                      case_setup.plate_to_air(), names, operator_options);
    entries.push_back(plate_entry);
    return entries;
}

inline AphiGMRESResult runTeam7NativeReloadEmSolve(AphiTeam7NativeReloadContactCase &case_setup,
                                                   const AphiTeam7NativeReloadSourceDrivenSmokeSpec &spec,
                                                   const AphiVariableNames &names,
                                                   const AphiJouleHeatingFieldNames &joule_names)
{
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    initializeTeam7NativeReloadAphiFields(case_setup, names, spec.plate_sigma);

    const AphiGMRESWorkspaceNames gmres_workspace = buildAphiGMRESWorkspaceNames(spec.restart_dimension);
    auto register_solver_fields = [&](SPHBody &body) {
        RegisterAphiGMRESWorkspaceCK register_workspace(body, spec.restart_dimension);
        (void)register_workspace;
        RegisterAphiJouleHeatingFieldsCK register_joule(body, joule_names);
        (void)register_joule;
    };
    register_solver_fields(case_setup.air_body());
    register_solver_fields(case_setup.coil_body());
    register_solver_fields(case_setup.plate_body());

    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_rhs(case_setup.air_body(), names.rhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_rhs(case_setup.plate_body(), names.rhs);
    StateDynamics<MainExecutionPolicy, AssignUniformImpressedCurrentRhsCK> assign_uniform_coil_source(
        case_setup.coil_body(), names.rhs, coil_current_real, coil_current_imag, spec.impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_air_solution(case_setup.air_body(), names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_coil_solution(case_setup.coil_body(), names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_plate_solution(case_setup.plate_body(), names.solution);

    const AphiTeam7NativeReloadBoundingBox coil_bbox = hostParticleBoundingBox(case_setup.coil_body());
    const Real team7_j0 = prepareTeam7NativeCoilPathRhs(case_setup, coil_bbox, spec.coil_path_source);
    StateDynamics<MainExecutionPolicy, AssignTeam7CoilPathImpressedCurrentRhsCK> assign_team7_coil_source(
        case_setup.coil_body(), names.rhs, team7_j0);

    AphiLhsAssemblyOptions options;
    options.omega = spec.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = spec.phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(spec.tolerance, spec.restart_dimension, spec.max_outer_iterations);

    zero_air_rhs.exec();
    zero_plate_rhs.exec();
    if (spec.use_team7_coil_path_rhs)
    {
        assign_team7_coil_source.exec();
    }
    else
    {
        assign_uniform_coil_source.exec();
    }
    zero_air_solution.exec();
    zero_coil_solution.exec();
    zero_plate_solution.exec();
    case_setup.updateRelations();

    const StdVec<AphiMultiBodyContactEntry> entries = buildTeam7NativeGmresEntries(case_setup, names, options);
    AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy> solver(entries, names, gmres_workspace, options,
                                                                  solver_options.gmres);
    const AphiGMRESResult solver_result = solver.solve();
    case_setup.updateRelations();
    runTeam7NativeContactJoulePipeline(case_setup, names, joule_names, spec.omega);
    return solver_result;
}

/** SYCL/CUDA: pull solved fields to host before VTK write (otherwise VTP shows zeros). */
inline void team7NativeSyncPostSolveVtpFieldsToHost(SPHBody &body, const AphiVariableNames &names,
                                                    const AphiJouleHeatingFieldNames &joule_names, bool sync_joule)
{
    BaseParticles &particles = body.getBaseParticles();
    syncAphiBlockToHost(particles, names.solution);
    syncVariableToHost<Vecd>(particles, "NativeProbeBReal");
    syncVariableToHost<Vecd>(particles, "NativeProbeBImag");
    if (sync_joule)
    {
        syncVariableToHost<Vecd>(particles, joule_names.current_density_real);
        syncVariableToHost<Vecd>(particles, joule_names.current_density_imag);
        syncVariableToHost<Vecd>(particles, joule_names.electric_field_a_real);
    }
}

/**
 * Post-solve VTP for ParaView (TEAM7-style figures).
 * - Air/Coil/Plate: A, phi, curl-A B (NativeProbeB*), sigma
 * - Plate: induced current density J (real/imag) for Jy colormap + glyphs
 * - Coil: Team7CoilSourceTangent for source-current arrows
 */
inline void writeTeam7NativeReloadPostSolveVtp(AphiTeam7NativeReloadContactCase &case_setup,
                                               const AphiVariableNames &names,
                                               const AphiJouleHeatingFieldNames &joule_names)
{
    static constexpr const char *k_b_real = "NativeProbeBReal";
    static constexpr const char *k_b_imag = "NativeProbeBImag";
    team7NativeSyncPostSolveVtpFieldsToHost(case_setup.air_body(), names, joule_names, false);
    team7NativeSyncPostSolveVtpFieldsToHost(case_setup.coil_body(), names, joule_names, false);
    team7NativeSyncPostSolveVtpFieldsToHost(case_setup.plate_body(), names, joule_names, true);
    syncVariableToHost<Vecd>(case_setup.coil_body().getBaseParticles(), kTeam7CoilSourceTangentName);

    case_setup.air_body().setNewlyUpdated();
    case_setup.coil_body().setNewlyUpdated();
    case_setup.plate_body().setNewlyUpdated();

    BodyStatesRecordingToVtp write_states(case_setup.sph_system);
    const auto add_em_fields = [&](SPHBody &body) {
        write_states.addToWrite<Vecd>(body, names.solution.a_real);
        write_states.addToWrite<Vecd>(body, names.solution.a_imag);
        write_states.addToWrite<Real>(body, names.solution.phi_real);
        write_states.addToWrite<Real>(body, names.solution.phi_imag);
        write_states.addToWrite<Vecd>(body, k_b_real);
        write_states.addToWrite<Vecd>(body, k_b_imag);
        write_states.addToWrite<Real>(body, names.material.sigma);
    };
    add_em_fields(case_setup.air_body());
    add_em_fields(case_setup.coil_body());
    add_em_fields(case_setup.plate_body());
    write_states.addToWrite<Vecd>(case_setup.plate_body(), joule_names.current_density_real);
    write_states.addToWrite<Vecd>(case_setup.plate_body(), joule_names.current_density_imag);
    write_states.addToWrite<Vecd>(case_setup.plate_body(), joule_names.electric_field_a_real);
    write_states.addToWrite<Vecd>(case_setup.coil_body(), kTeam7CoilSourceTangentName);
    write_states.writeToFile(0);
}

inline AphiTeam7NativeReloadSourceDrivenSmokeMetrics runTeam7NativeReloadSourceDrivenSmoke(
    int ac, char *av[], const AphiTeam7NativeReloadSourceDrivenSmokeSpec &spec)
{
    AphiTeam7NativeReloadSourceDrivenSmokeMetrics metrics;
    metrics.reload_files_present = allNativeReloadXmlExist();
    if (!metrics.reload_files_present)
    {
        return metrics;
    }

    AphiTeam7NativeReloadContactCase case_setup(ac, av);
    metrics.air_particles = case_setup.air_body().getBaseParticles().TotalRealParticles();
    metrics.plate_particles = case_setup.plate_body().getBaseParticles().TotalRealParticles();
    metrics.coil_particles = case_setup.coil_body().getBaseParticles().TotalRealParticles();

    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    const AphiGMRESResult solver_result = runTeam7NativeReloadEmSolve(case_setup, spec, names, joule_names);

    metrics.converged = gmresConvergencePassed(solver_result, spec.tolerance);
    metrics.final_true_rel = solver_result.final_true_relative_residual;
    metrics.coil_rhs_norm = hostBodyRhsBlockNorm(case_setup.coil_body(), names.rhs);

    metrics.plate_joule_integral = hostBodyJouleIntegral(case_setup.plate_body(), joule_names.joule_heat_source);
    metrics.air_joule_integral = hostBodyJouleIntegral(case_setup.air_body(), joule_names.joule_heat_source);
    metrics.coil_joule_integral = hostBodyJouleIntegral(case_setup.coil_body(), joule_names.joule_heat_source);

    BaseParticles &plate_particles = case_setup.plate_body().getBaseParticles();
    const size_t plate_count = plate_particles.TotalRealParticles();
    metrics.finite_fields =
        hostVecdFieldFinite(plate_particles, names.solution.a_real, plate_count) &&
        hostVecdFieldFinite(plate_particles, joule_names.electric_field_a_real, plate_count) &&
        hostScalarFieldFinite(plate_particles, joule_names.joule_heat_source, plate_count);

    if (spec.write_vtp)
    {
        BodyStatesRecordingToVtp write_states(case_setup.sph_system);
        write_states.addToWrite<Vecd>(case_setup.plate_body(), names.solution.a_real);
        write_states.addToWrite<Vecd>(case_setup.plate_body(), names.solution.a_imag);
        write_states.addToWrite<Real>(case_setup.plate_body(), names.solution.phi_real);
        write_states.addToWrite<Real>(case_setup.plate_body(), names.solution.phi_imag);
        write_states.addToWrite<Vecd>(case_setup.plate_body(), joule_names.electric_field_a_real);
        write_states.addToWrite<Vecd>(case_setup.plate_body(), joule_names.current_density_real);
        write_states.addToWrite<Real>(case_setup.plate_body(), joule_names.joule_heat_source);
        write_states.addToWrite<Real>(case_setup.plate_body(), names.material.sigma);
        write_states.writeToFile(0);
    }
    return metrics;
}

inline bool team7NativeReloadSourceDrivenSmokePassed(const AphiTeam7NativeReloadSourceDrivenSmokeMetrics &metrics,
                                                     const AphiTeam7NativeReloadSourceDrivenSmokeSpec &spec)
{
    const Real air_to_conductor =
        metrics.air_joule_integral / (metrics.plate_joule_integral + TinyReal);
    const bool solver_ok = metrics.converged || metrics.final_true_rel <= spec.diagnostic_max_true_rel;
    const bool joule_ok = spec.vacuum_plate_nonconductive ||
                          (metrics.plate_joule_integral > 0.0 &&
                           air_to_conductor <= spec.max_air_to_conductor_joule_ratio);
    return metrics.reload_files_present && metrics.air_particles > 0 && metrics.plate_particles > 0 &&
           metrics.coil_particles > 0 && solver_ok && metrics.finite_fields && metrics.coil_rhs_norm > 0.0 && joule_ok;
}

inline bool writeTeam7NativeJoulePostSummaryCsv(const std::string &path,
                                                const AphiTeam7NativeReloadSourceDrivenSmokeMetrics &metrics)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    const Real air_to_conductor = metrics.air_joule_integral / (metrics.plate_joule_integral + TinyReal);
    output << std::setprecision(10);
    output << "metric,value\n";
    output << "plate_joule_integral," << metrics.plate_joule_integral << "\n";
    output << "air_joule_integral," << metrics.air_joule_integral << "\n";
    output << "coil_joule_integral," << metrics.coil_joule_integral << "\n";
    output << "air_to_conductor_joule," << air_to_conductor << "\n";
    output << "final_true_rel," << metrics.final_true_rel << "\n";
    output << "converged," << (metrics.converged ? 1 : 0) << "\n";
    return true;
}

inline void printTeam7NativeReloadSourceDrivenSmoke(const char *test_name,
                                                    const AphiTeam7NativeReloadSourceDrivenSmokeMetrics &metrics,
                                                    const AphiTeam7NativeReloadSourceDrivenSmokeSpec &spec, bool passed)
{
    const Real air_to_conductor =
        metrics.air_joule_integral / (metrics.plate_joule_integral + TinyReal);
    std::cout << test_name << " passed=" << (passed ? 1 : 0)
              << " reload_files_present=" << (metrics.reload_files_present ? 1 : 0)
              << " converged=" << (metrics.converged ? 1 : 0)
              << " diagnostic_converged=" << (metrics.final_true_rel <= spec.diagnostic_max_true_rel ? 1 : 0)
              << " final_true_rel=" << metrics.final_true_rel
              << " coil_rhs_norm=" << metrics.coil_rhs_norm
              << " plate_joule_integral=" << metrics.plate_joule_integral
              << " air_joule_integral=" << metrics.air_joule_integral
              << " coil_joule_integral=" << metrics.coil_joule_integral
              << " air_to_conductor_joule=" << air_to_conductor
              << " finite_fields=" << (metrics.finite_fields ? 1 : 0) << " air_particles=" << metrics.air_particles
              << " plate_particles=" << metrics.plate_particles << " coil_particles=" << metrics.coil_particles
              << std::endl;
}

struct AphiTeam7NativeSourceRhsAuditSummary
{
    bool reload_files_present = false;
    size_t coil_particles = 0;
    Real coil_volume = 0.0;
    Real path_length = 0.0;
    Real a_eff_particle = 0.0;
    Real a_eff_fallback = 0.0;
    Real j0_used = 0.0;
    Real ni_target = 2742.0;
    Real i_eff_formula = 0.0;
    Real i_eff_particle = 0.0;
    Real i_eff_ratio_formula = 0.0;
    Real i_eff_ratio_particle = 0.0;
    Real rhs_l2_team7 = 0.0;
    Real rhs_l2_uniform_legacy = 0.0;
    Real tangent_norm_min = 0.0;
    Real tangent_norm_max = 0.0;
    Real centerline_min_distance = 0.0;
    Real centerline_mean_distance = 0.0;
    bool audit_csv_written = false;
};

inline AphiTeam7NativeSourceRhsAuditSummary runTeam7NativeSourceRhsAudit(int ac, char *av[],
                                                                         const Team7CoilPathSourceSpec &spec)
{
    AphiTeam7NativeSourceRhsAuditSummary summary;
    summary.reload_files_present = allNativeReloadXmlExist();
    summary.ni_target = spec.totalAmpereTurns();
    summary.a_eff_fallback = spec.fallback_cross_section_area;
    if (!summary.reload_files_present)
    {
        return summary;
    }

    AphiTeam7NativeReloadContactCase case_setup(ac, av);
    summary.coil_particles = case_setup.coil_body().getBaseParticles().TotalRealParticles();
    const AphiTeam7NativeReloadBoundingBox coil_bbox = hostParticleBoundingBox(case_setup.coil_body());
    const Team7CoilBoundingBox bbox_si = team7CoilBoundingBoxFromReload(coil_bbox);
    const StdVec<Vecd> centerline = buildTeam7CoilCenterlinePolyline(bbox_si, spec);
    summary.path_length = closedPolylinePathLength(centerline);
    summary.coil_volume = hostCoilVolume(case_setup.coil_body());
    summary.a_eff_particle = spec.use_particle_volume_cross_section && summary.path_length > TinyReal
                                 ? summary.coil_volume / summary.path_length
                                 : spec.fallback_cross_section_area;
    summary.j0_used = spec.source_scale * summary.ni_target / (summary.a_eff_particle + TinyReal);
    summary.i_eff_formula = summary.j0_used * summary.a_eff_particle;
    summary.i_eff_ratio_formula = summary.i_eff_formula / summary.ni_target;

    AphiVariableNames names;
    const Real nu_si = team7NativeSiMagneticReluctivity();
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_coil_for_audit(
        case_setup.coil_body(), AphiTeam7NativeGeometryConstants::sigma_coil_source, nu_si, names);
    initialize_coil_for_audit.exec();

    hostStoreTeam7CoilSourceTangents(case_setup.coil_body(), centerline, spec);
    syncVariableToHost<Vecd>(case_setup.coil_body().getBaseParticles(), kTeam7CoilSourceTangentName);
    const Vecd *tangents = case_setup.coil_body().getBaseParticles().getVariableDataByName<Vecd>(kTeam7CoilSourceTangentName);
    const size_t count = case_setup.coil_body().getBaseParticles().TotalRealParticles();
    syncVariableToHost<Vecd>(case_setup.coil_body().getBaseParticles(), "Position");
    const Vecd *positions = case_setup.coil_body().getBaseParticles().getVariableDataByName<Vecd>("Position");

    summary.tangent_norm_min = std::numeric_limits<Real>::max();
    summary.tangent_norm_max = 0.0;
    Real distance_sum = 0.0;
    summary.centerline_min_distance = std::numeric_limits<Real>::max();
    for (size_t i = 0; i != count; ++i)
    {
        const Real tangent_norm = tangents[i].norm();
        summary.tangent_norm_min = std::min(summary.tangent_norm_min, tangent_norm);
        summary.tangent_norm_max = std::max(summary.tangent_norm_max, tangent_norm);
        Real best_distance_squared = std::numeric_limits<Real>::max();
        for (size_t k = 0; k < centerline.size(); ++k)
        {
            const Vecd &a = centerline[k];
            const Vecd &b = centerline[(k + 1) % centerline.size()];
            const Vecd ab = b - a;
            const Real ab_squared = ab.squaredNorm();
            if (ab_squared <= TinyReal)
            {
                continue;
            }
            Real segment_parameter = (positions[i] - a).dot(ab) / ab_squared;
            segment_parameter = std::max(Real(0), std::min(Real(1), segment_parameter));
            const Vecd closest = a + segment_parameter * ab;
            best_distance_squared = std::min(best_distance_squared, (positions[i] - closest).squaredNorm());
        }
        const Real nearest = std::sqrt(best_distance_squared);
        summary.centerline_min_distance = std::min(summary.centerline_min_distance, nearest);
        distance_sum += nearest;
    }
    if (count > 0)
    {
        summary.centerline_mean_distance = distance_sum / static_cast<Real>(count);
    }

    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_rhs(case_setup.coil_body(), names.rhs);
    zero_rhs.exec();
    StateDynamics<MainExecutionPolicy, AssignTeam7CoilPathImpressedCurrentRhsCK> assign_team7(
        case_setup.coil_body(), names.rhs, summary.j0_used);
    assign_team7.exec();
    summary.rhs_l2_team7 = hostBodyRhsBlockNorm(case_setup.coil_body(), names.rhs);
    summary.i_eff_particle = hostCoilIntegratedCurrentFromRhs(case_setup.coil_body(), names.rhs, summary.path_length);
    summary.i_eff_ratio_particle = summary.i_eff_particle / summary.ni_target;

    zero_rhs.exec();
    const Vecd legacy_direction(0.0, 0.0, 1.0);
    StateDynamics<MainExecutionPolicy, AssignUniformImpressedCurrentRhsCK> assign_legacy(
        case_setup.coil_body(), names.rhs, legacy_direction, Vecd::Zero(), 8.0);
    assign_legacy.exec();
    summary.rhs_l2_uniform_legacy = hostBodyRhsBlockNorm(case_setup.coil_body(), names.rhs);

    std::error_code mkdir_error;
    std::filesystem::create_directories("team7_native_source_rhs_audit", mkdir_error);
    std::ofstream audit_csv("team7_native_source_rhs_audit/summary.csv");
    if (audit_csv)
    {
        audit_csv << std::setprecision(10);
        audit_csv << "metric,value\n";
        audit_csv << "coil_particles," << summary.coil_particles << "\n";
        audit_csv << "coil_volume," << summary.coil_volume << "\n";
        audit_csv << "path_length," << summary.path_length << "\n";
        audit_csv << "a_eff_particle," << summary.a_eff_particle << "\n";
        audit_csv << "a_eff_fallback," << summary.a_eff_fallback << "\n";
        audit_csv << "j0_used," << summary.j0_used << "\n";
        audit_csv << "i_eff_formula," << summary.i_eff_formula << "\n";
        audit_csv << "i_eff_particle," << summary.i_eff_particle << "\n";
        audit_csv << "i_eff_ratio_formula," << summary.i_eff_ratio_formula << "\n";
        audit_csv << "i_eff_ratio_particle," << summary.i_eff_ratio_particle << "\n";
        audit_csv << "rhs_l2_team7," << summary.rhs_l2_team7 << "\n";
        audit_csv << "rhs_l2_uniform_legacy," << summary.rhs_l2_uniform_legacy << "\n";
        audit_csv << "tangent_norm_min," << summary.tangent_norm_min << "\n";
        audit_csv << "tangent_norm_max," << summary.tangent_norm_max << "\n";
        audit_csv << "centerline_min_distance_m," << summary.centerline_min_distance << "\n";
        audit_csv << "centerline_mean_distance_m," << summary.centerline_mean_distance << "\n";
        summary.audit_csv_written = true;
    }
    return summary;
}

inline bool team7NativeSourceRhsAuditPassed(const AphiTeam7NativeSourceRhsAuditSummary &summary)
{
    constexpr Real ratio_low = 0.9;
    constexpr Real ratio_high = 1.1;
    const bool current_ok = summary.i_eff_ratio_formula > ratio_low && summary.i_eff_ratio_formula < ratio_high &&
                            summary.i_eff_ratio_particle > ratio_low && summary.i_eff_ratio_particle < ratio_high;
    const bool rhs_ok = summary.rhs_l2_team7 > 0.0;
    const bool tangent_ok = summary.tangent_norm_min > 0.5 && summary.tangent_norm_max < 1.5;
    return summary.reload_files_present && summary.coil_particles > 0 && summary.path_length > TinyReal &&
           summary.coil_volume > TinyReal && current_ok && rhs_ok && tangent_ok && summary.audit_csv_written;
}

inline void printTeam7NativeSourceRhsAudit(const char *test_name, const AphiTeam7NativeSourceRhsAuditSummary &summary,
                                           bool passed)
{
    std::cout << test_name << " passed=" << (passed ? 1 : 0)
              << " reload_files_present=" << (summary.reload_files_present ? 1 : 0)
              << " coil_particles=" << summary.coil_particles << " coil_volume=" << summary.coil_volume
              << " path_length=" << summary.path_length << " a_eff_particle=" << summary.a_eff_particle
              << " j0_used=" << summary.j0_used << " i_eff_formula=" << summary.i_eff_formula
              << " i_eff_particle=" << summary.i_eff_particle
              << " i_eff_ratio_formula=" << summary.i_eff_ratio_formula
              << " i_eff_ratio_particle=" << summary.i_eff_ratio_particle << " rhs_l2_team7=" << summary.rhs_l2_team7
              << " rhs_l2_uniform_legacy=" << summary.rhs_l2_uniform_legacy
              << " tangent_norm_min=" << summary.tangent_norm_min
              << " tangent_norm_max=" << summary.tangent_norm_max
              << " centerline_min_distance_m=" << summary.centerline_min_distance
              << " centerline_mean_distance_m=" << summary.centerline_mean_distance << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_RELOAD_GEOMETRY_HELPERS_H
