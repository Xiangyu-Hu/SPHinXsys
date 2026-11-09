/**
 * @file  two_phase_dambreak.cpp
 * @brief GPU paper A4: 2D two-phase dambreak (extra_source_and_tests copy of CK source).
 * @details Duration V2 End_Time=10 (tau=10 with H=g=1); probe + gated VTP; no regression gate.
 * @note    Geometry is current CK source (DL=5.3H, DH=2H, probe y=0.2H), not paper-matched.
 */
#include "gpu_paper_recording.h"
#include "near_wall_boundary.hpp"
#include "sphinxsys.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace SPH;
using namespace sphinxsys_gpu_paper;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 10.0; /**< Duration V2 default; override via GPU_PAPER_END_TIME */
Real DL = 5.3;                   /**< Tank length. */
Real DH = 2.0;                   /**< Tank height. */
Real LL = 2.0;                   /**< Liquid column length. */
Real LH = 1.0;                   /**< Liquid column height. */
Real lh_over_dp = 20.0;          /**< L1=20→dp=0.05, L2=40, L3=80 */
Real particle_spacing_ref = LH / lh_over_dp;
Real BW = particle_spacing_ref * 4;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real rho0_a = 0.001;
Real gravity_g = 1.0;
Real U_ref = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_ref;

void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        lh_over_dp = Real(std::atof(v));
    if (lh_over_dp < Real(5.0))
        lh_over_dp = Real(20.0);
    particle_spacing_ref = LH / lh_over_dp;
    BW = particle_spacing_ref * 4;
    system_domain_bounds = BoundingBoxd(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    observation_location = {Vecd(DL, 0.2)};
}

namespace SPH
{
namespace fluid_dynamics
{

template <typename...>
class BuoyancyForce;

template <typename... Parameters>
class BuoyancyForce<Contact<WithUpdate, Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using BaseDynamicsType = Interaction<Contact<Parameters...>>;

  public:
    BuoyancyForce(Contact<Parameters...> &water_air_relation, const Gravity &gravity, Real threshold = 0.3)
        : BaseDynamicsType(water_air_relation), ForcePriorCK(this->particles_, "BuoyancyForce"),
          gravity_(gravity), threshold_by_dimensions_(threshold * Dimensions),
          dv_mass(this->particles_->template getVariableByName<Real>("Mass")),
          dv_isolation_factor_(nullptr),
          dv_buoyancy_force_(this->getCurrentForce())
    {
        dv_isolation_factor_ = this->contact_particles_->template getVariableByName<Real>("IsolationFactor");
    };
    virtual ~BuoyancyForce() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser),
              gravity_(encloser.gravity_),
              threshold_by_dimensions_(encloser.threshold_by_dimensions_),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              mass_(encloser.dv_mass->DelegatedData(ex_policy)),
              isolation_factor_(encloser.dv_isolation_factor_->DelegatedData(ex_policy)),
              buoyancy_force_(encloser.dv_buoyancy_force_->DelegatedData(ex_policy)) {};

        void interact(size_t index_i, Real dt = 0.0)
        {
            Vecd buoyancy_i = Vecd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                UnsignedInt index_j = this->neighbor_index_[n];
                Real W_ij_Vi = this->W_ij(index_i, index_j) * Vol_[index_i];
                Vecd buoyancy_j = -SMAX(1.0 - isolation_factor_[index_j] / threshold_by_dimensions_, 0.0) *
                                  2.0 * gravity_.InducedAcceleration(Vecd::Zero(), 0.0);
                buoyancy_i += W_ij_Vi * buoyancy_j;
            }
            buoyancy_force_[index_i] = buoyancy_i * mass_[index_i];
        };

      protected:
        Gravity gravity_;
        Real threshold_by_dimensions_;
        Real *Vol_, *mass_, *isolation_factor_;
        Vecd *buoyancy_force_;
    };

  protected:
    const Gravity gravity_;
    Real threshold_by_dimensions_;
    DiscreteVariable<Real> *dv_mass, *dv_isolation_factor_;
    DiscreteVariable<Vecd> *dv_buoyancy_force_;
};
} // namespace fluid_dynamics
} // namespace SPH
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    applyResolutionFromEnv();

    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    GeometricShapeBox water_body_shape(BoundingBoxd(Vecd(0.0, 0.0), Vecd(LL, LH)), "WaterBody");
    GeometricShapeBox inner_wall_shape(BoundingBoxd(Vecd(0.0, 0.0), Vecd(DL, DH)), "InnerWall");
    GeometricShapeBox outer_wall_shape(BoundingBoxd(Vecd(-BW, -BW), Vecd(DL + BW, DH + BW)), "OuterWall");

    ComplexShape air_body_shape("AirBody");
    air_body_shape.add(&inner_wall_shape);
    air_body_shape.subtract(&water_body_shape);

    ComplexShape wall_boundary_shape("WallBoundary");
    wall_boundary_shape.add(&outer_wall_shape);
    wall_boundary_shape.subtract(&inner_wall_shape);

    FluidBody water_body(sph_system, water_body_shape);
    water_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_body.generateParticles<BaseParticles, Lattice>();

    FluidBody air_body(sph_system, air_body_shape);
    air_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_a, c_f);
    air_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, wall_boundary_shape);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody water_observer(sph_system, "FluidObserver");
    water_observer.generateParticles<ObserverParticles>(observation_location);

    auto &water_inner = sph_system.addInnerRelation(water_body);
    auto &water_air_contact = sph_system.addContactRelation(water_body, air_body);
    auto &water_wall_contact = sph_system.addContactRelation(water_body, wall_boundary);
    auto &air_inner = sph_system.addInnerRelation(air_body);
    auto &air_water_contact = sph_system.addContactRelation(air_body, water_body);
    auto &air_wall_contact = sph_system.addContactRelation(air_body, wall_boundary);
    auto &water_observer_contact = sph_system.addContactRelation(water_observer, water_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    auto &host_methods = sph_solver.getHostMethodContainer();

    host_methods.addStateDynamics<NormalFromSubShapeAndOpCK>(wall_boundary, "InnerWall").exec();

    StdVec<FluidBody *> fluid_bodies = {&water_body, &air_body};

    ParticleDynamicsGroup update_fluid_cell_linked_list = main_methods.addCellLinkedListDynamics(fluid_bodies);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    ParticleDynamicsGroup update_fluid_relation;
    update_fluid_relation.add(&main_methods.addRelationDynamics(water_inner, water_air_contact, water_wall_contact));
    update_fluid_relation.add(&main_methods.addRelationDynamics(air_inner, air_water_contact, air_wall_contact));
    ParticleDynamicsGroup update_fluid_configuration = update_fluid_cell_linked_list + update_fluid_relation;
    auto &water_observer_contact_relation = main_methods.addRelationDynamics(water_observer_contact);
    ParticleDynamicsGroup fluid_particle_sorting = main_methods.addSortDynamics(fluid_bodies);

    Gravity gravity(Vecd(0.0, -gravity_g));
    ParticleDynamicsGroup constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(fluid_bodies, gravity);
    ParticleDynamicsGroup advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(fluid_bodies);
    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(fluid_bodies);

    ParticleDynamicsGroup fluid_linear_correction_matrix;
    fluid_linear_correction_matrix.add(
        &main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(water_inner)
             .addPostContactInteraction(water_air_contact)
             .addPostContactInteraction(water_wall_contact));
    fluid_linear_correction_matrix.add(
        &main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(air_inner)
             .addPostContactInteraction(air_water_contact)
             .addPostContactInteraction(air_wall_contact));

    ParticleDynamicsGroup acoustic_step_1st_half;
    acoustic_step_1st_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact));
    acoustic_step_1st_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(air_wall_contact));

    ParticleDynamicsGroup acoustic_step_2nd_half;
    acoustic_step_2nd_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact));
    acoustic_step_2nd_half.add(
        &main_methods.addInteractionDynamicsOneLevel<
                         fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<AcousticRiemannSolverCK, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(air_wall_contact));

    ParticleDynamicsGroup density_regularization;
    density_regularization.add(
        &main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(water_inner)
             .addPostContactInteraction(water_wall_contact)
             .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeSurface>(water_body));
    density_regularization.add(
        &main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(air_inner)
             .addPostContactInteraction(air_water_contact)
             .addPostContactInteraction(air_wall_contact)
             .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, Internal>(air_body));

    ParticleDynamicsGroup transport_correction;
    transport_correction.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, LinearCorrectionCK>(air_inner)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(air_water_contact)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(air_wall_contact)
             .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter>(air_body));
    transport_correction.add(
        &main_methods.addInteractionDynamics<KernelGradientIntegral, LinearCorrectionCK>(water_inner)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(water_air_contact)
             .addPostContactInteraction<Boundary, LinearCorrectionCK>(water_wall_contact)
             .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter>(water_body));

    auto &air_near_wall_bounding = main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::BoundingFromWallCK>(air_wall_contact);

    ParticleDynamicsGroup buoyancy_force;
    buoyancy_force.add(&main_methods.addInteractionDynamics<fluid_dynamics::IsolationFactor>(air_inner));
    buoyancy_force.add(&main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::BuoyancyForce>(water_air_contact, gravity, 0.3));

    ReduceDynamicsGroup advection_time_step =
        main_methods.addReduceDynamics<ReduceMin<Real>, fluid_dynamics::AdvectionTimeStepCK>(fluid_bodies, U_ref);
    ReduceDynamicsGroup acoustic_time_step =
        main_methods.addReduceDynamics<ReduceMin<Real>, fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(fluid_bodies);

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(air_body, "NormalDirection");
    body_state_recorder.addToWrite<Real>(air_body, "SignedDistance");
    body_state_recorder.addToWrite<Real>(air_body, "IsolationFactor");
    body_state_recorder.addToWrite<Vecd>(water_body, "BuoyancyForce");

    auto &observe_pressure = main_methods.addObserveRecorder<Real>(water_observer_contact, "Pressure");
    auto &water_mechanical_energy = main_methods.addReduceDynamics<TotalMechanicalEnergyCK>(water_body, gravity);

    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    if (const char *end_env = std::getenv("GPU_PAPER_END_TIME"))
        total_physical_time = Real(std::atof(end_env));

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "A4_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "validate";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/A4_two_phase_dambreak/schedules/"
              "A4_validate_T10_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "A4", mode);
    recorder.setMeta("case_name", "2d_gpu_paper_a4_two_phase_dambreak");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("particle_spacing_ref", std::to_string(particle_spacing_ref));
    recorder.setMeta("LH_over_dp", std::to_string(lh_over_dp));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("end_time", std::to_string(total_physical_time));
    recorder.setMeta("DL", std::to_string(DL));
    recorder.setMeta("DH", std::to_string(DH));
    recorder.setMeta("LL", std::to_string(LL));
    recorder.setMeta("LH", std::to_string(LH));
    recorder.setMeta("rho0_f", std::to_string(rho0_f));
    recorder.setMeta("rho0_a", std::to_string(rho0_a));
    recorder.setMeta("gravity_g", std::to_string(gravity_g));
    recorder.setMeta("observer_xy", std::to_string(DL) + "," + std::to_string(0.2));
    recorder.setMeta("geometry_note", "CK_source_DL5.3_DH2_probe0.2_not_paper_matched");
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "no_explicit_transfer_API_in_this_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected",
                     "detail deferred to profile short window policy");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected",
                                         "only_on_B2_B3_profile_window");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "needs_valid_telemetry_coverage");
    recorder.registerMissing({"paper_geometry_match", "not_applicable",
                              "source DL=5.3H DH=2H probe y=0.2H; paper ~5.366H/2.5H/0.19H",
                              "case_a4", "freeze matched config before quantitative literature error"});
    recorder.registerMissing({"h1_h2_paper_gauge_locations", "not_collected",
                              "recording strip-max heights at x=LL/2 and x=DL/2 as provisional gauges",
                              "case_a4", "replace with paper gauge x if digitization defines them"});
    recorder.registerMissing({"air_mechanical_energy", "not_collected",
                              "only water TotalMechanicalEnergyCK recorded",
                              "case_a4", "add air energy reduce if required"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        for (int j = 0; j < 100; ++j)
            probe_times.push_back(double(j) * double(total_physical_time) / 99.0);
        recorder.registerMissing({"probe_schedule_json_unreadable", "not_collected",
                                  "failed to parse schedule; fell back to uniform100 of end_time",
                                  schedule_path, "fix GPU_PAPER_SCHEDULE_JSON path"});
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    // Provisional gauges: vertical strip half-width ~ 2 dp
    const Real strip_half = Real(2.0) * particle_spacing_ref;
    const Real h1_x = LL * Real(0.5);
    const Real h2_x = DL * Real(0.5);

    auto host_phase_stats = [&](FluidBody &body, Real rho0, Real gauge_x,
                                Real &mass_sum, Real &x_front, Real &h_gauge, Real &h_max,
                                Real &vmax, Real &rho_dev_max)
    {
        mass_sum = 0.0;
        x_front = -1.0e30;
        h_gauge = 0.0;
        h_max = 0.0;
        vmax = 0.0;
        rho_dev_max = 0.0;
        auto &bp = body.getBaseParticles();
        size_t n = bp.TotalRealParticles();
        auto pos = bp.getVariableDataByName<Vecd>("Position");
        auto vel = bp.getVariableDataByName<Vecd>("Velocity");
        auto mass = bp.getVariableDataByName<Real>("Mass");
        auto rho = bp.getVariableDataByName<Real>("Density");
        for (size_t i = 0; i < n; ++i)
        {
            mass_sum += mass[i];
            x_front = std::max(x_front, Real(pos[i][0]));
            h_max = std::max(h_max, Real(pos[i][1]));
            vmax = std::max(vmax, Real(vel[i].norm()));
            rho_dev_max = std::max(rho_dev_max, std::abs(Real(rho[i] / rho0) - Real(1.0)));
            if (std::abs(Real(pos[i][0]) - gauge_x) <= strip_half)
                h_gauge = std::max(h_gauge, Real(pos[i][1]));
        }
        if (n == 0)
            x_front = 0.0;
    };

    auto sample_probes = [&](double t, double requested)
    {
        water_observer_contact_relation.exec();
        observe_pressure.writeToFile(0);
        auto p_obs = water_observer.getBaseParticles().getVariableDataByName<Real>("Pressure");
        Real pressure = p_obs[0];
        Real mech_e = water_mechanical_energy.exec();

        Real water_mass = 0, air_mass = 0;
        Real x_front = 0, h1 = 0, h2 = 0, h_max = 0, vmax = 0, rho_dev = 0;
        Real air_x = 0, air_h = 0, air_hmax = 0, air_vmax = 0, air_rho_dev = 0;
        host_phase_stats(water_body, rho0_f, h1_x, water_mass, x_front, h1, h_max, vmax, rho_dev);
        // second gauge: recompute strip at h2_x only for water height
        {
            auto &bp = water_body.getBaseParticles();
            size_t n = bp.TotalRealParticles();
            auto pos = bp.getVariableDataByName<Vecd>("Position");
            for (size_t i = 0; i < n; ++i)
                if (std::abs(Real(pos[i][0]) - h2_x) <= strip_half)
                    h2 = std::max(h2, Real(pos[i][1]));
        }
        host_phase_stats(air_body, rho0_a, h2_x, air_mass, air_x, air_h, air_hmax, air_vmax, air_rho_dev);

        auto emit = [&](const std::string &pid, const std::string &body, const std::string &qid,
                        Real value, const std::string &unit, Real x, Real y)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = pid;
            ps.body = body;
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = unit;
            ps.x = double(x);
            ps.y = double(y);
            recorder.addProbe(ps);
        };
        emit("wall_pressure", "FluidObserver", "pressure", pressure, "Pa_code", DL, 0.2);
        emit("water", "WaterBody", "mechanical_energy", mech_e, "J_code", 0.0, 0.0);
        emit("water", "WaterBody", "mass", water_mass, "kg_code", 0.0, 0.0);
        emit("air", "AirBody", "mass", air_mass, "kg_code", 0.0, 0.0);
        emit("water", "WaterBody", "x_front", x_front, "m", x_front, 0.0);
        emit("h1_strip", "WaterBody", "h1", h1, "m", h1_x, h1);
        emit("h2_strip", "WaterBody", "h2", h2, "m", h2_x, h2);
        emit("water", "WaterBody", "h_max", h_max, "m", 0.0, h_max);
        emit("water", "WaterBody", "vmax", vmax, "m/s", 0.0, 0.0);
        emit("water", "WaterBody", "rho_dev_max", rho_dev, "1", 0.0, 0.0);
        emit("air", "AirBody", "vmax", air_vmax, "m/s", 0.0, 0.0);
        emit("air", "AirBody", "rho_dev_max", air_rho_dev, "1", 0.0, 0.0);

        recorder.addDiagnostic(t, "WaterBody", "N_active",
                               double(water_body.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "AirBody", "N_active",
                               double(air_body.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "WaterBody", "pressure_observer", double(pressure), "Pa_code", "measured");
        recorder.addDiagnostic(t, "WaterBody", "x_front", double(x_front), "m", "measured");
        recorder.addDiagnostic(t, "WaterBody", "mechanical_energy", double(mech_e), "J_code", "measured");
    };

    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep();
    int screening_interval = 100;

    wall_cell_linked_list.exec();
    update_fluid_configuration.exec();
    water_observer_contact_relation.exec();

    constant_gravity.exec();
    density_regularization.exec();
    advection_step_setup.exec();
    fluid_linear_correction_matrix.exec();
    buoyancy_force.exec();

    auto wall0 = std::chrono::steady_clock::now();
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;

    if (enable_vtp)
        body_state_recorder.writeToFile();
    {
        double req0 = probe_times.empty() ? 0.0 : probe_times[0];
        sample_probes(time_stepper.getPhysicalTime(), req0);
        if (!probe_times.empty())
            probe_cursor = 1;
        if (enable_vtp && !field_times.empty() && field_times[0] <= 0.0)
            field_cursor = 1;
    }

    TickCount t0 = TickCount::now();
    int dt_fail = 0;
    size_t acoustic_steps = 0;

    while (!time_stepper.isEndTime(total_physical_time))
    {
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        if (!(acoustic_dt > Real(0.0)) || !std::isfinite(double(acoustic_dt)))
        {
            std::cout << "ERROR: non-positive/non-finite acoustic_dt=" << acoustic_dt
                      << " at advection_N=" << advection_steps
                      << " Time=" << time_stepper.getPhysicalTime() << std::endl;
            dt_fail = 1;
            break;
        }
        acoustic_step_1st_half.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        acoustic_steps++;

        time_instance = TickCount::now();
        double t_now = time_stepper.getPhysicalTime();
        double requested = 0.0;
        if (dueForSample(t_now, probe_times, probe_cursor, &requested))
            sample_probes(t_now, requested);
        if (enable_vtp && dueForSample(t_now, field_times, field_cursor, &requested))
            body_state_recorder.writeToFile();
        interval_output += TickCount::now() - time_instance;

        if (advection_step(advection_time_step))
        {
            advection_steps++;
            update_particle_position.exec();
            air_near_wall_bounding.exec();

            time_instance = TickCount::now();
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }
            interval_output += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            if (advection_steps % 100)
            {
                fluid_particle_sorting.exec();
            }
            update_fluid_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            density_regularization.exec();
            advection_step_setup.exec();
            transport_correction.exec();
            fluid_linear_correction_matrix.exec();
            buoyancy_force.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }

    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    // Final sample at end
    sample_probes(time_stepper.getPhysicalTime(), double(total_physical_time));

    recorder.stages().add("acoustic_step", interval_acoustic_step.seconds());
    recorder.stages().add("advection_physics", interval_advection_step.seconds());
    recorder.stages().add("configuration", interval_updating_configuration.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("acoustic", int(acoustic_steps), 0, 0, 0, 0, 0, 0, "fluid_ck", "measured");
    recorder.addStepBlock("advection", int(advection_steps), 0, 0, 0, 0, 0, 0, "fluid_ck", "measured");
    recorder.addBodyStats(time_stepper.getPhysicalTime(), "WaterBody",
                          int(water_body.getBaseParticles().TotalRealParticles()),
                          int(water_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");
    recorder.addBodyStats(time_stepper.getPhysicalTime(), "AirBody",
                          int(air_body.getBaseParticles().TotalRealParticles()),
                          int(air_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    const char *status = dt_fail ? "timestep_collapse" : "ok";
    recorder.finalize(status, 0.0, time_stepper.getPhysicalTime(), wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval_output.seconds(),
                      0.0, interval_output.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return dt_fail;
}
