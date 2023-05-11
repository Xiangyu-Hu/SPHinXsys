/**
 * @file 	poiseuille_flow.cpp
 * @brief 	3D poiseuille flow example
 * @details This is the one of the basic test cases for validating viscous flow,
 with using emitter as inlet and disposer as outlet.
                        // TODO: include shell in next patch.
 * @author
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "base_data_type.h"
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

class ObersverAxial : public ObserverParticleGenerator {
public:
  ObersverAxial(SPHBody &sph_body, double full_length,
                Vec3d translation = Vec3d(0.0, 0.0, 0.0))
      : ObserverParticleGenerator(sph_body) {
    int ny = 51;
    for (int i = 0; i < ny; i++) {
      double y = full_length / (ny - 1) * i;
      Vec3d point_coordinate(0.0, y, 0.0);
      positions_.emplace_back(point_coordinate + translation);
    }
  }
};

class ObserverRadial : public ObserverParticleGenerator {
public:
  ObserverRadial(SPHBody &sph_body, double full_length, double diameter,
                 int number_of_particles,
                 Vec3d translation = Vec3d(0.0, 0.0, 0.0))
      : ObserverParticleGenerator(sph_body) {

    int n = number_of_particles + 1;
    double y = full_length / 2.0;
    for (int i = 0; i < n - 1;
         i++) // we leave out the point clsoe to the boundary as the
              // interpolation there is incorrect
    {
      double z = diameter / 2.0 * i / double(n);
      positions_.emplace_back(Vec3d(0.0, y, z) + translation);
      positions_.emplace_back(Vec3d(0.0, y, -z) + translation);
    }
  }
};

/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */

const Real scale = 0.001;
const Real diameter = 6.35 * scale * 1;
const Real fluid_radius = 0.5 * diameter;
const Real full_length = 15 * fluid_radius;
const int number_of_particles = 10;
const Real resolution_ref = diameter / number_of_particles;
const Real inflow_length = resolution_ref * 20.0; // Inflow region
const Real wall_thickness = resolution_ref * 4.0;
const Real shell_thickness = resolution_ref * 2.0;
const int simtk_resolution = 20;
const Vec3d translation_fluid(0., full_length * 0.5, 0.);
/**
 * @brief Geometry parameters for shell.
 */
const Real radius_mid_surface = fluid_radius + shell_thickness * 0.5;
const int particle_number_mid_surface =
    int(2.0 * radius_mid_surface * Pi / resolution_ref);
const int particle_number_height =
    2 * int((full_length * 0.5 + wall_thickness) / resolution_ref);
const int BWD =
    1; /** Width of the boundary layer measured by number of particles. */

/**
 * @brief Geometry parameters for boundary condition.
 */
const Vec3d emitter_halfsize =
    Vec3d(fluid_radius, resolution_ref * 2, fluid_radius);
const Vec3d emitter_translation = Vec3d(0., resolution_ref * 2, 0.);
const Vec3d emitter_buffer_halfsize =
    Vec3d(fluid_radius, inflow_length * 0.5, fluid_radius);
const Vec3d emitter_buffer_translation = Vec3d(0., inflow_length * 0.5, 0.);
const Vec3d disposer_halfsize =
    Vec3d(fluid_radius * 2, resolution_ref * 2, fluid_radius * 2);
const Vec3d disposer_translation =
    Vec3d(0., full_length, 0.) - Vec3d(0., disposer_halfsize[1], 0.);

/** Domain bounds of the system. */
const BoundingBox system_domain_bounds(Vec3d(-diameter, 0, -diameter) -
                                           Vec3d(wall_thickness, wall_thickness,
                                                 wall_thickness),
                                       Vec3d(diameter, full_length, diameter) +
                                           Vec3d(wall_thickness, wall_thickness,
                                                 wall_thickness));
/**
 * @brief Material properties of the fluid.
 */
const Real rho0_f = 1000.0; /**< Reference density of fluid. */
const Real mu_f = 6.5e-3;   /**< Viscosity. */
const Real Re = 10;
// const Real U_f = 30.0e-6 / 60. / (M_PI / 4.0 * diameter * diameter); /**<
// Characteristic velocity. Average velocity */
const Real U_f = Re / diameter * mu_f;
const Real U_max = 2.0 * U_f;  // parabolic inflow, Thus U_max = 2*U_f
const Real c_f = 10.0 * U_max; /**< Reference sound speed. */
/**
 * @brief Define water shape
 */
class WaterBlock : public ComplexShape {
public:
  explicit WaterBlock(const std::string &shape_name)
      : ComplexShape(shape_name) {
    add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius,
                                   full_length * 0.5, simtk_resolution,
                                   translation_fluid);
  }
};
/**
 * @brief Define wall shape
 */
class ShellBoundary : public SurfaceParticleGenerator {
public:
  explicit ShellBoundary(SPHBody &sph_body)
      : SurfaceParticleGenerator(sph_body){};
  void initializeGeometricVariables() override {
    for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++) {
      for (int j = 0; j < particle_number_height; j++) {
        Real x = radius_mid_surface *
                 cos(Pi + (i - BWD + 0.5) * 2 * Pi /
                              (Real)particle_number_mid_surface);
        Real y = (j - particle_number_height / 2) * resolution_ref +
                 resolution_ref * 0.5 + full_length * 0.5;
        Real z = radius_mid_surface *
                 sin(Pi + (i - BWD + 0.5) * 2 * Pi /
                              (Real)particle_number_mid_surface);
        initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
                                               resolution_ref * resolution_ref);
        Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
        initializeSurfaceProperties(n_0, shell_thickness);
      }
    }
  }
};
/**
 * @brief Inflow velocity
 */
struct InflowVelocity {
  Real u_ref_, t_ref_;
  AlignedBoxShape &aligned_box_;
  Vec3d halfsize_;

  template <class BoundaryConditionType>
  InflowVelocity(BoundaryConditionType &boundary_condition)
      : u_ref_(U_f), t_ref_(2.0),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {}

  Vec3d operator()(Vec3d &position, Vec3d &velocity) {
    Vec3d target_velocity = Vec3d(0, 0, 0);
    if (aligned_box_.checkInBounds(0, position)) {
      target_velocity[1] =
          2.0 * U_f *
          (1.0 - (position[0] * position[0] + position[2] * position[2]) /
                     fluid_radius / fluid_radius);
    }
    return target_velocity;
  }
};

/**
 * @brief 	Main program starts here.
 */
int main() {
  /**
   * @brief Build up -- a SPHSystem --
   */
  SPHSystem system(system_domain_bounds, resolution_ref);
  /** Set the starting time. */
  GlobalStaticVariables::physical_time_ = 0.0;
  /**
   * @brief Material property, particles and body creation of fluid.
   */
  FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
  water_block
      .defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(
          rho0_f, c_f, mu_f);
  water_block.generateParticles<ParticleGeneratorLattice>();
  /**
   * @brief 	Particle and body creation of wall boundary.
   */
  SolidBody shell_boundary(system, makeShared<DefaultShape>("Shell"));
  shell_boundary.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(
      rho0_f, 1e3, 0.45);
  shell_boundary.generateParticles<ShellBoundary>();
  /** topology */
  InnerRelation water_inner(water_block);
  ComplexRelation water_block_complex(water_block, {&shell_boundary});
  /**
   * @brief 	Define all numerical methods which are used in this case.
   */
  /**
   * @brief 	Methods used for time stepping.
   */
  /** Define external force. */
  SimpleDynamics<NormalDirectionFromBodyShape> shell_boundary_normal_direction(
      shell_boundary);
  /** Initialize particle acceleration. */
  SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(
      water_block, makeShared<Gravity>(Vec3d(0.0, 0.0, 0.0)));
  /**
   * @brief 	Algorithms of fluid dynamics.
   */
  /** time-space method to detect surface particles. (Important for
   * DensitySummationFreeSurfaceComplex work correctly.)*/
  InteractionWithUpdate<
      fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
      inlet_outlet_surface_particle_indicator(water_block_complex);
  /** Evaluation of density by summation approach. */
  InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex>
      update_density_by_summation(water_block_complex);
  /** Time step size without considering sound wave speed. */
  ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize>
      get_fluid_advection_time_step_size(water_block, U_max);
  /** Time step size with considering sound wave speed. */
  ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(
      water_block);
  /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
  Dynamics1Level<fluid_dynamics::FluidShellIntegration1stHalfRiemann>
      pressure_relaxation(water_block_complex);
  /** Pressure relaxation algorithm by using position verlet time stepping. */
  Dynamics1Level<fluid_dynamics::FluidShellIntegration2ndHalfRiemann>
      density_relaxation(water_block_complex);
  /** Computing viscous acceleration. */
  InteractionDynamics<fluid_dynamics::ViscousAccelerationWithShell>
      viscous_acceleration(water_block_complex);
  /** Impose transport velocity. */
  InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
      transport_velocity_correction(water_block_complex);
  /**
   * @brief 	Boundary conditions. Inflow & Outflow in Y-direction
   */
  BodyAlignedBoxByParticle emitter(
      water_block,
      makeShared<AlignedBoxShape>(Transform3d(Vec3d(emitter_translation)),
                                  emitter_halfsize));
  SimpleDynamics<fluid_dynamics::EmitterInflowInjection,
                 BodyAlignedBoxByParticle>
      emitter_inflow_injection(emitter, 10, 1);
  /** Emitter buffer inflow condition. */
  BodyAlignedBoxByCell emitter_buffer(
      water_block, makeShared<AlignedBoxShape>(
                       Transform3d(Vec3d(emitter_buffer_translation)),
                       emitter_buffer_halfsize));
  SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>,
                 BodyAlignedBoxByCell>
      emitter_buffer_inflow_condition(emitter_buffer);
  BodyAlignedBoxByCell disposer(
      water_block,
      makeShared<AlignedBoxShape>(Transform3d(Vec3d(disposer_translation)),
                                  disposer_halfsize));
  SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion, BodyAlignedBoxByCell>
      disposer_outflow_deletion(disposer, 1);

  /**
   * @brief Output.
   */
  IOEnvironment io_environment(system);
  water_block.addBodyStateForRecording<int>("SurfaceIndicator");

  shell_boundary.addBodyStateForRecording<double>("MassiveMeasure");
  shell_boundary.addBodyStateForRecording<double>("Density");
  shell_boundary.addBodyStateForRecording<double>("VolumetricMeasure");
  shell_boundary.addBodyStateForRecording<double>("Thickness");
  /** Output the body states. */
  BodyStatesRecordingToVtp body_states_recording(io_environment,
                                                 system.real_bodies_);
  /**
   * @brief OBSERVER.
   */
  ObserverBody observer_axial(system, "fluid_observer_axial");
  observer_axial.generateParticles<ObersverAxial>(full_length);
  ObserverBody observer_radial(system, "fluid_observer_radial");
  observer_radial.generateParticles<ObserverRadial>(full_length, diameter,
                                                    number_of_particles);
  ContactRelation observer_contact_axial(observer_axial, {&water_block});
  ContactRelation observer_contact_radial(observer_radial, {&water_block});
  ObservedQuantityRecording<Vec3d> write_fluid_velocity_axial(
      "Velocity", io_environment, observer_contact_axial);
  ObservedQuantityRecording<Vec3d> write_fluid_velocity_radial(
      "Velocity", io_environment, observer_contact_radial);
  /**
   * @brief Setup geometry and initial conditions.
   */
  system.initializeSystemCellLinkedLists();
  system.initializeSystemConfigurations();
  shell_boundary_normal_direction.parallel_exec();
  /** Output the start states of bodies. */
  body_states_recording.writeToFile(0);
  /**
   * @brief 	Basic parameters.
   */
  size_t number_of_iterations = system.RestartStep();
  int screen_output_interval = 100;
  Real end_time = 10.0;   /**< End time. */
  Real Output_Time = 0.1; /**< Time stamps for output of body states. */
  Real dt = 0.0;          /**< Default acoustic time step sizes. */
  /** statistics for computing CPU time. */
  tick_count t1 = tick_count::now();
  tick_count::interval_t interval;
  tick_count::interval_t interval_computing_time_step;
  tick_count::interval_t interval_computing_pressure_relaxation;
  tick_count::interval_t interval_updating_configuration;
  tick_count time_instance;
  /**
   * @brief 	Main loop starts here.
   */
  while (GlobalStaticVariables::physical_time_ < end_time) {
    Real integration_time = 0.0;
    /** Integrate time (loop) until the next output time. */
    while (integration_time < Output_Time) {
      /** Acceleration due to viscous force and gravity. */
      time_instance = tick_count::now();
      initialize_a_fluid_step.parallel_exec();
      Real Dt = get_fluid_advection_time_step_size.parallel_exec();
      inlet_outlet_surface_particle_indicator.parallel_exec();
      update_density_by_summation.parallel_exec();
      viscous_acceleration.parallel_exec();
      transport_velocity_correction.parallel_exec();
      interval_computing_time_step += tick_count::now() - time_instance;
      /** Dynamics including pressure relaxation. */
      time_instance = tick_count::now();
      Real relaxation_time = 0.0;
      while (relaxation_time < Dt) {
        dt = SMIN(get_fluid_time_step_size.parallel_exec(),
                  Dt - relaxation_time);
        pressure_relaxation.parallel_exec(dt);
        emitter_buffer_inflow_condition.parallel_exec();
        density_relaxation.parallel_exec(dt);
        relaxation_time += dt;
        integration_time += dt;
        GlobalStaticVariables::physical_time_ += dt;
      }
      interval_computing_pressure_relaxation +=
          tick_count::now() - time_instance;
      if (number_of_iterations % screen_output_interval == 0) {
        std::cout << std::fixed << std::setprecision(9)
                  << "N=" << number_of_iterations
                  << "	Time = " << GlobalStaticVariables::physical_time_
                  << "	Dt = " << Dt << "	dt = " << dt << "\n";
        body_states_recording.writeToFile();
      }
      number_of_iterations++;
      /** Update cell linked list and configuration. */
      time_instance = tick_count::now();

      /** Water block configuration and periodic condition. */
      emitter_inflow_injection.parallel_exec();
      disposer_outflow_deletion.parallel_exec();

      water_block.updateCellLinkedListWithParticleSort(100);
      water_block_complex.updateConfiguration();
      interval_updating_configuration += tick_count::now() - time_instance;
    }
    tick_count t2 = tick_count::now();
    tick_count t3 = tick_count::now();
    interval += t3 - t2;

    /** Update observer and write output of observer. */
    observer_contact_axial.updateConfiguration();
    observer_contact_radial.updateConfiguration();
    write_fluid_velocity_axial.writeToFile(number_of_iterations);
    write_fluid_velocity_radial.writeToFile(number_of_iterations);
  }
  tick_count t4 = tick_count::now();

  tick_count::interval_t tt;
  tt = t4 - t1 - interval;
  std::cout << "Total wall time for computation: " << tt.seconds()
            << " seconds." << std::endl;
  std::cout << std::fixed << std::setprecision(9)
            << "interval_computing_time_step ="
            << interval_computing_time_step.seconds() << "\n";
  std::cout << std::fixed << std::setprecision(9)
            << "interval_computing_pressure_relaxation = "
            << interval_computing_pressure_relaxation.seconds() << "\n";
  std::cout << std::fixed << std::setprecision(9)
            << "interval_updating_configuration = "
            << interval_updating_configuration.seconds() << "\n";

  /**
   * @brief 	Gtest strat from here.
   */
  /* Define analytical solution of the inflow velocity.*/
  std::function<Vec3d(Vec3d)> inflow_velocity = [&](Vec3d pos) {
    return Vec3d(0.0,
                 2.0 * U_f *
                     (1.0 - (pos[0] * pos[0] + pos[2] * pos[2]) /
                                (diameter * 0.5) / (diameter * 0.5)),
                 0.0);
  };
  /* Compare all simulation to the anayltical solution. */
  // Axial direction.
  for (size_t i = 0; i < observer_axial.getBaseParticles().pos_.size(); i++) {
    EXPECT_NEAR(inflow_velocity(observer_axial.getBaseParticles().pos_[i])[1],
                observer_axial.getBaseParticles().vel_[i][1],
                U_max * 10e-2); // it's below 5% but 10% for CI
  }
  // Radial direction
  for (size_t i = 0; i < observer_radial.getBaseParticles().pos_.size(); i++) {
    EXPECT_NEAR(inflow_velocity(observer_radial.getBaseParticles().pos_[i])[1],
                observer_radial.getBaseParticles().vel_[i][1],
                U_max * 10e-2); // it's below 5% but 10% for CI
  }

  return 0;
}
