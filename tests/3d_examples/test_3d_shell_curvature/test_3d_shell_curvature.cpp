/**
 * @file 	channel_flow_shell.cpp
 * @brief 	This is a test of fluid-shell interaction.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real radius = 1.0;
const Real resolution_ref = 0.1;
const Real shell_thickness = resolution_ref;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius - resolution_ref, -radius - resolution_ref, -radius - resolution_ref), Vec3d(radius + resolution_ref, radius + resolution_ref, radius + resolution_ref));
//----------------------------------------------------------------------
//	Material Parameters
//----------------------------------------------------------------------
Real rho0_s = 10.0;
Real Youngs_modulus = 1e3;
Real poisson = 0.3;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ParticleGeneratorDirect final : public SPH::ParticleGenerator<Base>
{
    // Particle generation fully performed by SPHBody::generateParticles(...) in base_body.h
    // Enable construction only by a SPHBody to prevent misuse
    friend class SPH::SPHBody;

    const SPH::StdLargeVec<SPH::Vecd> &positions;
    const SPH::StdLargeVec<SPH::Real> &volumes;

    ParticleGeneratorDirect(
        SPH::SPHBody &body,
        const SPH::StdLargeVec<SPH::Vecd> &positions,
        const SPH::StdLargeVec<SPH::Real> &volumes) : ParticleGenerator(body),
                                                      positions(positions),
                                                      volumes(volumes)
    {
    }

    void initializeGeometricVariables() override
    {
        for (size_t i = 0; i < positions.size(); ++i)
        {
            initializePositionAndVolumetricMeasure(positions[i], volumes[i]);
        }
    };
};
class ShellParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real dphi = resolution_ref / radius;
        int particle_number_phi = int(2 * Pi / dphi);
        for (int i = 0; i < particle_number_phi; i++)
        {
            Real phi = dphi * i;
            Real radius_phi = radius * sin(phi);
            Real z = radius * cos(phi);
            int particle_number_theta = int(2 * Pi * radius_phi / resolution_ref);
            for (int j = 0; j < particle_number_theta; j++)
            {
                Real theta = (j + 0.5) * 2 * Pi / (Real)particle_number_theta;
                Real x = radius_phi * cos(theta);
                Real y = radius_phi * sin(theta);
                initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
                                                       resolution_ref * resolution_ref);
                Vec3d n_0 = Vec3d(x / radius, y / radius, z / radius);
                initializeSurfaceProperties(n_0, shell_thickness);
            }
        }
    }
};
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ControlDisplacement : public LocalDynamics, public thin_structure_dynamics::ShellDataSimple
{
  public:
    ControlDisplacement(SPHBody &sph_body, Real moving_v)
        : LocalDynamics(sph_body),
          thin_structure_dynamics::ShellDataSimple(sph_body),
          n_(particles_->n_), vel_(particles_->vel_), moving_v_(moving_v){};

  protected:
    StdLargeVec<Vecd> &n_;
    StdLargeVec<Vecd> &vel_;
    Real moving_v_;

    inline void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = moving_v_ * n_[index_i];
    };
};

int main(int ac, char *av[])
{
    const Real end_time = 10.0;
    const Real target_radius_increase = 0.5 * radius;
    const Real moving_v = target_radius_increase / end_time;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //---------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson); // dummy material parameters
    shell.generateParticles<ShellParticleGenerator>();

    StdLargeVec<Vecd> water_positions({Vec3d(0, 0, 0)});
    StdLargeVec<Real> water_volumes({resolution_ref * resolution_ref * resolution_ref});
    FluidBody dummy_fluid(sph_system, makeShared<DefaultShape>("DummyFluid"));
    dummy_fluid.defineAdaptation<SPHAdaptation>(1.3, radius / (radius + target_radius_increase));
    dummy_fluid.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1.0, 0, 0);
    dummy_fluid.generateParticles<ParticleGeneratorDirect>(water_positions, water_volumes);
    //----------------------------------------------------------------------
    // Contact
    //----------------------------------------------------------------------
    InnerRelation shell_inner_contact(shell);
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, dummy_fluid);
    //----------------------------------------------------------------------
    //  Solid Algorithms
    //----------------------------------------------------------------------
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
        corrected_configuration(shell_inner_contact);
    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell);
    /** stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_inner_contact);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_inner_contact);
    /** Control the displacement. */
    SimpleDynamics<ControlDisplacement> dis_control(shell, moving_v);
    /** Constrain the Boundary. */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
        cylinder_position_damping(0.2, shell_inner_contact, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
        cylinder_rotation_damping(0.2, shell_inner_contact, "AngularVelocity", physical_viscosity);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
    auto update_shell_volume = [&]()
    {
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                shell.getBaseParticles().Vol_[index_i] = shell.getBaseParticles().mass_[index_i] / shell.getBaseParticles().rho_[index_i];
            });
    };
    /** Compute curvature. */
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_initial_curvature(shell_inner_contact);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    // compute initial curvature after B_ is computed
    shell_initial_curvature.compute_initial_curvature();
    shell_average_curvature.exec();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("Average2ndPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("2ndPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    write_real_body_states.writeToFile();
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;

    /** Setup time stepping control parameters. */
    int ite = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_period)
        {
            shell_update_normal.exec();
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            dis_control.exec(dt);
            cylinder_position_damping.exec(dt);
            cylinder_rotation_damping.exec(dt);
            dis_control.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        TickCount t2 = TickCount::now();
        shell_initial_curvature.exec();
        update_shell_volume();
        shell.updateCellLinkedList();
        shell_curvature_inner.updateConfiguration();
        shell_average_curvature.exec();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}