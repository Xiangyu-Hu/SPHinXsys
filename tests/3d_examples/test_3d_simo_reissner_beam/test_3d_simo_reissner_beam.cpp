#include "beam_wrapper.h"
#include <gtest/gtest.h>

void arc_steady_test(int res_factor);
void elbow_dynamic_test(int res_factor);

TEST(SimoReisnnerBeam, steady)
{
    arc_steady_test(4);
}

TEST(SimoReisnnerBeam, dynamic)
{
    elbow_dynamic_test(4);
}

int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    testing::GTEST_FLAG(filter) = "SimoReisnnerBeam.dynamic";
    return RUN_ALL_TESTS();
}

void arc_steady_test(int res_factor)
{
    // Global parameters
    const Real Radius = 100;
    const Real arc_angle = Pi / 4.0; // 45 degrees
    const Real length = arc_angle * Radius;
    const Real width = 1.0;
    const Real thickness = 1.0;

    // resolutions
    const Real dtheta = arc_angle / (10.0 * res_factor); // angle between particles
    const Real dp = length / (10.0 * res_factor);        // distance between particles

    BarParameters bar_params;
    bar_params.body_name = "Beam";
    bar_params.hourglass_control = false;

    // bar geometry parameters
    {
        bar_params.geometry_params.dp_ = dp;
        Real theta = 0.5 * dtheta;
        while (theta < arc_angle)
        {
            bar_params.geometry_params.pos_.emplace_back(Radius * sin(theta), Radius * (1 - cos(theta)), 0.0);
            bar_params.geometry_params.n0_.emplace_back(Vec3d::UnitZ());
            bar_params.geometry_params.b_n0_.emplace_back(-sin(theta), cos(theta), 0.0);
            bar_params.geometry_params.width_.emplace_back(width);
            bar_params.geometry_params.thickness_.emplace_back(thickness);
            theta += dtheta;
        }
    }
    size_t particle_number = bar_params.geometry_params.pos_.size();

    // material
    const Real rho0_s = 1.0;
    const Real Youngs_modulus = 1e7;
    const Real poisson = 0.0;

    // Material models
    bar_params.material = std::make_shared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    Real G = bar_params.material->ShearModulus();

    // Material properties
    Real area = get_area_rectangular_shape(width, thickness);
    Real I_y = get_moment_of_inertia_rectangular_shape(width, thickness);
    Real I_z = get_moment_of_inertia_rectangular_shape(thickness, width);
    Real J = I_y + I_z; // torsional constant for rectangular section
    Mat3d I_rho = rho0_s * Vec3d(J, I_y, I_z).asDiagonal();
    Mat3d C_N = Vec3d(Youngs_modulus * area, G * area, G * area).asDiagonal();
    Mat3d C_M = Vec3d(G * J, Youngs_modulus * I_y, Youngs_modulus * I_z).asDiagonal();
    bar_params.material_params.rho_A_.resize(particle_number, rho0_s * area);
    bar_params.material_params.I_rho_l_.resize(particle_number, I_rho);
    bar_params.material_params.C_N_.resize(particle_number, C_N);
    bar_params.material_params.C_M_.resize(particle_number, C_M);
    bar_params.physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, thickness);

    // load
    Real force = 600;

    // Setup the system
    bar_simulation simulation(dp);

    // Create a bar body
    simulation.add_bar_object(bar_params);
    auto &bar_object = simulation.objects[0];
    auto &bar_particles = bar_object.bar_body->getBaseParticles();

    // System initialization
    simulation.initialize_system();

    // Boundary conditions
    // clamping condition
    Vec3d fixed_pos = Vec3d::Zero();
    IndexVector fixed_ids = [&]()
    {
        const auto *pos = bar_particles.getVariableDataByName<Vec3d>("Position");
        IndexVector ids;
        for (size_t i = 0; i < bar_particles.TotalRealParticles(); ++i)
        {
            if ((pos[i] - fixed_pos).norm() < 1.3 * dp)
                ids.push_back(i);
        }
        return ids;
    }();
    BodyPartByParticle fixed_part(*bar_object.bar_body);
    fixed_part.body_part_particles_ = fixed_ids;
    SimpleDynamics<slender_structure_dynamics::ConstrainBarBodyRegion> fix_bc(fixed_part);

    // loading condition
    Vec3d tip_pos = Radius * Vec3d(sin(arc_angle), 1 - cos(arc_angle), 0.0);
    auto tip_id = get_closest_particle(bar_particles, tip_pos);
    auto *force_prior = bar_particles.registerStateVariableData<Vec3d>("ForcePrior");
    auto loading_bc = [&]()
    {
        force_prior[tip_id] = force * Vec3d::UnitZ();
    };

    // Output
    bar_particles.addVariableToWrite<Vec3d>("DrDs");
    bar_particles.addVariableToWrite<Vec3d>("Velocity");
    bar_particles.addVariableToWrite<Vec3d>("ForcePrior");
    bar_particles.addVariableToWrite<Vec3d>("PriorAngularAcceleration");
    bar_particles.addVariableToWrite<Vec3d>("Force");
    bar_particles.addVariableToWrite<Vec3d>("AngularVelocity");
    bar_particles.addVariableToWrite<Vec3d>("AngularAcceleration");
    bar_particles.addVariableToWrite<Vec3d>("PseudoNormal");
    bar_particles.addVariableToWrite<Vec3d>("PseudoBinormal");
    BodyStatesRecordingToVtp vtp_output(simulation.system);
    vtp_output.writeToFile(0);

    // Observer
    ObserverBody tip(simulation.system, "tip");
    tip.generateParticles<ObserverParticles>(std::vector{tip_pos});
    ContactRelation observer_relation(tip, {bar_object.bar_body.get()});
    observer_relation.updateConfiguration();
    InteractionDynamics<CorrectInterpolationKernelWeights>{observer_relation}.exec();
    ObservedQuantityRecording<Vec3d> obs_tip_displacement("Position", observer_relation);
    auto *obs_tip_pos = tip.getBaseParticles().getVariableDataByName<Vec3d>("Position");

    simulation.output_function = [&](size_t ite)
    {
        vtp_output.writeToFile(ite + 1);
        obs_tip_displacement.writeToFile(ite + 1);
    };
    simulation.acceleration_bc = [&](Real dt)
    {
        loading_bc();
    };
    simulation.velocity_bc = [&](Real dt)
    {
        fix_bc.exec();
    };
    Real dt_ref = length / (bar_params.material->ReferenceSoundSpeed() / 20.0);
    Real end_time = 50 * dt_ref;
    simulation.output_number = 100;
    try
    {
        simulation.run_until(end_time, true);
    }
    catch (const std::runtime_error &e)
    {
        vtp_output.writeToFile(simulation.output_ite + 1);
        std::cout << "Simulation failed: " << e.what() << std::endl;
    }

    // post-processing and test
    // reference value
    const Vec3d disp_ref = Vec3d(47.23, 15.79, 53.37) - tip_pos;
    const Vec3d disp = obs_tip_pos[0] - tip_pos;
    const Real rel_error = (disp - disp_ref).norm() / disp_ref.norm();
    std::cout << "Tip displacement: " << disp.transpose() << ", relative error: " << rel_error * 100 << "%" << std::endl;
    EXPECT_LE(rel_error, 10e-2); // 10%
}

void elbow_dynamic_test(int res_factor)
{
    // Global parameters
    const Real length = 10.0;
    const Real radius = 0.063;
    const Real end_time = 30.0;

    // resolutions
    const Real dp = length / Real(10 * res_factor);

    // material parameters
    const Real EA = 1e6;
    const Real GA = EA;
    const Real EI = 1e3;
    const Real GJ = EI;
    // The inertia properties are scaled, not the real ones
    const Real rho_J2 = 10;
    const Real rho_J1 = 2 * rho_J2;
    const Real rho_A = 1;

    // bar geometry parameters
    BarParameters bar_params;
    bar_params.body_name = "Elbow";
    bar_params.hourglass_control = true;
    {
        Vec3d g2_1 = Vec3d::UnitY();
        Vec3d g2_2 = -Vec3d::UnitX();
        Vec3d g3 = Vec3d::UnitZ();

        bar_params.geometry_params.dp_ = dp;

        // add the elbow point
        {
            bar_params.geometry_params.pos_.emplace_back(Vec3d::Zero());
            bar_params.geometry_params.n0_.emplace_back(g3);
            bar_params.geometry_params.b_n0_.emplace_back((g2_1 + g2_2).normalized());
            bar_params.geometry_params.width_.emplace_back(radius);
            bar_params.geometry_params.thickness_.emplace_back(radius);
        }

        // add the horizontal part
        Real x = -dp;
        while (x > -length)
        {
            Vec3d position = x * Vec3d::UnitX();
            bar_params.geometry_params.pos_.emplace_back(position);
            bar_params.geometry_params.n0_.emplace_back(g3);
            bar_params.geometry_params.b_n0_.emplace_back(g2_1);
            bar_params.geometry_params.width_.emplace_back(radius);
            bar_params.geometry_params.thickness_.emplace_back(radius);
            x -= dp;
        }

        // add the vertical part
        Real y = dp;
        while (y < length)
        {
            Vec3d position = y * Vec3d::UnitY();
            bar_params.geometry_params.pos_.emplace_back(position);
            bar_params.geometry_params.n0_.emplace_back(g3);
            bar_params.geometry_params.b_n0_.emplace_back(g2_2);
            bar_params.geometry_params.width_.emplace_back(radius);
            bar_params.geometry_params.thickness_.emplace_back(radius);
            y += dp;
        }
    }

    // material
    const Real rho0_s = 83;
    const Real Youngs_modulus = 8.3e7;
    const Real poisson = 0;

    // Material models
    bar_params.material = std::make_shared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);

    // Material properties
    size_t particle_number = bar_params.geometry_params.pos_.size();
    Mat3d I_rho = Vec3d(rho_J1, rho_J2, rho_J2).asDiagonal();
    Mat3d C_N = Vec3d(EA, GA, GA).asDiagonal();
    Mat3d C_M = Vec3d(GJ, EI, EI).asDiagonal();
    bar_params.material_params.rho_A_.resize(particle_number, rho_A);
    bar_params.material_params.I_rho_l_.resize(particle_number, I_rho);
    bar_params.material_params.C_N_.resize(particle_number, C_N);
    bar_params.material_params.C_M_.resize(particle_number, C_M);
    bar_params.physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, radius);

    // load
    const Vec3d force_direction = Vec3d::UnitZ();
    auto force = [F_max = 50.0, t_ref = 1.0](Real time)
    {
        if (time < t_ref)
            return F_max * time / t_ref;
        else if (time < 2.0 * t_ref)
            return F_max * (2.0 - time / t_ref);
        else
            return 0.0;
    };

    // Setup the system
    bar_simulation simulation(dp);

    // Create a bar body
    simulation.add_bar_object(bar_params);
    auto &bar_object = simulation.objects[0];
    auto &bar_particles = bar_object.bar_body->getBaseParticles();
    // reset mass
    {
        auto *mass = bar_particles.getVariableDataByName<Real>("Mass");
        for (size_t i = 0; i < bar_particles.TotalRealParticles(); ++i)
        {
            mass[i] = rho_A * dp;
        }
    }

    // System initialization
    simulation.initialize_system();

    // Boundary conditions
    // clamping condition
    Vec3d fixed_pos = -length * Vec3d::UnitX();
    IndexVector fixed_ids = [&]()
    {
        const auto *pos = bar_particles.getVariableDataByName<Vec3d>("Position");
        IndexVector ids;
        for (size_t i = 0; i < bar_particles.TotalRealParticles(); ++i)
        {
            if ((pos[i] - fixed_pos).norm() < 1.3 * dp)
                ids.push_back(i);
        }
        return ids;
    }();
    BodyPartByParticle fixed_part(*bar_object.bar_body);
    fixed_part.body_part_particles_ = fixed_ids;
    SimpleDynamics<slender_structure_dynamics::ConstrainBarBodyRegion> fix_bc(fixed_part);
    simulation.velocity_bc = [&](Real)
    {
        fix_bc.exec();
    };

    // loading condition
    Vec3d elbow_pos = Vec3d::Zero();
    auto tip_id = get_closest_particle(bar_particles, elbow_pos);
    auto *f_prior = bar_particles.getVariableDataByName<Vec3d>("ForcePrior");
    const Real &physical_time = *simulation.system.getSystemVariableDataByName<Real>("PhysicalTime");
    simulation.acceleration_bc = [&](Real)
    {
        f_prior[tip_id] = force(physical_time) * force_direction;
    };

    // Output
    bar_particles.addVariableToWrite<Vec3d>("Velocity");
    bar_particles.addVariableToWrite<Vec3d>("ForcePrior");
    bar_particles.addVariableToWrite<Vec3d>("AngularAcceleration");
    bar_particles.addVariableToWrite<Vec3d>("AngularVelocity");
    bar_particles.addVariableToWrite<Vec3d>("PseudoNormal");
    bar_particles.addVariableToWrite<Vec3d>("PseudoBinormal");
    bar_particles.addVariableToWrite<Mat3d>("NablaDisplacement");
    bar_particles.addVariableToWrite<Vec3d>("DrDs");
    bar_particles.addVariableToWrite<Vec3d>("PositionHourglassForce");
    BodyStatesRecordingToVtp vtp_output(simulation.system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(*bar_object.bar_body);
    vtp_output.writeToFile(0);
    simulation.output_function = [&](size_t ite)
    {
        vtp_output.writeToFile(ite);
    };

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    const Real dt_ref = bar_object.get_time_step_size();
    std::cout << "Reference time step size: " << dt_ref << std::endl;

    simulation.output_number = 100;
    try
    {
        simulation.run_until(end_time, false);
    }
    catch (std::exception &)
    {
        vtp_output.writeToFile();
    }
}