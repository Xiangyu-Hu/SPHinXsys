#include "beam_wrapper.h"

void pure_bending_test(int res_factor);
void beam_arc_rotation(int res_factor);

int main(int ac, char *av[])
{
    // beam_arc_rotation(4);
    pure_bending_test(4);
}

void pure_bending_test(int res_factor)
{
    // Global parameters
    const double length = 10.0;
    const double EI = 100.0;
    const double GA = 5000.0;
    const double width = sqrt(6.0 * EI / GA); // radius is determined by the ratio of EI and GA

    const auto g2 = Vec3d::UnitY();
    const auto g3 = Vec3d::UnitZ();

    // resolutions
    const double dp = length / double(10 * res_factor);

    // bar geometry parameters
    BarParameters bar_params;
    bar_params.body_name = "Beam";
    bar_params.hourglass_control = true;
    {
        bar_params.geometry_params.dp_ = dp;
        double x = 0.5 * dp;
        while (x < length)
        {
            Vec3d position = x * Vec3d::UnitX();
            bar_params.geometry_params.pos_.emplace_back(position);
            bar_params.geometry_params.n0_.emplace_back(g3);
            bar_params.geometry_params.b_n0_.emplace_back(g2);
            bar_params.geometry_params.width_.emplace_back(width);
            bar_params.geometry_params.thickness_.emplace_back(width);
            x += dp;
        }
    }

    // material
    const double rho0_s = 1.0;
    const double Youngs_modulus = 7.95e4; // Pa
    const double poisson = 0;

    // Material models
    bar_params.material = std::make_shared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    bar_params.physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, width);
    double c_s = bar_params.material->ReferenceSoundSpeed();

    // load
    const double t_ref = length / (c_s / 20.0);
    const double M = 2.5 * Pi; // Nm

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
    simulation.velocity_bc = [&](Real)
    {
        fix_bc.exec();
    };

    // loading condition
    Vec3d tip_pos = length * Vec3d::UnitX();
    auto tip_id = get_closest_particle(bar_particles, tip_pos);
    auto *m_prior = bar_particles.getVariableDataByName<Vec3d>("ExternalTorquePerUnitLength");
    simulation.acceleration_bc = [&](Real)
    {
        double m_t = M / dp;
        m_prior[tip_id] = m_t * Vec3d::UnitZ();
    };

    // Output
    bar_particles.addVariableToWrite<Vec3d>("Velocity");
    bar_particles.addVariableToWrite<Vec3d>("PriorAngularAcceleration");
    bar_particles.addVariableToWrite<Vec3d>("Force");
    bar_particles.addVariableToWrite<Vec3d>("AngularAcceleration");
    bar_particles.addVariableToWrite<Vec3d>("AngularVelocity");
    bar_particles.addVariableToWrite<Vec3d>("PseudoNormal");
    bar_particles.addVariableToWrite<Vec3d>("PseudoBinormal");
    BodyStatesRecordingToVtp vtp_output(simulation.system);
    vtp_output.writeToFile(0);
    simulation.output_function = [&](size_t ite)
    {
        vtp_output.writeToFile(ite);
    };

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    simulation.output_number = 100;
    double end_time = 50 * t_ref;
    simulation.run_until(end_time, false);
}

void beam_arc_rotation(int res_factor)
{
    // Global parameters
    const double Radius = 100;
    const double arc_angle = Pi / 4.0; // 45 degrees
    const double length = arc_angle * Radius;
    const double width = 1.0;
    const double thickness = 1.0;

    // resolutions
    const double dtheta = arc_angle / (10.0 * res_factor); // angle between particles
    const double dp = length / (10.0 * res_factor);        // distance between particles

    BarParameters bar_params;
    bar_params.body_name = "Beam";
    bar_params.hourglass_control = true;

    // bar geometry parameters
    {
        bar_params.geometry_params.dp_ = dp;
        double theta = 0.5 * dtheta;
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
    const double rho0_s = 1.0;
    const double Youngs_modulus = 1e7;
    const double poisson = 0.0;

    // Material models
    bar_params.material = std::make_shared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    bar_params.physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, thickness);

    // load
    std::vector<Real> force_steps{300, 450, 600};
    size_t loading_steps = force_steps.size();

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
    auto loading_bc = [&](double force)
    {
        force_prior[tip_id] = force * Vec3d::UnitZ();
    };

    // Output
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

    // Threshold for steady state
    const double vel_threshold = 1e-2;
    ReduceDynamics<MaximumSpeed> max_vel(*bar_object.bar_body);

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    const double dt_ref = bar_object.get_time_step_size();
    double dt = 0;

    // Loop over loading steps
    auto loop = [&](double loading_force)
    {
        loading_bc(loading_force);
        dt = bar_object.get_time_step_size();
        if (dt < dt_ref / 1.1)
        {
            throw std::runtime_error("Time step size too small, simulation stopped.");
        }
        bar_object.stress_relaxation_first_half(dt);
        fix_bc.exec();
        bar_object.damping(dt);
        fix_bc.exec();
        bar_object.stress_relaxation_second_half(dt);
    };

    // reference value
    std::vector<Vec3d> pos_ref{
        Vec3d{58.84, 22.33, 40.08},
        Vec3d{52.32, 18.62, 48.39},
        Vec3d{47.23, 15.79, 53.37}};
    auto log_result = [&](size_t n_step)
    {
        const Vec3d pos_ref_n = pos_ref[n_step];
        const Vec3d disp_ref_n = pos_ref_n - tip_pos;
        obs_tip_displacement.writeToFile(n_step + 1);
        const Vec3d disp = obs_tip_pos[0] - tip_pos;
        const double rel_error = (disp - disp_ref_n).norm() / disp_ref_n.norm();

        std::cout << "Reference displacement: " << disp_ref_n.x() << ", " << disp_ref_n.y() << ", " << disp_ref_n.z() << std::endl;
        std::cout << "Computed displacement: " << disp.x() << ", " << disp.y() << ", " << disp.z() << std::endl;
        std::cout << "Relative error in displacement: " << rel_error * 100 << "%" << std::endl;
    };

    auto run_simulation = [&]()
    {
        Real &physical_time = *simulation.system.getSystemVariableDataByName<Real>("PhysicalTime");

        TickCount t1 = TickCount::now();
        TimeInterval interval;

        const size_t min_iterations = 10;  // at least run 10 iterations
        const size_t max_iterations = 1e6; // maximum iterations to avoid infinite loop

        for (size_t n_step = 0; n_step < loading_steps; n_step++)
        {
            bool steady_state = false;
            size_t ite = 0;

            double load_n = force_steps[n_step];

            std::cout << "Loading step: " << n_step + 1 << " starts" << std::endl;
            std::cout << "Applied force: " << load_n << " N" << std::endl;

            while (!steady_state && ite < max_iterations)
            {
                loop(load_n);

                ite++;
                physical_time += dt;
                if (ite > min_iterations)
                    steady_state = max_vel.exec() < vel_threshold;
            }

            TickCount t2 = TickCount::now();

            if (ite >= max_iterations)
            {
                std::cout << "Loading step: " << n_step + 1 << " did not reach steady state after " << ite
                          << " iterations." << std::endl;
                std::cout << "Final maximum velocity: " << max_vel.exec() << " m/s" << std::endl;
                throw std::runtime_error("Loading step did not reach steady state");
            }
            std::cout << "Loading step: " << n_step + 1 << " reached steady state after " << ite << " iterations." << std::endl;

            vtp_output.writeToFile(n_step + 1);
            log_result(n_step);

            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::runtime_error &e)
    {
        vtp_output.writeToFile(loading_steps + 1);
        std::cout << "Simulation failed: " << e.what() << std::endl;
    }
}