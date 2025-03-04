#include "test_3d_shell_to_solid_coupling.h"

void run_shell_to_solid_coupling(bool two_sided, int res_factor, Real stiffness_ratio, bool use_interpolation, int load_type = 0);

int main(int ac, char *av[])
{
    run_shell_to_solid_coupling(false, 1, 0.1, false, 0);
    // run_solid(1, 0.1, 1);
}

void run_shell_to_solid_coupling(bool two_sided, int res_factor, Real stiffness_ratio, bool use_interpolation, int load_type)
{
    // parameters
    plate_parameters params;
    Real dp = params.height / (10.0 * res_factor);

    // Import input variables
    auto [solid_inputs, shell_inputs] =
        two_sided ? two_sided_problem<SaintVenantKirchhoffSolid>{}(dp, stiffness_ratio) : one_sided_problem<SaintVenantKirchhoffSolid>{}(dp, stiffness_ratio);

    // System
    simulation_system simu_sys(dp);

    // Body
    for (const auto &solid_input : solid_inputs)
        simu_sys.add_solid_object<SaintVenantKirchhoffSolid>(solid_input, use_interpolation);
    for (auto &shell_input : shell_inputs)
        simu_sys.add_shell_object<SaintVenantKirchhoffSolid>(shell_input);

    // Boundary condition
    simu_sys.add_bcs(params.maximum_elongation, params.speed, load_type);

    // Coupling conditions
    simu_sys.add_coupling_algs();
    simu_sys.initialise_coupling_algs(use_interpolation);

    // Initialization
    simu_sys.system_initialize();
    simu_sys.init_config();

    // Output
    simu_sys.add_variable_to_write<int>("IsCoupled");
    simu_sys.add_variable_to_write<Vec3d>("Velocity");
    simu_sys.add_variable_to_write<Vec3d>("Force");
    simu_sys.add_variable_to_write<Vec3d>("ForcePrior");
    BodyStatesRecordingToVtp vtp_output(simu_sys.system);
    for (auto &solid : simu_sys.solid_objects_)
    {
        vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(solid->body_);
        vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(solid->body_);
        vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(solid->body_);
    }
    for (auto &shell : simu_sys.shell_objects_)
    {
        shell->body_.getBaseParticles().addVariableToWrite<Vecd>("CouplingForce");
        vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell->body_);
        vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(shell->body_);
        vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(shell->body_);
    }
    vtp_output.writeToFile(0);

    // check coupling algorithm
    if (!use_interpolation)
    {
        // Check for solid
        auto check_solid = [&](solid_coupling_algs &alg)
        {
            const auto *pos_solid = alg.part.getBaseParticles().getVariableDataByName<Vec3d>("Position");
            for (const auto &index_i : alg.part.body_part_particles_)
            {
                const auto &pos_i = pos_solid[index_i];
                auto [k, j] = alg.vel_bc_nearest_neighbor.get_nearest_id(index_i);
                const auto *pos_shell = alg.contact_relation.contact_particles_[k]->getVariableDataByName<Vec3d>("Position");
                const auto &pos_j = pos_shell[j];
                Vec3d displacement = pos_j - pos_i;
                if (displacement.norm() > Eps)
                {
                    std::cout << alg.part.getSPHBody().getName() << " coupling id is wrong!" << std::endl;
                    exit(0);
                }
            }
        };
        // Check for shells
        auto check_shell_coupling = [&](shell_coupling_algs &alg)
        {
            const auto *pos_shell = alg.part.getBaseParticles().getVariableDataByName<Vec3d>("Position");
            for (const auto &index_i : alg.part.body_part_particles_)
            {
                const auto &pos_i = pos_shell[index_i];
                for (size_t k = 0; k < alg.contact_relation.contact_configuration_.size(); k++)
                {
                    const size_t j = alg.force_bc_nearest_neighbor.get_nearest_id(index_i, k);
                    const auto &pos_j = alg.contact_relation.contact_particles_[k]->getVariableDataByName<Vec3d>("Position")[j];
                    Vec3d displacement = pos_j - pos_i;
                    if (displacement.norm() > Eps)
                    {
                        std::cout << alg.part.getSPHBody().getName() << " coupling id is wrong!" << std::endl;
                        exit(0);
                    }
                }
            }
        };
        for (auto &solid : simu_sys.solid_objects_)
            check_solid(*solid->coupling_algs_);
        for (auto &shell : simu_sys.shell_objects_)
            check_shell_coupling(*shell->coupling_algs_);
    }

    // Simulation
    const Real end_time = 2.0;
    Real &physical_time = *simu_sys.system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = simu_sys.get_time_step_size();
    auto run_simulation = [&]()
    {
        while (physical_time < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                dt = simu_sys.get_time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                // solid 1st half
                for (auto &solid : simu_sys.solid_objects_)
                    solid->algs_->stress_relaxation_first(dt);

                // compute shell coupling force
                for (auto &shell : simu_sys.shell_objects_)
                    shell->coupling_algs_->exec(use_interpolation);

                // update shell
                for (auto &shell : simu_sys.shell_objects_)
                {
                    shell->algs_->stress_relaxation_first(dt);
                    shell->bcs_->exec();
                    shell->algs_->damping_exec(dt);
                    shell->bcs_->exec();
                    shell->algs_->stress_relaxation_second(dt);
                }

                // update solid kinematic constraint and 2nd half
                for (auto &solid : simu_sys.solid_objects_)
                {
                    solid->coupling_algs_->exec(use_interpolation);
                    solid->bcs_->exec();
                    solid->algs_->damping_exec(dt);
                    solid->coupling_algs_->exec(use_interpolation);
                    solid->bcs_->exec();
                    solid->algs_->stress_relaxation_second(dt);
                }

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        vtp_output.writeToFile(ite_output + 1);
        exit(0);
    }
}