#include "test_3d_shell_to_solid_coupling.h"
#include <gtest/gtest.h>

void run_shell_to_solid_coupling(bool two_sided, int res_factor_solid, int res_factor_shell, Real stiffness_ratio, bool use_relaxation, bool use_interpolation, int load_type = 0);
void run_solid(bool two_sided, int res_factor_1, Real stiffness_ratio, int load_type, bool relax = false, bool reload = false);

int main(int ac, char *av[])
{
    for (const auto &two_sided : {false, true})
        for (const auto &load_type : {0, 1})
            for (const auto &res_factor : {1, 2, 3, 4, 5})
                run_shell_to_solid_coupling(two_sided, res_factor, res_factor, 0.1, true, true, load_type);
    // run_mr_solid(false, 1, 1, 0.1, 0);
}

void run_shell_to_solid_coupling(bool two_sided, int res_factor_solid, int res_factor_shell, Real stiffness_ratio, bool use_relaxation, bool use_interpolation, int load_type)
{
    // parameters
    plate_parameters params;
    Real dp_solid = params.height / (10.0 * res_factor_solid);
    Real dp_shell = params.height / (10.0 * res_factor_shell);

    // Import input variables
    auto [solid_inputs, shell_inputs] =
        two_sided ? two_sided_problem<SaintVenantKirchhoffSolid>{}(dp_solid, dp_shell, stiffness_ratio)
                  : one_sided_problem<SaintVenantKirchhoffSolid>{}(dp_solid, dp_shell, stiffness_ratio);

    // System
    simulation_system simu_sys(dp_solid);
    // change output path
    {
        std::string folder = two_sided ? "./output_two_sided" : "./output_one_sided";
        load_type == 0 ? folder += "_extension" : folder += "_shear";
        std::string name = "solid_res_x" + std::to_string(res_factor_solid) + "_shell_res_x" + std::to_string(res_factor_shell);
        std::string path = fs::path(folder) / name;
        fs::remove_all(path);
        fs::create_directories(path);
        simu_sys.io_environment.output_folder_ = path;
    }

    // Body
    for (const auto &solid_input : solid_inputs)
        simu_sys.add_solid_object<SaintVenantKirchhoffSolid>(solid_input, use_relaxation);
    for (auto &shell_input : shell_inputs)
        simu_sys.add_shell_object<SaintVenantKirchhoffSolid>(shell_input);

    // Boundary condition
    simu_sys.add_bcs(params.maximum_elongation, params.speed, load_type);

    // Coupling conditions
    simu_sys.add_coupling_algs(use_interpolation);

    // Initialization
    simu_sys.system_initialize();
    simu_sys.init_config();
    simu_sys.total_weight_initialization();

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
        solid->body_.getBaseParticles().addVariableToWrite<Real>("TotalWeight");
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
    if (!use_interpolation && res_factor_shell == res_factor_solid)
    {
        // Check for solid
        auto check_pos = [&](auto &alg)
        {
            const auto *pos_source = alg.part.getBaseParticles().template getVariableDataByName<Vec3d>("Position");
            for (const auto &index_i : alg.part.body_part_particles_)
            {
                const auto &pos_i = pos_source[index_i];
                size_t neighbor_number = 0;
                for (size_t k = 0; k < alg.contact_relation.contact_particles_.size(); k++)
                {
                    auto neighbor = alg.contact_relation.contact_configuration_[k][index_i];
                    neighbor_number += neighbor.current_size_;
                    for (size_t n = 0; n < neighbor.current_size_; n++)
                    {
                        const auto *pos_target = alg.contact_relation.contact_particles_[k]->template getVariableDataByName<Vec3d>("Position");
                        size_t j = neighbor.j_[n];
                        const auto &pos_j = pos_target[j];
                        Vec3d displacement = pos_j - pos_i;
                        ASSERT_LT(displacement.norm(), Eps);
                    }
                }
                ASSERT_EQ(neighbor_number, 1);
            }
        };
        for (auto &solid : simu_sys.solid_objects_)
            check_pos(*solid->coupling_algs_);
        for (auto &shell : simu_sys.shell_objects_)
            check_pos(*shell->coupling_algs_);
    }

    // check force consistency
    std::vector<std::unique_ptr<ReducedQuantityRecording<InterfaceTotalForce>>> solid_total_force_recorders;
    std::vector<std::unique_ptr<ReducedQuantityRecording<InterfaceTotalForcePrior>>> shell_total_force_recorders;
    for (auto &solid : simu_sys.solid_objects_)
        solid_total_force_recorders.emplace_back(std::make_unique<ReducedQuantityRecording<InterfaceTotalForce>>(solid->body_));
    for (auto &shell : simu_sys.shell_objects_)
        shell_total_force_recorders.emplace_back(std::make_unique<ReducedQuantityRecording<InterfaceTotalForcePrior>>(shell->body_));

    auto record_force = [&](size_t ite)
    {
        for (auto &recorder : solid_total_force_recorders)
            recorder->writeToFile(ite);

        for (auto &recorder : shell_total_force_recorders)
            recorder->writeToFile(ite);
    };

    // Observer
    StdVec<Vecd> observation_locations_x;
    StdVec<Vecd> observation_locations_y;
    {
        // number of particles per line
        size_t number_of_observer_particles = 21;
        // z: position of the lower surface
        const Real z = -0.5 * params.height;
        // add along x direction
        for (size_t i = 0; i < number_of_observer_particles; i++)
        {
            Real x = -0.5 * params.length + i * params.length / (number_of_observer_particles - 1);
            observation_locations_x.emplace_back(x, 0, z);
        }
        // add along y direction
        for (size_t i = 0; i < number_of_observer_particles; i++)
        {
            Real y = -0.5 * params.width + i * params.width / (number_of_observer_particles - 1);
            observation_locations_y.emplace_back(0, y, z);
        }
    }
    observer_object observer_x(simu_sys.get_all_solid_bodies(), observation_locations_x, "ObserverX");
    observer_object observer_y(simu_sys.get_all_solid_bodies(), observation_locations_y, "ObserverY");

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
                    shell->coupling_algs_->exec();

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
                    solid->coupling_algs_->exec();
                    solid->bcs_->exec();
                    solid->algs_->damping_exec(dt);
                    solid->coupling_algs_->exec();
                    solid->bcs_->exec();
                    solid->algs_->stress_relaxation_second(dt);
                }

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);

            record_force(ite_output);
            observer_x.record(ite_output);
            observer_y.record(ite_output);
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

void run_solid(bool two_sided, int res_factor, Real stiffness_ratio, int load_type, bool relax, bool reload)
{
    // parameters
    plate_parameters params;
    Real dp = params.thickness_shell / (4.0 * res_factor);

    // Material properties
    const Real youngs_modulus_solid = params.youngs_modulus_shell * stiffness_ratio;

    // Import model
    auto get_one_sided_mesh = [&]()
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, params.height + 2 * params.thickness_shell);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    };
    auto get_two_sided_mesh = [&]()
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, params.height + params.thickness_shell);
        return makeShared<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize, "solid");
    };
    auto mesh = two_sided ? get_two_sided_mesh() : get_one_sided_mesh();

    // Shell shape
    auto get_shell_shape_one_sides = [&]()
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, params.thickness_shell);
        const Vec3d translation = 0.5 * (params.height + params.thickness_shell) * Vec3d::UnitZ();
        auto shape = makeShared<ComplexShape>("refinement_shape");
        shape->add<TransformShape<GeometricShapeBox>>(Transform(translation), halfsize);
        shape->add<TransformShape<GeometricShapeBox>>(Transform(-translation), halfsize);
        return shape;
    };
    auto get_shell_shape_two_sides = [&]()
    {
        const Vec3d halfsize = 0.5 * Vec3d(params.length, params.width, params.thickness_shell);
        auto shape = makeShared<ComplexShape>("refinement_shape");
        shape->add<TransformShape<GeometricShapeBox>>(Transform(Vec3d::Zero()), halfsize);
        return shape;
    };
    auto shell_mesh = two_sided ? get_shell_shape_two_sides() : get_shell_shape_one_sides();

    // System
    auto bbox = mesh->getBounds();
    SPHSystem system(bbox, dp);
    system.setReloadParticles(reload);
    IOEnvironment io_environment(system);

    // Body
    SolidBody solid_body(system, mesh, "solid");
    solid_body.defineBodyLevelSetShape()->cleanLevelSet(0);
    solid_body.defineMaterial<CompositeSolid>(params.rho);
    system.ReloadParticles() ? solid_body.generateParticles<BaseParticles, Reload>(solid_body.getName()) : solid_body.generateParticles<BaseParticles, Lattice>();
    auto &material = *dynamic_cast<CompositeSolid *>(&solid_body.getBaseMaterial());
    material.add<SaintVenantKirchhoffSolid>(params.rho, youngs_modulus_solid, params.poisson_ratio);
    material.add<SaintVenantKirchhoffSolid>(params.rho, params.youngs_modulus_shell, params.poisson_ratio);

    // Algorithm
    Real length_scale = two_sided ? params.height + params.thickness_shell : params.height + 2 * params.thickness_shell;
    Real eta = get_physical_viscosity_general(params.rho, youngs_modulus_solid, length_scale);
    solid_algs algs_solid(solid_body, eta);

    // Relaxation
    if (relax)
        relax_solid(algs_solid.inner_relation);
    if (!reload)
    {
        ReloadParticleIO write_particle_reload_files(solid_body);
        write_particle_reload_files.writeToFile(0);
    }

    // assign material ids
    SimpleDynamics<SolidMaterialInitialization> material_id_initialization(solid_body, [&](Vec3d &pos) -> bool
                                                                           { return shell_mesh->checkContain(pos); });
    material_id_initialization.exec();

    // Boundary condition
    boundary_condition bc(solid_body, params.maximum_elongation, params.speed, load_type);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    algs_solid.corrected_config();

    // Output
    solid_body.getBaseParticles().addVariableToWrite<int>("MaterialID");
    solid_body.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStress>>(solid_body);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<VonMisesStrain>>(solid_body);
    vtp_output.writeToFile(0);

    // Observer
    StdVec<Vecd> observation_locations_x;
    StdVec<Vecd> observation_locations_y;
    {
        // number of particles per line
        size_t number_of_observer_particles = 21;
        // z: position of the lower surface
        const Real z = two_sided ? -0.5 * (params.height + params.thickness_shell) : -0.5 * params.height;
        // add along x direction
        for (size_t i = 0; i < number_of_observer_particles; i++)
        {
            Real x = -0.5 * params.length + i * params.length / (number_of_observer_particles - 1);
            observation_locations_x.emplace_back(x, 0, z);
        }
        // add along y direction
        for (size_t i = 0; i < number_of_observer_particles; i++)
        {
            Real y = -0.5 * params.width + i * params.width / (number_of_observer_particles - 1);
            observation_locations_y.emplace_back(0, y, z);
        }
    }
    observer_object observer_x({&solid_body}, observation_locations_x, "ObserverX");
    observer_object observer_y({&solid_body}, observation_locations_y, "ObserverY");

    // Simulation
    const Real end_time = 2.0;
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    int ite_output = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    const Real dt_ref = algs_solid.time_step_size();
    auto run_simulation = [&]()
    {
        while (physical_time < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                dt = algs_solid.time_step_size();
                if (dt < dt_ref / 1e2)
                    throw std::runtime_error("time step decreased too much");

                algs_solid.stress_relaxation_first(dt);
                bc.exec();
                algs_solid.damping_exec(dt);
                bc.exec();
                algs_solid.stress_relaxation_second(dt);

                ++ite;
                integral_time += dt;
                physical_time += dt;
            }

            ite_output++;
            vtp_output.writeToFile(ite_output);
            observer_x.record(ite_output);
            observer_y.record(ite_output);
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