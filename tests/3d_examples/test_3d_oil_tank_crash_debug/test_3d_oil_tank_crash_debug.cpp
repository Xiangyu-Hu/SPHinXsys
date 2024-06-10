/**
 * @file 	test_3d_roof_analytical.cpp
 * @brief 	Shell verification  incl. refinement study
 * @details Roof shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH069: Barrel Vault Roof Under Self Weight
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <numeric>

using namespace SPH;

void oil_tank_fsi(int res_factor_tank, int res_factor_oil, bool relax = false, bool damping_on = false, bool gravity_on = true, bool fsi_on = true);

TEST(oil_tank_debug, no_relax_no_damping)
{
    oil_tank_fsi(1, 1, false, false);
}

TEST(oil_tank_debug, relax_no_damping)
{
    oil_tank_fsi(1, 1, true, false);
}

TEST(oil_tank_debug, relax_damping)
{
    oil_tank_fsi(1, 1, true, true);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(filter) = "oil_tank_debug.relax_damping";
    return RUN_ALL_TESTS();
}

class ShellTankParticleGenerator : public ParticleGenerator<Surface>
{
    const StdVec<Vec3d> &pos_0_;
    const StdVec<Vec3d> &normal_;
    const StdVec<Real> &particle_area_;
    const Real shell_thickness_;

  public:
    explicit ShellTankParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d> &pos_0, const StdVec<Vec3d> &normal, const StdVec<Real> &particle_area, Real thickness)
        : ParticleGenerator<Surface>(sph_body),
          pos_0_(pos_0),
          normal_(normal),
          particle_area_(particle_area),
          shell_thickness_(thickness){};
    void initializeGeometricVariables() override
    {
        ASSERT_EQ(pos_0_.size(), normal_.size());
        ASSERT_EQ(pos_0_.size(), particle_area_.size());
        for (size_t index_i = 0; index_i < pos_0_.size(); ++index_i)
        {
            initializePositionAndVolumetricMeasure(pos_0_[index_i], particle_area_[index_i]);
            initializeSurfaceProperties(normal_[index_i].normalized(), shell_thickness_);
        }
    }
};

void relax_body(RealBody &body, InnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 200 == 0)
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
}

template <typename VectorType>
BoundingBox get_particles_bounding_box(const VectorType &pos_0)
{
    Vec3d lower(pos_0[0]);
    Vec3d upper(pos_0[0]);
    for (const auto &pos : pos_0)
    {
        for (int i = 0; i < 3; i++)
        {
            if (lower[i] > pos[i])
                lower[i] = pos[i];
            if (upper[i] < pos[i])
                upper[i] = pos[i];
        }
    }
    return BoundingBox(lower, upper);
}

StdVec<Vec3d> read_obj_vertices(const std::string &file_name)
{
    std::cout << "read_obj_vertices started" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_obj_vertices: file doesn't exist");

    StdVec<Vec3d> pos_0;
    Vec3d particle(0);
    unsigned int count = 0;
    Real value = 0;

    while (my_file >> value)
    {
        particle[count] = value;
        ++count;
        if (count % 3 == 0)
        {
            count = 0;
            pos_0.push_back(particle);
        }
    }

    std::cout << "read_obj_vertices finished" << std::endl;
    return pos_0;
}

StdVec<Real> read_obj_scalars(const std::string &file_name)
{
    std::cout << "read_obj_scalars started" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_obj_scalars: file doesn't exist");

    StdVec<Real> area;
    Real value = 0;

    while (my_file >> value)
    {
        area.push_back(value);
    }

    std::cout << "read_obj_scalars finished" << std::endl;
    return area;
}

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

void oil_tank_fsi(int res_factor_tank, int res_factor_oil, bool relax, bool damping_on, bool gravity_on, bool fsi_on)
{
    // main geometric parameters
    const Real tank_diameter = 2.5;
    const Real tank_length = 2 * tank_diameter;
    const Real tank_thickness = 11e-3; // 11mm
    const Real oil_height = 0.5 * tank_diameter;

    // resolutions
    const Real res_ref = tank_diameter / 10.0;
    const Real res_tank = res_ref / res_factor_tank;
    const Real res_oil = res_ref / res_factor_oil;

    const Real oil_diameter = tank_diameter - res_tank;
    const Real oil_length = tank_length - 2 * res_tank;

    // end time
    const Real end_time = 0.05;

    // material
    const Real rho0_s = 7850.0;
    const Real poisson = 0.28;
    const Real Youngs_modulus = 206e9;
    const Real physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, tank_thickness);
    SaintVenantKirchhoffSolid material_tank(rho0_s, Youngs_modulus, poisson);

    const Real rho0_f = 700.0;
    const Real mu_f = 0.1;
    const Real gravity_g = 9.8;
    const Vec3d velocity_truck = 20 * Vec3d::UnitZ(); // 20 m/s
    const Real U_f = std::max(Real(2.0 * sqrt(gravity_g * oil_height)), velocity_truck.norm());
    const Real c_f = 10.0 * U_f;

    // Import meshes
    const auto tank_translation = Vec3d::Zero();
    auto mesh_tank = makeShared<TriangleMeshShapeCylinder>(
        SimTK::UnitVec3(xAxis), 0.5 * tank_diameter, 0.5 * tank_length, 10, tank_translation, "Tank");
    const auto oil_cylinder_translation = Vec3d::Zero();
    const Vec3d oil_block_halfsize = 0.5 * Vec3d(oil_length, oil_diameter - oil_height, oil_diameter) + res_oil * Vec3d::Ones();
    const auto oil_block_translation = Vec3d(0, -0.5 * oil_diameter + oil_height + oil_block_halfsize.y(), 0);
    auto mesh_oil = makeShared<ComplexShape>("oil");
    mesh_oil->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(xAxis), 0.5 * oil_diameter, 0.5 * oil_length, 10, oil_cylinder_translation);
    mesh_oil->subtract<TriangleMeshShapeBrick>(oil_block_halfsize, 10, oil_block_translation);

    // system bounding box
    const BoundingBox bb_system = mesh_tank->getBounds();

    // System
    SPHSystem system(bb_system, res_tank);
    IOEnvironment io_environment(system);

    // Create objects
    // tank
    SolidBody tank_body(system, mesh_tank);
    tank_body.defineParticlesWithMaterial<ShellParticles>(&material_tank);
    // generating particles from predefined positions from obj file
    StdVec<Vec3d> obj_vertices = read_obj_vertices("input/oil_tank_x" + std::to_string(res_factor_tank) + "_pos.txt");
    StdVec<Vec3d> obj_normals = read_obj_vertices("input/oil_tank_x" + std::to_string(res_factor_tank) + "_normal.txt");
    StdVec<Real> obj_areas = read_obj_scalars("input/oil_tank_x" + std::to_string(res_factor_tank) + "_area.txt");
    ShellTankParticleGenerator tank_particle_generator(tank_body, obj_vertices, obj_normals, obj_areas, tank_thickness);
    tank_body.generateParticles(tank_particle_generator);
    auto particles_tank = dynamic_cast<ShellParticles *>(&tank_body.getBaseParticles());

    // oil
    FluidBody oil_block(system, mesh_oil);
    oil_block.defineBodyLevelSetShape();
    oil_block.defineAdaptationRatios(1.15, res_tank / res_oil);
    oil_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    oil_block.generateParticles<Lattice>();

    // Inner relations
    InnerRelation tank_inner(tank_body);
    InnerRelation oil_block_inner(oil_block);

    // relax
    if (relax)
        relax_body(oil_block, oil_block_inner);

    // fluid-shell contact
    ContactRelationToShell oil_tank_contact(oil_block, {&tank_body}, {false});
    ContactRelationFromShell tank_oil_contact(tank_body, {&oil_block}, {false});
    ShellInnerRelationWithContactKernel tank_curvature_inner(tank_body, oil_block);
    ComplexRelation oil_block_complex(oil_block_inner, {&oil_tank_contact});

    // Methods
    // tank
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration_tank(tank_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size_tank(tank_body);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half_tank(tank_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half_tank(tank_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal_tank(tank_body);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> tank_average_curvature(tank_curvature_inner);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        velocity_damping(0.2, tank_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        rotation_damping(0.2, tank_inner, "AngularVelocity", physical_viscosity);

    // fluid
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity(oil_block, gravity);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(oil_block_inner, oil_tank_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(oil_block_inner, oil_tank_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(oil_block_inner, oil_tank_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(oil_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(oil_block);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(oil_block_inner, oil_tank_contact);

    // FSI
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_tank(tank_oil_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_on_tank(tank_oil_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(tank_body);

    // Boundary condition
    StdLargeVec<int> constraint_indicator;
    particles_tank->registerVariable(constraint_indicator, "ConstraintIndicator");
    auto fix_id = [&]()
    {
        const Real dx = (tank_length - 4 * res_tank) / 3.0;
        const Real y_min = -0.25 * tank_diameter;
        std::vector<Real> pos_fix_x = {-1.5 * dx, -0.5 * dx, 0.5 * dx, 1.5 * dx};
        IndexVector ids;
        for (size_t index_i = 0; index_i < particles_tank->total_real_particles_; index_i++)
        {
            for (const Real xi : pos_fix_x)
            {
                if (std::abs(particles_tank->pos_[index_i].x() - xi) < res_tank &&
                    particles_tank->pos_[index_i].y() < y_min)
                {
                    ids.push_back(index_i);
                    constraint_indicator[index_i] = 1;
                }
            }
        }
        return ids;
    }();
    auto fix_bc = [&]()
    {
        particle_for(
            execution::ParallelPolicy(),
            fix_id,
            [&](size_t index_i)
            {
                particles_tank->vel_[index_i] = Vec3d::Zero();
                particles_tank->angular_vel_[index_i] = Vec3d::Zero();
            });
    };

    // output
    oil_block.addBodyStateForRecording<Real>("Pressure");
    particles_tank->addVariableToWrite<Real>("Average1stPrincipleCurvature");
    particles_tank->addVariableToWrite<Real>("Average2ndPrincipleCurvature");
    particles_tank->addVariableToWrite<int>("ConstraintIndicator");
    BodyStatesRecordingToVtp vtp_output(system.real_bodies_);

    // Initialization
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration_tank.exec();
    tank_average_curvature.exec();
    oil_block_complex.updateConfiguration();
    tank_oil_contact.updateConfiguration();
    if (gravity_on)
        constant_gravity.exec();

    vtp_output.writeToFile(0);

    // Simulation
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    int ite_output = 0;
    int screen_output_interval = 1;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    Real dt_s = 0.0;
    TickCount t1 = TickCount::now();
    const Real Dt_ref = get_fluid_advection_time_step_size.exec();
    const Real dt_ref = get_fluid_time_step_size.exec();
    const Real dt_s_ref = computing_time_step_size_tank.exec();
    auto run_simulation = [&]()
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                Real Dt = get_fluid_advection_time_step_size.exec();
                if (Dt < Dt_ref / 1e2)
                    throw std::runtime_error("Advective time step decreased too much");
                update_density_by_summation.exec();
                viscous_acceleration.exec();
                if (fsi_on)
                    viscous_force_on_tank.exec();

                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    pressure_relaxation.exec(dt);
                    if (fsi_on)
                        pressure_force_on_tank.exec();
                    density_relaxation.exec(dt);

                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        dt_s = computing_time_step_size_tank.exec();
                        if (dt_s < dt_s_ref / 1e2)
                            throw std::runtime_error("time step decreased too much");

                        stress_relaxation_first_half_tank.exec(dt_s);
                        fix_bc();
                        if (damping_on)
                        {
                            velocity_damping.exec(dt_s);
                            rotation_damping.exec(dt_s);
                            fix_bc();
                        }
                        stress_relaxation_second_half_tank.exec(dt_s);

                        dt_s_sum += dt_s;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);

                    dt = get_fluid_time_step_size.exec();
                    if (dt < dt_ref / 1e2)
                        throw std::runtime_error("Acoustic time step decreased too much");

                    relaxation_time += dt;
                    integral_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                }

                if (ite % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
                }
                ++ite;

                tank_body.updateCellLinkedList();
                oil_block.updateCellLinkedListWithParticleSort(100);

                update_normal_tank.exec();
                tank_curvature_inner.updateConfiguration();
                tank_average_curvature.exec();
                tank_oil_contact.updateConfiguration();
                oil_block_complex.updateConfiguration();
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
        tank_body.setNewlyUpdated();
        vtp_output.writeToFile(ite_output + 1);
        exit(0);
    }
}
