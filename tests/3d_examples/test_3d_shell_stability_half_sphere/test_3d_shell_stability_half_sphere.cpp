/**
 * @file 	test_3d_sphere_compression.cpp
 * @brief 	Shell verification  incl. refinement study
 * @details Circular plastic shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH051: Bending of a Circular Plate Using Axisymmetric Elements
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <numeric>

using namespace SPH;

class ShellSphereParticleGenerator : public SurfaceParticleGenerator
{
    const StdVec<Vec3d> &pos_0_;
    const Vec3d center_;
    const Real particle_area_;
    const Real thickness_;

  public:
    explicit ShellSphereParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d> &pos_0, const Vec3d &center, Real particle_area, Real thickness)
        : SurfaceParticleGenerator(sph_body),
          pos_0_(pos_0),
          center_(center),
          particle_area_(particle_area),
          thickness_(thickness){};
    virtual void initializeGeometricVariables() override
    {
        for (const auto &pos : pos_0_)
        {
            Vec3d center_to_pos = pos - center_;
            initializePositionAndVolumetricMeasure(pos, particle_area_);
            initializeSurfaceProperties(center_to_pos.normalized(), thickness_);
        }
    }
};

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

    std::ifstream myfile(file_name, std::ios_base::in);
    if (!myfile.is_open())
        throw std::runtime_error("read_obj_vertices: file doesn't exist");

    StdVec<Vec3d> pos_0;
    Vec3d particle(0);
    unsigned int count = 0;
    Real value = 0;

    while (myfile >> value)
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

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

void sphere_compression(int dp_ratio, Real pressure, Real gravity_z)
{
    // main geometric parameters
    Real scale = 1;
    Real unit_mm = 1e-3;
    Real radius = 50 * scale;   // 50 mm
    Real thickness = 1 * scale; // 1 mm
    Vec3d center(0);
    // resolution
    Real dp = dp_ratio * thickness;
    Real total_area = 0.5 * 4 * Pi * radius * radius;
    std::cout << "total_area: " << total_area << std::endl;
    // material
    Real rho = 1e3 * pow(unit_mm, 3);
    Real E = 5e7 * pow(unit_mm, 2);
    Real mu = 0.3;
    auto material = makeShared<SaintVenantKirchhoffSolid>(rho, E, mu); // NeoHookean always locks one particle with excessive strain
    Real physical_viscosity = 7e3;
    std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
    physical_viscosity = get_physical_viscosity_general(rho, E, thickness);
    std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
    // pressure
    Vec3d gravity = gravity_z * Vec3d(1, 0, 0) / unit_mm;
    // system bounding box
    BoundingBox bb_system;

    // generating particles from predefined positions from obj file
    StdVec<Vec3d> obj_vertices = read_obj_vertices("input/shell_sphere_half_" + std::to_string(dp_ratio) + ".txt");
    std::for_each(obj_vertices.begin(), obj_vertices.end(), [&](Vec3d &vec)
                  { vec *= scale; });
    Real particle_area = total_area / obj_vertices.size();
    // find out BoundingBox
    bb_system = get_particles_bounding_box(obj_vertices);
    std::cout << "bb_system.first_: " << bb_system.first_ << std::endl;
    std::cout << "bb_system.second_: " << bb_system.second_ << std::endl;

    // shell
    auto shell_shape = makeShared<ComplexShape>("shell_shape" + std::to_string(dp_ratio)); // keep all data for parameter study

    // starting the actual simulation
    SPHSystem system(bb_system, dp);
    SolidBody shell_body(system, shell_shape);
    shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
    shell_body.generateParticles<ShellSphereParticleGenerator>(obj_vertices, center, particle_area, thickness);
    auto shell_particles = dynamic_cast<ShellParticles *>(&shell_body.getBaseParticles());
    // output
    IOEnvironment io_env(system, false);
    shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
    BodyStatesRecordingToVtp vtp_output(io_env, {shell_body});
    vtp_output.writeToFile(0);

    // methods
    InnerRelation shell_body_inner(shell_body);
    SimpleDynamics<TimeStepInitialization> initialize_external_force(shell_body, makeShared<Gravity>(gravity));
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(shell_body_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell_body);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_body_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_body_inner);

    // pressure boundary condition
    auto apply_pressure = [&]()
    {
        Real pressure_MPa = pressure * pow(unit_mm, 2);
        for (size_t i = 0; i < shell_particles->acc_prior_.size(); ++i)
        {
            // opposite to normals
            shell_particles->acc_prior_[i] -= pressure_MPa * shell_particles->Vol_[i] / shell_particles->ParticleMass(i) * shell_particles->n_[i];
        }
    };

    BodyPartByParticle constrained_edges(shell_body, "constrained_edges");
    auto constrained_edge_ids = [&]() { // brute force finding the edges
        IndexVector ids;
        for (size_t i = 0; i < shell_body.getBaseParticles().pos_.size(); ++i)
            if (shell_body.getBaseParticles().pos_[i][2] < 0.67 * dp)
                ids.push_back(i);
        return ids;
    }();
    constrained_edges.body_part_particles_ = constrained_edge_ids;

    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(constrained_edges);

    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_velocity_damping(0.2, shell_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_rotation_damping(0.2, shell_body_inner, "AngularVelocity", physical_viscosity);

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();

    {     // tests on initialization
        { // checking particle distances - avoid bugs of reading file
            Real min_rij = Infinity;
            Real max_rij = 0;
            for (size_t i = 0; i < shell_particles->pos0_.size(); ++i)
            {
                Neighborhood &inner_neighborhood = shell_body_inner.inner_configuration_[i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    if (r_ij < min_rij)
                        min_rij = r_ij;
                    if (r_ij > max_rij)
                        max_rij = r_ij;
                }
            }
            std::cout << "min_rij: " << min_rij << std::endl;
            std::cout << "max_rij: " << max_rij << std::endl;
            EXPECT_GT(min_rij, dp / 2);
        }

        // test volume
        Real total_volume = std::accumulate(shell_particles->Vol_.begin(), shell_particles->Vol_.end(), 0.0);
        std::cout << "total_volume: " << total_volume << std::endl;
        Real total_mass = std::accumulate(shell_particles->mass_.begin(), shell_particles->mass_.end(), 0.0);
        std::cout << "total_mass: " << total_mass << std::endl;
        EXPECT_FLOAT_EQ(total_volume, total_area);
        EXPECT_FLOAT_EQ(total_mass, total_area * rho);
    }

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    Real end_time = 0.5; // 1 is better
    Real output_period = end_time / 25.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    /**
     * Main loop
     */
    Real max_dt = 0.0;
    // recording - not pushed to GitHub due to lack of matplotlib there
    StdVec<Real> time, max_displacement, center_deflection;
    try
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                initialize_external_force.exec(dt);
                if (pressure > TinyReal)
                    apply_pressure();

                dt = computing_time_step_size.exec();
                { // checking for excessive time step reduction
                    if (dt > max_dt)
                        max_dt = dt;
                    if (dt < max_dt / 1e3)
                        throw std::runtime_error("time step decreased too much, iteration: " + std::to_string(ite));
                }

                stress_relaxation_first_half.exec(dt);
                constrain_holder.exec();
                shell_velocity_damping.exec(dt);
                shell_rotation_damping.exec(dt);
                constrain_holder.exec();
                stress_relaxation_second_half.exec(dt);

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                // shell_body.updateCellLinkedList();

                { // checking if any position has become nan
                    for (const auto &pos : shell_body.getBaseParticles().pos_)
                        if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                            throw std::runtime_error("position has become nan, iteration: " + std::to_string(ite));
                }
            }
            { // output data
                std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
                vtp_output.writeToFile(ite);
            }
            { // recording - not pushed to GitHub due to lack of matplotlib there
                time.push_back(GlobalStaticVariables::physical_time_);
                max_displacement.push_back(shell_particles->getMaxDisplacement());
            }
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
        vtp_output.writeToFile(ite);
        throw std::runtime_error(e.what());
    }
}

TEST(sphere_compression, half_sphere)
{
    fs::remove_all("output");
    fs::create_directory("output");

    int dp_ratio = 2;         // 4, 2, 1
    Real pressure = 0;        // Pa
    Real gravity_z = -9.8066; // m/s
    EXPECT_NO_THROW(sphere_compression(dp_ratio, pressure, gravity_z));
}

// TEST(sphere_compression, dp_1)
// {
// 	fs::remove_all("output");
// 	fs::create_directory("output");

// 	int dp_ratio = 2;
// 	Real pressure = 1e5; // Pa
// 	Real gravity_z = 0; // m/s
// 	EXPECT_NO_THROW(sphere_compression(dp_ratio, pressure, gravity_z));
// }

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
