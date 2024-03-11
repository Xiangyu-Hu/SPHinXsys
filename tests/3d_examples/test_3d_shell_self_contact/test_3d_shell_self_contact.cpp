/**
 * @file 	test_3d_shell_self_contact.cpp
 * @brief 	Shell self contact test
 * @author 	Weiyi Kong
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <numeric>

using namespace SPH;

class ShellParticleGenerator : public ParticleGeneratorSurface
{
    const StdVec<Vec3d> &pos_0_;
    const StdVec<Vec3d> &n0_;
    const StdVec<Real> &area_;
    const Real thickness_;
    const bool reset_normal_;

  public:
    explicit ShellParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d> &pos_0, const StdVec<Vec3d> &n0, const StdVec<Real> &area, Real thickness, bool reset_normal)
        : ParticleGeneratorSurface(sph_body),
          pos_0_(pos_0),
          n0_(n0),
          area_(area),
          thickness_(thickness),
          reset_normal_(reset_normal){};
    void initializeGeometricVariables() override
    {
        for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
        {
            auto n = n0_[index_i];
            if (reset_normal_)
                n = n.normalized();
            if (n.norm() < std::numeric_limits<Real>::epsilon())
                std::cout << "Warning: defect particle detected!" << std::endl;
            else
            {
                initializePositionAndVolumetricMeasure(pos_0_[index_i], area_[index_i]);
                initializeSurfaceProperties(n, thickness_);
            };
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

StdVec<Vec3d> read_obj_vectors(const std::string &file_name)
{
    std::cout << "read_obj_vectors" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_obj_vectors: file doesn't exist");

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

    std::cout << "read_obj_vectors finished" << std::endl;
    return pos_0;
}

StdVec<Real> read_obj_scalars(const std::string &file_name)
{
    std::cout << "read_obj_scalars started" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_obj_scalars: file doesn't exist");

    StdVec<Real> area_0;
    Real value = 0;

    while (my_file >> value)
    {
        area_0.push_back(value);
    }

    std::cout << "read_obj_scalars finished" << std::endl;
    return area_0;
}

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

void shell_self_contact(Real shell_resolution, Real shell_thickness, bool reset_damping, bool reset_normal)
{
    // Currently running it in these units
    Real unit_mm = 0.001; // mm, MPa, N, g

    // Global parameters
    Real end_time = 1.0;
    Real scale = 1.0; // in mm currently

    // resolutions
    Real shell_res = shell_resolution * scale;
    Real thickness = shell_thickness * scale;

    // material
    Real rho0_s = 1000.0 * std::pow(unit_mm, 2);
    Real poisson = 0.3;
    Real Youngs_modulus = 1e6 * std::pow(unit_mm, 2);
    Real physical_viscosity = get_physical_viscosity_general(rho0_s, Youngs_modulus, thickness);

    auto material = makeShared<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);

    // system bounding box
    BoundingBox bb_system;

    // generating particles from predefined positions from obj file
    const auto obj_vertices = read_obj_vectors("input/shell_vertices.txt");
    const auto obj_normals = read_obj_vectors("input/shell_normals.txt");
    const auto obj_areas = read_obj_scalars("input/shell_areas.txt");

    // find out BoundingBox
    bb_system = get_particles_bounding_box(obj_vertices); // store this

    // starting the actual simulation
    SPHSystem system(bb_system, shell_res);
    IOEnvironment io_environment(system);

    // solid body
    SolidBody shell_body(system, makeShared<DefaultShape>("ShellBody"));
    shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
    shell_body.generateParticles<ShellParticleGenerator>(obj_vertices, obj_normals, obj_areas, thickness, reset_normal);

    auto shell_particles = dynamic_cast<ShellParticles *>(&shell_body.getBaseParticles());
    bb_system = get_particles_bounding_box(shell_particles->pos_);
    system.system_domain_bounds_ = bb_system;

    // reset to 0
    if (reset_damping)
    {
        for (auto &mat : shell_particles->numerical_damping_scaling_)
            mat(Dimensions - 1, Dimensions - 1) = 0;
    }

    // output
    shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
    BodyStatesRecordingToVtp vtp_output({shell_body});
    vtp_output.writeToFile(0);

    // methods
    InnerRelation shell_body_inner(shell_body);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(shell_body_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell_body);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_body_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_body_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> update_normal(shell_body);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_curvature(shell_body_inner);
    // bc
    auto bc_bot = [&]()
    {
        particle_for(
            par,
            shell_particles->total_real_particles_,
            [&](auto i)
            {
                if (shell_particles->pos0_[i].x() <= bb_system.first_.x() + 0.6 * shell_res)
                    shell_particles->vel_[i] = {2. / 3., 0., 0.};
            });
    };
    auto bc_top = [&]()
    {
        particle_for(
            par,
            shell_particles->total_real_particles_,
            [&](auto i)
            {
                if (shell_particles->pos0_[i].x() >= bb_system.second_.x() - 0.6 * shell_res)
                    shell_particles->vel_[i] = {-2. / 3., 0., 0.};
            });
    };
    // damping
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_velocity_damping(0.2, shell_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_rotation_damping(0.2, shell_body_inner, "AngularVelocity", physical_viscosity);
    // self contact
    ShellSelfContactRelation self_contact_relation(shell_body);
    InteractionDynamics<solid_dynamics::ShellSelfContactDensityUsingDummyParticles> shell_self_contact_density(self_contact_relation);
    InteractionWithUpdate<solid_dynamics::SelfContactForce> shell_self_contact_force(self_contact_relation);

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();
    shell_curvature.compute_initial_curvature();
    self_contact_relation.updateConfiguration();

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    /**
     * Main loop
     */
    Real max_dt = 0.0;
    try
    {
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integral_time = 0.0;
            while (integral_time < output_period)
            {
                if (ite % 10 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                shell_self_contact_density.exec();
                shell_self_contact_force.exec();

                dt = computing_time_step_size.exec();
                { // checking for excessive time step reduction
                    if (dt > max_dt)
                        max_dt = dt;
                    if (dt < max_dt / 1e3)
                        throw std::runtime_error("time step decreased too much");
                }

                stress_relaxation_first_half.exec(dt);
                bc_top();
                bc_bot();
                shell_velocity_damping.exec(dt);
                shell_rotation_damping.exec(dt);
                bc_top();
                bc_bot();
                stress_relaxation_second_half.exec(dt);

                ++ite;
                integral_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                update_normal.exec();
                shell_curvature.exec();
                shell_body.updateCellLinkedList();
                self_contact_relation.updateConfiguration();

                { // checking if any position has become nan
                    for (const auto &pos : shell_body.getBaseParticles().pos_)
                        if (std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]))
                            throw std::runtime_error("position has become nan");
                }
            }

            vtp_output.writeToFile(ite);
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    }
    catch (const std::exception &e)
    {
        shell_body.setNewlyUpdated();
        vtp_output.writeToFile(1e6);
        std::cerr << e.what() << '\n';
    }
}

TEST(shell_self_contact, reset_damping)
{ // for CI
    Real shell_resolution = 0.1;
    Real shell_thickness = 0.1;
    shell_self_contact(shell_resolution, shell_thickness, true, false);
}

TEST(shell_self_contact, reset_normal)
{ // for CI
    Real shell_resolution = 0.1;
    Real shell_thickness = 0.1;
    shell_self_contact(shell_resolution, shell_thickness, false, true);
}

TEST(shell_self_contact, default_case)
{ // for CI
    Real shell_resolution = 0.1;
    Real shell_thickness = 0.1;
    shell_self_contact(shell_resolution, shell_thickness, false, false);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(filter) = "shell_self_contact.default_case";
    return RUN_ALL_TESTS();
}
