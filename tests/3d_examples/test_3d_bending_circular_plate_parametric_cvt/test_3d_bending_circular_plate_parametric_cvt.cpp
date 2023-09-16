/**
 * @file 	test_3d_bending_circular_plate.cpp
 * @brief 	Shell verification  incl. refinement study
 * @details Circular plastic shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH051: Bending of a Circular Plate Using Axisymmetric Elements
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <numeric>

using namespace SPH;

static const Real psi_to_pa = 6894.75729;
static const Real inch_to_m = 0.0254;

class ShellCircleParticleGenerator : public SurfaceParticleGenerator
{
    const StdVec<Vec3d> &pos_0_;
    const Vec3d normal_;
    const Real particle_area_;
    const Real thickness_;

  public:
    explicit ShellCircleParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d> &pos_0, const Vec3d &normal, Real particle_area, Real thickness)
        : SurfaceParticleGenerator(sph_body),
          pos_0_(pos_0),
          normal_(normal),
          particle_area_(particle_area),
          thickness_(thickness){};
    virtual void initializeGeometricVariables() override
    {
        for (const auto &pos : pos_0_)
        {
            initializePositionAndVolumetricMeasure(pos, particle_area_);
            initializeSurfaceProperties(normal_, thickness_);
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
        throw std::runtime_error("read_obj_vertices: file doesn't exist: " + file_name);

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

template <typename VariableType>
VariableType interpolate_observer(
    ShellParticles &particles,
    const IndexVector &neighbor_ids,
    const Vec3d &observer_pos_0,
    std::function<VariableType(size_t)> get_variable_value)
{
    Kernel *kernel_ptr = particles.getSPHBody().sph_adaptation_->getKernel();
    Real smoothing_length = particles.getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
    VariableType variable_sum = VariableType::Zero();
    Real kernel_sum = 0;
    for (auto id : neighbor_ids)
    {
        Real distance = (particles.pos0_[id] - observer_pos_0).norm();
        Real kernel = kernel_ptr->W_3D(distance / smoothing_length);
        kernel_sum += kernel;
        variable_sum += kernel * (get_variable_value(id));
    }
    variable_sum /= (kernel_sum + TinyReal);
    return variable_sum;
}

struct observer_point_shell
{
    Vec3d pos_0;
    IndexVector neighbor_ids;
    Vec3d pos_n;
    Vec3d displacement;
    Vec3d global_shear_stress;
    Mat3d global_stress;
    Mat3d def_gradient;
    Mat3d pk2_stress;
    Mat3d cauchy_stress;

    void interpolate(ShellParticles &particles)
    {
        pos_n = interpolate_observer<Vec3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                            { return (*particles.getVariableByName<Vec3d>("Position"))[id]; });
        displacement = interpolate_observer<Vec3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                                   { return (*particles.getVariableByName<Vec3d>("Displacement"))[id]; });
        global_shear_stress = interpolate_observer<Vec3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                                          { return (*particles.getVariableByName<Vec3d>("GlobalShearStress"))[id]; });
        global_stress = interpolate_observer<Mat3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                                    { return (*particles.getVariableByName<Mat3d>("GlobalStress"))[id]; });
        def_gradient = interpolate_observer<Mat3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                                   { return (*particles.getVariableByName<Mat3d>("DeformationGradient"))[id]; });
        pk2_stress = interpolate_observer<Mat3d>(particles, neighbor_ids, pos_0, [&](size_t id)
                                                 {
			Mat3d F = (*particles.getVariableByName<Mat3d>("DeformationGradient"))[id];
			return particles.elastic_solid_.StressPK2(F, id); });
        cauchy_stress = (1.0 / def_gradient.determinant()) * def_gradient * pk2_stress * def_gradient.transpose();
    }

    void write_data() const
    {
        Vec3d z_dir(0, 0, 1);
        std::cout << std::endl
                  << "===================================================" << std::endl;
        std::cout << "pos_n: " << pos_n << std::endl;
        std::cout << "displacement: " << displacement << std::endl;
        std::cout << "global_shear_stress: " << global_shear_stress << std::endl;
        std::cout << "global_stress: " << global_stress << std::endl;
        std::cout << "pk2_stress: " << pk2_stress << std::endl;
        std::cout << "pk2_z_dir: " << z_dir.dot(pk2_stress * z_dir) << std::endl;
        std::cout << "cauchy_stress: " << cauchy_stress << std::endl;
        std::cout << "cauchy_z_dir: " << z_dir.dot(cauchy_stress * z_dir) << std::endl;
        std::cout << "===================================================" << std::endl
                  << std::endl;
    }
};

struct return_data
{
    Real deflection;
    Real stress_max;

    void write_data_to_txt(const std::string &file_name) const
    {
        std::ofstream myfile;
        myfile.open(file_name);
        myfile << "deflection; stress_max\n";
        myfile << deflection << "; " << stress_max << "\n";
        myfile.close();
    }
};

return_data bending_circular_plate(Real dp_ratio)
{
    // main geometric parameters
    Vec3d sym_vec(0, 0, 1);
    unsigned int sym_axis = 2;
    Real radius = 40 * inch_to_m; // 1.016 [m]
    Real thickness = 1 * inch_to_m;
    // observer point
    observer_point_shell point_center;
    point_center.pos_0 = Vec3d::Zero();
    // resolution
    Real dp = dp_ratio * thickness;
    Real total_area = (radius + dp / 2) * (radius + dp / 2) * Pi; // accounting for particles being on the edges
    std::cout << "total_area: " << total_area << std::endl;
    // material
    Real rho = 1; // unit rho, not specified in test case description
    Real E = 3e7 * psi_to_pa;
    Real mu = 0.3;
    auto material = makeShared<LinearElasticSolid>(rho, E, mu);
    Real physical_viscosity = 7e3; // where is this value coming from?
    // pressure
    Real pressure = 6 * psi_to_pa;
    Vec3d gravity = -pressure / (thickness * rho) * sym_vec; // force/mass simplified by area
    // system bounding box
    BoundingBox bb_system;

    // generating particles from predefined positions from obj file
    StdVec<Vec3d> obj_vertices = read_obj_vertices("input/shell_circle_" + std::to_string(int(dp_ratio * 1e3)) + ".txt");
    { // modifying the vertices such that the most outer particle will be on the edge of the circle due to the CVT algorithm
        Real r_max = 0;
        for (const auto &vertex : obj_vertices)
        {
            Real r_i = vertex.norm();
            if (r_max < r_i)
                r_max = r_i;
        }
        std::cout << "r_max: " << r_max << std::endl;
        for (auto &vertex : obj_vertices)
        {
            vertex *= radius / r_max;
        }
        // also update dp and total area
        dp *= radius / r_max;
        total_area = (radius + dp / 2) * (radius + dp / 2) * Pi; // accounting for particles being on the edges
        std::cout << "total_area new: " << total_area << std::endl;
    }
    Real particle_area = total_area / obj_vertices.size();
    // find out BoundingBox
    bb_system = get_particles_bounding_box(obj_vertices);
    std::cout << "bb_system.first_: " << bb_system.first_ << std::endl;
    std::cout << "bb_system.second_: " << bb_system.second_ << std::endl;

    // shell
    auto shell_shape = makeShared<ComplexShape>("shell_shape" + std::to_string(int(dp_ratio * 1e3))); // keep all data for parameter study

    // starting the actual simulation
    SPHSystem system(bb_system, dp);
    SolidBody shell_body(system, shell_shape);
    shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
    shell_body.generateParticles<ShellCircleParticleGenerator>(obj_vertices, sym_vec, particle_area, thickness);
    auto shell_particles = dynamic_cast<ShellParticles *>(&shell_body.getBaseParticles());
    // output
    IOEnvironment io_env(system, false);
    shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
    BodyStatesRecordingToVtp vtp_output(io_env, {shell_body});
    vtp_output.writeToFile(0);
    // observer point
    point_center.neighbor_ids = [&]() { // full neighborhood
        IndexVector ids;
        Real smoothing_length = shell_particles->getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
        for (size_t i = 0; i < shell_particles->pos0_.size(); ++i)
        {
            if ((shell_particles->pos0_[i] - point_center.pos_0).norm() < 2 * smoothing_length)
                ids.push_back(i);
        }
        return ids;
    }();
    EXPECT_FALSE(point_center.neighbor_ids.empty());
    point_center.interpolate(*shell_particles);

    // methods
    InnerRelation shell_body_inner(shell_body);
    SimpleDynamics<TimeStepInitialization> initialize_external_force(shell_body, makeShared<Gravity>(gravity));
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(shell_body_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell_body);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_body_inner, 3, false);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_body_inner);

    BodyPartByParticle constrained_edges(shell_body, "constrained_edges");
    auto constrained_edge_ids = [&]() { // brute force finding the edges
        IndexVector ids;
        for (size_t i = 0; i < shell_body.getBaseParticles().pos_.size(); ++i)
            if (shell_body.getBaseParticles().pos_[i].norm() > radius - dp / 2)
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

    { // tests on initialization
        // checking particle distances - avoid bugs of reading file
        Real min_rij = Infinity;
        for (size_t index_i = 0; index_i < shell_particles->pos0_.size(); ++index_i)
        {
            Neighborhood &inner_neighborhood = shell_body_inner.inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                if (inner_neighborhood.r_ij_[n] < min_rij)
                    min_rij = inner_neighborhood.r_ij_[n];
        }
        EXPECT_GT(min_rij, dp / 2);

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
    Real end_time = 0.001;
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
                if (ite % 1000 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                initialize_external_force.exec(dt);

                dt = std::min(thickness / dp, Real(0.5)) * computing_time_step_size.exec();
                { // checking for excessive time step reduction
                    if (dt > max_dt)
                        max_dt = dt;
                    if (dt < max_dt / 1e3)
                        throw std::runtime_error("time step decreased too much");
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
                            throw std::runtime_error("position has become nan");
                }
            }
            { // output data
                // std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
                vtp_output.writeToFile(ite);
            }
        }
        TimeInterval tt = TickCount::now() - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        vtp_output.writeToFile(ite);
    }
    { // output data
        std::cout << "max displacement: " << shell_particles->getMaxDisplacement() << std::endl;
        point_center.interpolate(*shell_particles);
        point_center.write_data();
    }
    {                                               // testing final values
        Real deflection_ref = -0.08736 * inch_to_m; // -0.00221894
        std::cout << "deflection_ref: " << deflection_ref << std::endl;

        EXPECT_NEAR(std::abs(point_center.displacement[sym_axis]), std::abs(deflection_ref), std::abs(deflection_ref) * 5e-2); // 5%
                                                                                                                               // EXPECT_NEAR(point_center.stress_max, stress_max_ref, stress_max_ref*1e-2);
    }
    return_data data;
    data.deflection = point_center.displacement[sym_axis];
    // data.stress_max = point_center.stress_max;
    return data;
}

TEST(bending_circular_plate, dp_4)
{ // for CI - lowest resolution
    fs::remove_all("output");
    fs::create_directory("output");

    int dp_ratio = 4;
    auto data = bending_circular_plate(dp_ratio);
    data.write_data_to_txt("bending_circular_plate" + std::to_string(int(dp_ratio * 1e3)) + ".txt");
}

TEST(DISABLED_bending_circular_plate, parametric_dp)
{ // for proper parameter study
    fs::remove_all("output");
    fs::create_directory("output");

    StdVec<Real> dp_vec = {4, 2, 1}; // ,0.5,0.25 --> stop at thickness = dp --> dp_ratio = 1 is the highest resolution
    for (auto dp_ratio : dp_vec)
    {
        auto data = bending_circular_plate(dp_ratio);
        data.write_data_to_txt("bending_circular_plate" + std::to_string(int(dp_ratio * 1e3)) + ".txt");
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
