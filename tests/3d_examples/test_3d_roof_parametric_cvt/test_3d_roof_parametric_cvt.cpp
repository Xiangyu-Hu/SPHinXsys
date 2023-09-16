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

Real to_rad(Real angle) { return angle * Pi / 180; }

void relax_shell(RealBody &plate_body, Real thickness, Real level_set_refinement_ratio)
{
    // BUG: apparently only works if dp > thickness, otherwise ShellNormalDirectionPrediction::correctNormalDirection() throws error

    InnerRelation imported_model_inner(plate_body);
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(plate_body);
    relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(imported_model_inner, thickness, level_set_refinement_ratio);
    relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(imported_model_inner, thickness);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_step_inner.mid_surface_bounding_.exec();
    plate_body.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
        }
        relaxation_step_inner.exec();
        ite_p += 1;
    }
    shell_normal_prediction.exec();
    std::cout << "The physics relaxation process of imported model finish !" << std::endl;
}

class ShellRoofParticleGenerator : public SurfaceParticleGenerator
{
    const StdVec<Vec3d> &pos_0_;
    const Vec3d center_;
    const Real particle_area_;
    const Real thickness_;

  public:
    explicit ShellRoofParticleGenerator(SPHBody &sph_body, const StdVec<Vec3d> &pos_0, const Vec3d &center, Real particle_area, Real thickness)
        : SurfaceParticleGenerator(sph_body),
          pos_0_(pos_0),
          center_(center),
          particle_area_(particle_area),
          thickness_(thickness){};
    virtual void initializeGeometricVariables() override
    {
        for (const auto &pos : pos_0_)
        {
            // creating the normal direction - z coordinate is always zero
            Vec3d center_to_pos = pos - center_;
            center_to_pos[2] = 0;
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
    variable_sum /= kernel_sum;
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

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

class CylinderParticleGenerator : public SurfaceParticleGenerator
{
    Real particle_number_;

  public:
    explicit CylinderParticleGenerator(SPHBody &sph_body, Real particle_number = 16) : SurfaceParticleGenerator(sph_body),
                                                                                       particle_number_(particle_number){};
    virtual void initializeGeometricVariables() override
    {
        // Real radius = 24.875;								/** Radius of the inner boundary of the cylinder. */
        Real height = 50.0;             /** Height of the cylinder. */
        Real thickness = 0.25;          /** Thickness of the cylinder. */
        Real radius_mid_surface = 25.0; /** Radius of the mid surface. */
        /** Initial reference particle spacing. */
        Real particle_spacing_ref = 2.0 * radius_mid_surface * Pi * 80.0 / 360.0 / (Real)particle_number_;
        std::cout << "particle_spacing_ref: " << particle_spacing_ref << std::endl;
        int BWD = 1;                                /** Width of the boundary layer measured by number of particles. */
        Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
        // the cylinder and boundary
        for (int i = 0; i < particle_number_; i++)
        {
            for (int j = 0; j < (height / particle_spacing_ref + 2 * BWD - 1); j++)
            {
                Real x = radius_mid_surface * cos(50.0 / 180.0 * Pi + (i + 0.5) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_);
                Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
                Real z = radius_mid_surface * sin(50.0 / 180.0 * Pi + (i + 0.5) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_);
                initializePositionAndVolumetricMeasure(Vecd(x, z - radius_mid_surface, y - radius_mid_surface + 1), particle_spacing_ref * particle_spacing_ref);
                Vecd n_0 = Vec3d(x / radius_mid_surface, z / radius_mid_surface, 0.0);
                initializeSurfaceProperties(n_0, thickness);
            }
        }
    }
};

struct return_data
{
    Real displ_y_A;
    Real displ_x_A;
    Real total_area;
    Real total_mass;
    Real dp;

    void write_data_to_txt(const std::string &file_name) const
    {
        std::ofstream myfile;
        myfile.open(file_name);
        myfile << "displ_y_A; displ_x_A; total_area; total_mass; dp\n";
        myfile << displ_y_A << "; " << displ_x_A << "; " << total_area << "; " << total_mass << "; " << dp << "\n";
        myfile.close();
    }
};

return_data roof_under_self_weight(Real dp, bool cvt = true, int particle_number = 16)
{
    // main geometric parameters
    Vec3d tangential_vec(1, 0, 0);
    Vec3d radial_vec(0, 1, 0);
    Vec3d length_vec(0, 0, 1);
    unsigned int tangential_axis = 0;
    unsigned int radial_axis = 1;
    unsigned int length_axis = 2;
    Real radius = 25;
    Real length = 50;
    Real thickness = 0.25;
    Real theta = 40;
    Real theta_radian = to_rad(theta);
    Real arc = radius * theta_radian;
    Vec3d center(0, -radius, 0);
    // observer points A and B
    observer_point_shell point_A;
    observer_point_shell point_B;
    point_A.pos_0 = Vec3d(radius * std::sin(theta_radian), radius * std::cos(theta_radian) - radius, 0);
    point_B.pos_0 = Vec3d::Zero();
    // resolution
    const int dp_cm = dp * 100;
    Real total_area = length * 2 * arc; // accounting for particles being on the edges
    std::cout << "total_area: " << total_area << std::endl;
    // material
    Real rho = 36.7347;
    Real E = 4.32e8;
    Real mu = 0.3;
    auto material = makeShared<LinearElasticSolid>(rho, E, mu);
    Real physical_viscosity = 7e3;
    std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
    // physical_viscosity = 2*get_physical_viscosity_general(rho, E, thickness);
    // std::cout << "physical_viscosity: " << physical_viscosity << std::endl;
    // gravity
    Vec3d gravity = -9.8066 * radial_vec;
    // system bounding box
    BoundingBox bb_system;
    StdVec<Vec3d> obj_vertices;
    Real particle_area = -1; // initialized when CVT-based mesh is used

    if (cvt)
    {
        // generating particles from predefined positions from obj file
        obj_vertices = read_obj_vertices("input/shell_50mm_80d_" + std::to_string(dp_cm) + "cm.txt");
        particle_area = total_area / obj_vertices.size();
        // find out BoundingBox
        bb_system = get_particles_bounding_box(obj_vertices); // store this
    }

    // shell
    auto shell_shape = makeShared<ComplexShape>("shell_shape" + std::to_string(dp_cm)); // keep all data for parameter study

    // starting the actual simulation
    SPHSystem system(bb_system, dp);
    SolidBody shell_body(system, shell_shape);
    shell_body.defineParticlesWithMaterial<ShellParticles>(material.get());
    if (cvt)
        shell_body.generateParticles<ShellRoofParticleGenerator>(obj_vertices, center, particle_area, thickness);
    else
        shell_body.generateParticles<CylinderParticleGenerator>(particle_number);
    auto shell_particles = dynamic_cast<ShellParticles *>(&shell_body.getBaseParticles());
    bb_system = get_particles_bounding_box(shell_particles->pos_);
    system.system_domain_bounds_ = bb_system;
    std::cout << "bb_system.first_: " << bb_system.first_ << std::endl;
    std::cout << "bb_system.second_: " << bb_system.second_ << std::endl;
    { // recalculate the volume/area after knowing the particle positions
      // for (auto& vol: shell_particles->Vol_) vol = total_area / shell_particles->total_real_particles_;
      // for (auto& mass: shell_particles->mass_) mass = total_area*rho / shell_particles->total_real_particles_;
    }
    // output
    IOEnvironment io_env(system, false);
    shell_body.addBodyStateForRecording<Vec3d>("NormalDirection");
    BodyStatesRecordingToVtp vtp_output(io_env, {shell_body});
    vtp_output.writeToFile(0);
    // observer points A & B
    point_A.neighbor_ids = [&]() { // only neighbors on the edges
        IndexVector ids;
        Real smoothing_length = shell_particles->getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
        Real x_min = std::abs(point_A.pos_0[tangential_axis]) - dp / 2;
        for (size_t i = 0; i < shell_particles->pos0_.size(); ++i)
        {
            if ((shell_particles->pos0_[i] - point_A.pos_0).norm() < smoothing_length &&
                std::abs(shell_particles->pos0_[i][tangential_axis]) > x_min)
                ids.push_back(i);
        }
        return ids;
    }();
    point_B.neighbor_ids = [&]() { // full neighborhood
        IndexVector ids;
        Real smoothing_length = shell_particles->getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
        for (size_t i = 0; i < shell_particles->pos0_.size(); ++i)
        {
            if ((shell_particles->pos0_[i] - point_B.pos_0).norm() < smoothing_length)
                ids.push_back(i);
        }
        return ids;
    }();
    point_A.interpolate(*shell_particles);
    point_B.interpolate(*shell_particles);
    // point_A.write_data();
    // point_B.write_data();

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
            if (shell_body.getBaseParticles().pos_[i][length_axis] < bb_system.first_[length_axis] + dp / 2 ||
                shell_body.getBaseParticles().pos_[i][length_axis] > bb_system.second_[length_axis] - dp / 2)
                ids.push_back(i);
        return ids;
    }();
    constrained_edges.body_part_particles_ = constrained_edge_ids;

    SimpleDynamics<solid_dynamics::FixedInAxisDirection> constrain_holder(constrained_edges, length_vec);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_velocity_damping(0.2, shell_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vec3d>>>
        shell_rotation_damping(0.2, shell_body_inner, "AngularVelocity", physical_viscosity);

    /** Apply initial condition. */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    corrected_configuration.exec();

    // TESTS on initialization
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

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    int ite = 0;
    Real end_time = 3.0;
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
            vtp_output.writeToFile(ite);
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
        point_A.interpolate(*shell_particles);
        point_B.interpolate(*shell_particles);
        point_A.write_data();
        point_B.write_data();
    }
    { // testing final values
        Real displ_y_A = -0.3019;
        Real displ_x_A = -0.1593;

        EXPECT_NEAR(point_A.displacement[radial_axis], displ_y_A, std::abs(displ_y_A * 0.06));     // 10% - difficult accuracy due incomplete edge kernels
        EXPECT_NEAR(point_A.displacement[tangential_axis], displ_x_A, std::abs(displ_x_A * 0.11)); // 10% - difficult accuracy due incomplete edge kernels
        std::cout << "point_A y displacement: " << point_A.displacement[radial_axis] << std::endl;
        std::cout << "point_A x displacement: " << point_A.displacement[tangential_axis] << std::endl;
    }
    return_data data;
    data.displ_y_A = point_A.displacement[radial_axis];
    data.displ_x_A = point_A.displacement[tangential_axis];
    data.total_area = total_volume;
    data.total_mass = total_mass;
    data.dp = dp;
    return data;
}

TEST(roof_under_self_weight, dp_4)
{ // for CI
    fs::remove_all("output");
    fs::create_directory("output");

    Real dp = 4;
    auto data = roof_under_self_weight(dp);
    data.write_data_to_txt("roof_under_self_weight" + std::to_string(int(dp * 100)) + "cm.txt");
}

TEST(DISABLED_roof_under_self_weight, parametric_dp)
{ // for proper parameter study
    fs::remove_all("output");
    fs::create_directory("output");

    StdVec<Real> dp_vec = {4, 2, 1}; // ,0.5,0.25
    for (auto dp : dp_vec)
    {
        try
        {
            auto data = roof_under_self_weight(dp);
            data.write_data_to_txt("roof_under_self_weight" + std::to_string(int(dp * 100)) + "cm.txt");
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
    }
}

TEST(DISABLED_roof_under_self_weight, lattice)
{
    fs::remove_all("output");
    fs::create_directory("output");

    Real radius_mid_surface = 25.0;

    StdVec<int> particle_number_vec = {16, 32, 64, 128, 256};
    for (auto particle_number : particle_number_vec)
    {
        try
        {
            Real dp = 2.0 * radius_mid_surface * Pi * 80.0 / 360.0 / (Real)particle_number;
            auto data = roof_under_self_weight(dp, false, particle_number);
            data.write_data_to_txt("roof_under_self_weight_lattice_" + std::to_string(particle_number) + ".txt");
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
