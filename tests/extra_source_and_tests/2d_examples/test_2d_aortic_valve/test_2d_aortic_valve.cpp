/*
 * @file modified_T_shaped_pipe.cpp
 * @brief This is the benchmark test of multi -inlet and multi - outlet.
 * @details We consider a flow with one inlet and two outlets in a T - shaped pipe in 2D.
 * @author Xiangyu Hu,Shuoguo Zhang
 */
#include "test_2d_aortic_valve.h" // case file to setup the test case
#include "bidirectional_buffer.h"
#include "density_correction.h"
#include "density_correction.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include <unsupported/Eigen/Splines>

void run_fsi2();

int main(int ac, char *av[])
{
    run_fsi2();
}

namespace SPH
{
class Shell;
template <>
class ParticleGenerator<SurfaceParticles, Shell> : public ParticleGenerator<SurfaceParticles>
{
    const StdVec<Vec2d> &positions_;
    const StdVec<Vec2d> &normals_;
    Real dp_;
    Real thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               const StdVec<Vec2d> &positions,
                               const StdVec<Vec2d> &normals,
                               Real dp,
                               Real thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          positions_(positions), normals_(normals), dp_(dp), thickness_(thickness)
    {
        if (positions_.size() != normals_.size())
        {
            std::cout << "Error: In ParticleGenerator<Shell>, positions size is not equal to normals size!" << std::endl;
            exit(1);
        }
    };
    void prepareGeometricData() override
    {
        const auto particle_number = positions_.size();
        // generate particles for the elastic gate
        for (size_t i = 0; i < particle_number; i++)
        {
            addPositionAndVolumetricMeasure(positions_[i], dp_);
            addSurfaceProperties(normals_[i], thickness_);
        }
    }
};
} // namespace SPH

void run_relaxation(InnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random(inner.getSPHBody());
    RelaxationStepLevelSetCorrectionInner relaxation_step_complex(inner);
    ReloadParticleIO write_particle_reload_files(inner.getSPHBody());
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random.exec(0.25);
    relaxation_step_complex.SurfaceBounding().exec();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_complex.exec();
        ite_p += 1;
        if (ite_p % 200 == 0)
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
    }
    std::cout << "The physics relaxation process finish !" << std::endl;
    /** Output results. */
    write_particle_reload_files.writeToFile(0);
}

class ShellFluidMixtureMass : public LocalDynamics
{
  private:
    Real rho_f0_; // assume the density of fluid is constant for now
    Real dp_;     // initial particle spacing
    Real *thickness_;
    Real *mass_;
    Real *Vol_;

  public:
    ShellFluidMixtureMass(SPHBody &shell_body, Real rho_f0)
        : LocalDynamics(shell_body),
          rho_f0_(rho_f0),
          dp_(shell_body.getSPHAdaptation().ReferenceSpacing()),
          thickness_(shell_body.getBaseParticles().getVariableDataByName<Real>("Thickness")),
          mass_(shell_body.getBaseParticles().getVariableDataByName<Real>("Mass")),
          Vol_(shell_body.getBaseParticles().getVariableDataByName<Real>("VolumetricMeasure"))
    {
    }

    void update(size_t index_i, Real)
    {
        Real dp_m_t = dp_ - thickness_[index_i];
        if (dp_m_t < 0)
            throw std::runtime_error("Error: In ShellFluidMixtureMass, dp - thickness < 0!");
        Real V_f = Vol_[index_i] * dp_m_t; // fluid volume
        Real mass_f = V_f * rho_f0_;
        mass_[index_i] += mass_f;
    }
};

struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    explicit LeftInflowPressure(BoundaryConditionType &) {}

    Real operator()(Real p, Real)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    explicit RightInflowPressure(BoundaryConditionType &) {}

    Real operator()(Real p, Real)
    {
        return p;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
// Colume 1: time, column 2: flow rate
std::vector<std::vector<Real>> read_csv_data(const std::string &file_name)
{
    std::cout << "read_csv_data started" << std::endl;

    std::ifstream my_file(file_name, std::ios_base::in);
    if (!my_file.is_open())
        throw std::runtime_error("read_csv_data: file doesn't exist");

    std::vector<std::vector<Real>> data;

    std::string line;
    while (std::getline(my_file, line))
    {
        std::vector<Real> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ','))
            row.emplace_back(std::stod(cell));

        data.push_back(row);
    }
    my_file.close();
    // convert to column-wise data
    std::vector<std::vector<Real>> column_data;
    if (!data.empty())
    {
        size_t num_columns = data[0].size();
        column_data.resize(num_columns);
        for (const auto &row : data)
        {
            for (size_t i = 0; i < num_columns; ++i)
            {
                column_data[i].push_back(row[i]);
            }
        }
    }
    std::cout << "read_csv_data finished" << std::endl;
    return column_data;
}

struct InflowVelocity
{
    using Spline1d = Eigen::Spline<double, 1>;
    Spline1d Q_splines_{};

    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &boundary_condition)
    {
        auto get_flowrate_splines = []()
        {
            // Read data from file
            std::string filename = "input/flowrate.csv";
            auto csv_data = read_csv_data(filename);
            auto &time_data = csv_data[0];
            auto &flow_rate_data = csv_data[1];

            // Convert to Eigen
            Eigen::VectorXd X = Eigen::Map<const Eigen::VectorXd>(time_data.data(), time_data.size());
            Eigen::VectorXd flow_rate_vec = Eigen::Map<const Eigen::VectorXd>(flow_rate_data.data(), flow_rate_data.size());

            // Normalize the time from 0-1 sec (could also do this in the csv file)
            double x_min = X.minCoeff();
            double x_max = X.maxCoeff();
            Eigen::VectorXd t = (X.array() - x_min) / (x_max - x_min);

            // Fit cubic spline
            Spline1d spline = Eigen::SplineFitting<Spline1d>::Interpolate(
                flow_rate_vec.transpose(), // column vector
                3,                         // cubic
                t);

            return spline;
        };
        Q_splines_ = get_flowrate_splines();
    }

    Vecd operator()(Vecd &pos, Vecd &, Real time)
    {
        if (time < time_flow_init)
            return Vec2d::Zero();
        Real t = fmod(time - time_flow_init, time_cycle) / time_cycle;
        Real Q = Q_splines_(t)(0) * L_to_m3 / min_to_s;
        Real U_ave = Q / (width * height);
        // parabolic velocity profile
        Real y = pos.y();
        Real u = 1.5 * U_ave * (1.0 - 4.0 * y * y / (height * height));
        return u * Vec2d::UnitX();
    }
};

struct OutflowVelocity
{
    template <class BoundaryConditionType>
    explicit OutflowVelocity(BoundaryConditionType &boundary_condition)
    {
    }

    Vecd operator()(Vecd &, Vecd &vel, Real)
    {
        Vecd target_velocity = vel;
        target_velocity.y() = 0.0;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
void run_fsi2()
{
    std::cout << "dp_shell/thickness = " << dp_shell / shell_thickness << std::endl;
    std::cout << "U_max = " << U_max << std::endl;
    const Real Re = rho0_f * U_f * height / mu_f;
    std::cout << "The Reynolds number is " << Re << std::endl;

    // create shape
    // fluid shape
    std::cout << "Creating fluid shape ... " << std::endl;
    MultiPolygon fluid_shape_poly;
    Vec2d fluid_halfsize = 0.5 * Vec2d(length + 2 * buffer_length, height);
    Vec2d translation = 0.5 * Vec2d(length, height);
    fluid_shape_poly.addBox(Transform(translation), fluid_halfsize, GeometricOps::add);
    fluid_shape_poly.addCircle(circle_center, radius, 10, GeometricOps::add);
    MultiPolygonShape fluid_shape(fluid_shape_poly, "Fluid");
    std::cout << "Creating fluid shape - done. " << std::endl;

    // wall shape
    std::cout << "Creating wall shape ... " << std::endl;
    MultiPolygon wall_shape_poly;
    Vec2d wall_outer_halfsize = fluid_halfsize + Vec2d(dp_fluid, wall_thickness);
    Vec2d wall_inner_halfsize = fluid_halfsize + dp_fluid * Vec2d::UnitX();
    Real wall_radius = radius + wall_thickness;
    wall_shape_poly.addBox(Transform(translation), wall_outer_halfsize, GeometricOps::add);
    wall_shape_poly.addCircle(circle_center, wall_radius, 20, GeometricOps::add);
    wall_shape_poly.addBox(Transform(translation), wall_inner_halfsize, GeometricOps::sub);
    wall_shape_poly.addCircle(circle_center, radius, 20, GeometricOps::sub);
    MultiPolygonShape wall_shape(wall_shape_poly, "Wall");
    std::cout << "Creating wall shape - done. " << std::endl;

    // shell
    StdVec<Vec2d> positions_shell;
    {
        Real x0 = circle_center.x() - radius;
        Real y0 = circle_center.y();
        Real ksi = 0;
        while (ksi < shell_length)
        {
            Real x = x0 + ksi * std::cos(angle);
            Real y = y0 - ksi * std::sin(angle);
            positions_shell.emplace_back(x, y);
            ksi += dp_shell;
        }
    }
    StdVec<Vec2d> normals_shell(positions_shell.size(), Vec2d(std::sin(angle), std::cos(angle)));
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    const auto bbox = wall_shape.getBounds();
    SPHSystem sph_system(bbox, dp_fluid);
    sph_system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);      // Tag for computation with save particles distribution
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody blood_body(sph_system, fluid_shape);
    blood_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    blood_body.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    blood_body.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, dp_fluid / dp_solid);
    wall_boundary.defineBodyLevelSetShape().writeLevelSet();
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.Name())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody shell_body(sph_system, makeShared<DefaultShape>("shell"));
    shell_body.defineAdaptation<SPHAdaptation>(1.15, dp_fluid / dp_shell);
    shell_body.defineMatterMaterial<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
    shell_body.generateParticles<SurfaceParticles, Shell>(positions_shell, normals_shell, dp_shell, shell_thickness);

    // reset mass
    SimpleDynamics<ShellFluidMixtureMass> reset_shell_mass(shell_body, rho0_f);
    reset_shell_mass.exec();

    // run relaxation
    InnerRelation wall_boundary_inner(wall_boundary);
    if (sph_system.RunParticleRelaxation())
        run_relaxation(wall_boundary_inner);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation blood_inner(blood_body);
    InnerRelation shell_inner(shell_body);
    // Upper and lower walls are slip walls
    ContactRelation blood_wall_contact(blood_body, {&wall_boundary});
    ContactRelationFSI2 blood_shell_contact(blood_body, {&shell_body});
    ContactRelationSFI2 shell_blood_contact(shell_body, {&blood_body});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation blood_body_complex(blood_inner, {&blood_wall_contact, &blood_shell_contact});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first_half(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second_half(shell_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_computing_time_step_size(shell_body);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell_body);
    const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * shell_thickness * shell_thickness;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);

    auto clamp_id = [&]()
    {
        IndexVector ids;
        const auto *pos = shell_body.getBaseParticles().getVariableDataByName<Vec2d>("Position");
        auto fix_point = Vec2d(circle_center.x() - radius, circle_center.y());
        for (size_t i = 0; i < shell_body.getBaseParticles().TotalRealParticles(); ++i)
        {
            if ((pos[i] - fix_point).norm() < dp_shell)
                ids.push_back(i);
        }
        return ids;
    }();
    BodyPartByParticle shell_base(shell_body);
    shell_base.body_part_particles_ = clamp_id;
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_shell_base(shell_base);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<ComplexInteraction<NablaWV<Inner<>, Contact<>, Contact<>>>> kernel_summation(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> boundary_indicator(blood_inner, blood_wall_contact, blood_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(blood_inner, blood_wall_contact, blood_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SPHAdaptation, NoLimiter>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_correction(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity, NoKernelCorrection>> viscous_force(blood_inner, blood_wall_contact, blood_shell_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(blood_body, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(blood_body);

    //-----------Inflow and periodic condition --------//
    Vec2d bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_length, 1.2 * height);
    Vec2d left_bidirectional_translation = 0.5 * Vec2d(-buffer_length, height);
    Vec2d right_bidirectional_translation(length + 0.5 * buffer_length, 0.5 * height);
    OrientedBoxByCell left_emitter(blood_body, OrientedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    OrientedBoxByCell right_emitter(blood_body, OrientedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationPressureComplex<Inner<>, Contact<>, Contact<>>> update_fluid_density(blood_inner, blood_wall_contact, blood_shell_contact);

    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<OutflowVelocity>> outflow_velocity_condition(right_emitter);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_body);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(shell_blood_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(shell_blood_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(blood_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vec2d>(blood_body, "Velocity");
    write_real_body_states.addToWrite<Real>(blood_body, "Density");
    write_real_body_states.addToWrite<Real>(blood_body, "Pressure");
    write_real_body_states.addToWrite<int>(blood_body, "Indicator");
    write_real_body_states.addToWrite<int>(blood_body, "BufferIndicator");
    write_real_body_states.addToWrite<Vec2d>(wall_boundary, "NormalDirection");
    write_real_body_states.addToWrite<Vec2d>(shell_body, "NormalDirection");
    write_real_body_states.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    // buffer and surface
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    left_bidirection_buffer.deletion.exec();
    right_bidirection_buffer.deletion.exec();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    shell_corrected_configuration.exec();
    //----------------------------------------------------------------------
    // initial relaxation of fluid body
    //----------------------------------------------------------------------
    auto run_fluid_relaxation = [&]()
    {
        size_t relaxation_fluid_itr = 0;
        std::cout << "Fluid relaxation starting..." << std::endl;
        while (relaxation_fluid_itr < 100)
        {
            boundary_indicator.exec();

            transport_correction.exec();
            relaxation_fluid_itr++;

            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            blood_body.updateCellLinkedList();
            blood_body_complex.updateConfiguration();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        std::cout << "Fluid relaxation finished !" << std::endl;
    };
    run_fluid_relaxation();
    write_real_body_states.writeToFile();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = time_flow_init + time_cycle * num_cycles;
    Real output_interval = time_cycle / 100.0; // output 50 frames per cycle
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    auto run_fsi = [&]()
    {
        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                Real Dt = get_fluid_advection_time_step_size.exec();
                update_fluid_density.exec();
                viscous_force.exec();
                transport_correction.exec();

                if (physical_time > time_flow_init)
                {
                    /** FSI for viscous force. */
                    viscous_force_from_fluid.exec();
                    /** Update normal direction on elastic body.*/
                    shell_update_normal.exec();
                }

                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                    /** Fluid pressure relaxation */
                    pressure_relaxation.exec(dt);

                    // Boundary conditions
                    kernel_summation.exec();
                    inflow_velocity_condition.exec();
                    outflow_velocity_condition.exec();

                    if (physical_time > time_flow_init)
                    {
                        /** FSI for pressure force. */
                        pressure_force_from_fluid.exec();
                    }

                    /** Fluid density relaxation */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    if (physical_time > time_flow_init)
                    {
                        inner_ite_dt_s = 0;
                        Real dt_s_sum = 0.0;
                        average_velocity_and_acceleration.initialize_displacement_.exec();
                        while (dt_s_sum < dt)
                        {
                            Real dt_s = SMIN(shell_computing_time_step_size.exec(), dt - dt_s_sum);
                            shell_stress_relaxation_first_half.exec(dt_s);
                            constraint_shell_base.exec();
                            shell_position_damping.exec(dt_s);
                            shell_rotation_damping.exec(dt_s);
                            constraint_shell_base.exec();
                            shell_stress_relaxation_second_half.exec(dt_s);
                            dt_s_sum += dt_s;
                            inner_ite_dt_s++;
                        }
                        average_velocity_and_acceleration.update_averages_.exec(dt);
                    }

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                    inner_ite_dt++;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time
                              << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt
                              << "	dt / dt_s = " << inner_ite_dt_s
                              << "\n";
                }
                number_of_iterations++;

                /** Water block configuration and periodic condition. */
                // first do injection for all buffers
                left_bidirection_buffer.injection.exec();
                right_bidirection_buffer.injection.exec();
                // then do deletion for all buffers
                left_bidirection_buffer.deletion.exec();
                right_bidirection_buffer.deletion.exec();

                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                blood_body.updateCellLinkedList();
                if (physical_time > time_flow_init)
                    shell_body.updateCellLinkedList();
                blood_body_complex.updateConfiguration();
                shell_blood_contact.updateConfiguration();

                boundary_indicator.exec();
                left_bidirection_buffer.tag_buffer_particles.exec();
                right_bidirection_buffer.tag_buffer_particles.exec();
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            shell_body.updateCellLinkedList();
            write_real_body_states.writeToFile();
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
        run_fsi();
    }
    catch (const std::exception &e)
    {
        std::cerr << "An error occurred during the simulation: " << e.what() << std::endl;
        write_real_body_states.writeToFile();
    }
}
