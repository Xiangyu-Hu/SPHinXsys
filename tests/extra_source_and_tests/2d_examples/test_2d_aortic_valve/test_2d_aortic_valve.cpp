/*
 * @file modified_T_shaped_pipe.cpp
 * @brief This is the benchmark test of multi -inlet and multi - outlet.
 * @details We consider a flow with one inlet and two outlets in a T - shaped pipe in 2D.
 * @author Xiangyu Hu,Shuoguo Zhang
 */

#include "bidirectional_buffer.h"
#include "density_correction.h"
#include "density_correction.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;

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
    explicit LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

struct RightInflowPressure
{
    Real outlet_pressure = 0.0;

    template <class BoundaryConditionType>
    explicit RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = outlet_pressure;
        return pressure;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_ = 51.3 * 0.01;

    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &boundary_condition) {}

    Vecd operator()(Vecd &, Vecd &, Real)
    {
        return u_ref_ * Vec2d::UnitX();
    }
};

// delete particles too close to the contact body
void delete_particles(BaseContactRelation &contact_relation)
{
    auto &body = contact_relation.getSPHBody();
    auto &particles = body.getBaseParticles();
    Real h_1 = body.getSPHAdaptation().ReferenceSmoothingLength();
    std::cout << "h_1 = " << h_1 << std::endl;

    std::cout << "total_real_particles before = " << particles.TotalRealParticles() << std::endl;

    // looping over all source body particles
    // if there is a neighbor on the contact body, store the id
    // sort all the ids to remove
    // loop backwards and switch to buffer
    IndexVector ids_to_remove;
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        bool if_remove = false;
        for (size_t k = 0; k < contact_relation.contact_bodies_.size(); k++)
        {
            auto &body_k = *contact_relation.contact_bodies_[k];
            Real h_k = body_k.getSPHAdaptation().ReferenceSmoothingLength();
            Real distance = 0.5 * (h_1 + h_k);

            Neighborhood &contact_neighborhood = (contact_relation.contact_configuration_[k])[i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                if (contact_neighborhood.r_ij_[n] <= distance)
                {
                    ids_to_remove.push_back(i);
                    if_remove = true;
                    break;
                }
            }
            if (if_remove)
                break;
        }
    }
    // sort
    sort(ids_to_remove.begin(), ids_to_remove.end());
    // remove duplicates
    ids_to_remove.erase(unique(ids_to_remove.begin(), ids_to_remove.end()), ids_to_remove.end());
    for (size_t i = ids_to_remove.size(); i-- > 0;)
        particles.switchToBufferParticle(ids_to_remove[i]);
    std::cout << "total_real_particles after = " << particles.TotalRealParticles() << std::endl;
}

//----------------------------------------------------------------------
void run_fsi2()
{
    constexpr Real mm_to_m = 1e-3;
    constexpr Real L_to_m3 = 1e-3;
    constexpr Real min_to_s = 60.0;
    constexpr Real degree_to_rad = Pi / 180.0;

    // geometry
    const Real height = 20 * mm_to_m;
    const Real length = 6 * height;
    const Real width = 117 * mm_to_m;
    const Real radius = 20 * mm_to_m;
    const Real shell_length = 26 * mm_to_m;
    const Real shell_thickness = 0.16 * mm_to_m;
    const Real angle = 45 * degree_to_rad;

    const Real dp_fluid = height / 20.0;
    const Real dp_solid = dp_fluid;
    const Real dp_shell = dp_fluid;
    std::cout << "dp_shell/thickness = " << dp_shell / shell_thickness << std::endl;

    const Real buffer_length = dp_fluid * 3.0;
    const Real wall_thickness = dp_fluid * 4.0;

    const Vec2d circle_center(2.5 * height, height);

    // material
    // fluid
    const Real rho0_f = 1000;                      /**< Density. */
    const Real Q_max = 23.44 * L_to_m3 / min_to_s; /**< Maximum flow rate. */
    const Real U_f = Q_max / (width * height);     /**< Characteristic velocity. */
    const Real U_max = 1.5 * U_f;                  /**< Maximum velocity. */
    const Real c_f = 10.0 * U_max;                 /**< Speed of sound. */
    const Real mu_f = 4.3e-3;                      /**< Dynamics viscosity. */
    const Real Re = rho0_f * U_f * height / mu_f;
    std::cout << "The Reynolds number is " << Re << std::endl;

    // solid
    const Real rho0_s = rho0_f;        /**< Reference density.*/
    const Real poisson = 0.49;         /**< Poisson ratio.*/
    const Real Youngs_modulus = 1.5e6; /**< Youngs modulus.*/

    // create shape
    // fluid shape
    std::cout << "Creating fluid shape ... " << std::endl;
    MultiPolygon fluid_shape_poly;
    Vec2d fluid_halfsize(0.5 * length, 0.5 * height);
    Vec2d translation = fluid_halfsize;
    fluid_shape_poly.addBox(Transform(translation), fluid_halfsize, GeometricOps::add);
    fluid_shape_poly.addCircle(circle_center, radius, 10, GeometricOps::add);
    MultiPolygonShape fluid_shape(fluid_shape_poly, "Fluid");
    std::cout << "Creating fluid shape - done. " << std::endl;

    // wall shape
    std::cout << "Creating wall shape ... " << std::endl;
    MultiPolygon wall_shape_poly;
    Vec2d wall_outer_halfsize(0.5 * length + 0.5 * dp_fluid, 0.5 * height + wall_thickness);
    Vec2d wall_inner_halfsize(0.5 * length + 0.5 * dp_fluid, 0.5 * height);
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

    // remove fluid and square particles too close to the shell and run relaxation
    ContactRelation fluid_to_shell(blood_body, {&shell_body});
    sph_system.initializeSystemCellLinkedLists();
    fluid_to_shell.updateConfiguration();
    delete_particles(fluid_to_shell);
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

    auto clamp_id = [&]()
    {
        IndexVector ids;
        auto *pos = shell_body.getBaseParticles().getVariableDataByName<Vec2d>("Position");
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

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(blood_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(blood_body);

    //-----------Inflow and periodic condition --------//
    Vec2d bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_length, 1.2 * height);
    Vec2d left_bidirectional_translation = 0.5 * Vec2d(buffer_length, height);
    Vec2d right_bidirectional_translation(length - 0.5 * buffer_length, 0.5 * height);
    OrientedBoxByCell left_emitter(blood_body, OrientedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    OrientedBoxByCell right_emitter(blood_body, OrientedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationPressureComplex<Inner<>, Contact<>, Contact<>>> update_fluid_density(blood_inner, blood_wall_contact, blood_shell_contact);

    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
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
    write_real_body_states.addToWrite<Real>(blood_body, "Pressure");
    write_real_body_states.addToWrite<int>(blood_body, "Indicator");
    write_real_body_states.addToWrite<int>(blood_body, "BufferIndicator");
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
    Real end_time = 1.0;
    Real output_interval = end_time / 200.0;
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

                /** FSI for viscous force. */
                viscous_force_from_fluid.exec();
                /** Update normal direction on elastic body.*/
                shell_update_normal.exec();

                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                    /** Fluid pressure relaxation */
                    pressure_relaxation.exec(dt);
                    kernel_summation.exec();
                    // left_inflow_pressure_condition.exec(dt);
                    // right_inflow_pressure_condition.exec(dt);
                    // inflow_velocity_condition.exec();
                    /** FSI for pressure force. */
                    pressure_force_from_fluid.exec();
                    /** Fluid density relaxation */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    {
                        inner_ite_dt_s = 0;
                        Real dt_s_sum = 0.0;
                        average_velocity_and_acceleration.initialize_displacement_.exec();
                        while (dt_s_sum < dt)
                        {
                            Real dt_s = SMIN(shell_computing_time_step_size.exec(), dt - dt_s_sum);
                            shell_stress_relaxation_first_half.exec(dt_s);
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
                shell_body.updateCellLinkedList();
                blood_body_complex.updateConfiguration();
                shell_blood_contact.updateConfiguration();

                boundary_indicator.exec();
                left_bidirection_buffer.tag_buffer_particles.exec();
                right_bidirection_buffer.tag_buffer_particles.exec();
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
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
