/**
 * @file 	wedge_impact_water_surface.cpp
 * @brief 	high speed impact of aluminum wedge on water surface.
 * @details Test for fluid-shell interaction
 * @author 	Weiyi Kong
 * @version 0.1
 */
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real DL_half = 1.5;                                     /**< Tank half length. */
const Real DH = 1.0;                                          /**< Tank height. */
const Real wedge_length = 0.6;                                /**< Length of the wedge. */
const Real wedge_thickness = 0.04;                            /**< Thickness of the wedge. */
const Real particle_spacing_shell = wedge_thickness / 2.0;    /**< Initial reference particle spacing*/
const Real particle_spacing_ref = 2 * particle_spacing_shell; /**< Initial reference particle spacing*/
const Real BW = particle_spacing_ref * 4.0;                   /**< Wall boundary thickness*/
const Real initial_gap = wedge_thickness;
const Real angle = 10.0 * Pi / 180.0; /**< Wedge angle*/

const BoundingBox system_domain_bounds(Vec2d(-DL_half - BW, -BW), Vec2d(DL_half + BW, DH + initial_gap + wedge_length * tan(angle)));
//----------------------------------------------------------------------
//	Define observer initial position.
//----------------------------------------------------------------------
// observer location
const Vec2d point_A(0, DH + initial_gap + 0.5 * particle_spacing_shell);
const Vec2d point_C = 0.3 * Vec2d(1, tan(angle)) + Vec2d(0, DH + initial_gap);
const Vec2d point_D = (0.6 - 0.15) * Vec2d(1, tan(angle)) + Vec2d(0, DH + initial_gap + 0.5 * particle_spacing_shell);
const StdVec<Vec2d> observation_location_pressure_A = {point_A};
const StdVec<Vec2d> observation_location_pressure_D = {point_D};
const StdVec<Vec2d> observation_location_disp = {point_C};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
const Real rho0_f = 1000.0; /**< Reference density of fluid. */
const Real U_ref = 30;
const Real U_max = U_ref;      /**< Characteristic velocity. */
const Real c_f = 10.0 * U_max; /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Material properties of the elastic wedge.
//----------------------------------------------------------------------
const Real rho0_s = 2700.0; /**< Reference solid density. */
const Real poisson = 0.34;  /**< Poisson ratio. */
const Real Ae = 6.75e10;    /**< Normalized Youngs Modulus. */
const Real Youngs_modulus = Ae;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * wedge_thickness * wedge_thickness;
//----------------------------------------------------------------------
//	Geometry definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d water_halfsize(DL_half, DH * 0.5);
        Vec2d water_transform(0, 0.5 * DH);
        multi_polygon_.addABox(Transform(water_transform), water_halfsize, ShapeBooleanOps::add);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d wall_outer_halfsize(DL_half + BW, DH * 0.5 + BW);
        Vec2d wall_outer_transform(0, 0.5 * DH);
        Vec2d wall_inner_halfsize(DL_half, 0.5 * (DH + BW));
        Vec2d wall_inner_transform(0, 0.5 * (DH + BW));
        multi_polygon_.addABox(Transform(wall_outer_transform), wall_outer_halfsize, ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(wall_inner_transform), wall_inner_halfsize, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	wall body shape definition.
//----------------------------------------------------------------------
class WedgeParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit WedgeParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        int particle_number = int(wedge_length / cos(angle) / particle_spacing_shell) + 1;
        Real x = 0;
        Real y = DH + initial_gap;
        Real dx = particle_spacing_shell * cos(angle);
        Real dy = particle_spacing_shell * sin(angle);
        initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_shell);
        initializeSurfaceProperties(Vec2d(0, 1), wedge_thickness);
        Vec2d n1(-sin(angle), cos(angle));
        Vec2d n2(sin(angle), cos(angle));
        for (int i = 1; i < particle_number; i++)
        {
            x += dx;
            y += dy;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_shell);
            initializeSurfaceProperties(n1, wedge_thickness);
            initializePositionAndVolumetricMeasure(Vecd(-x, y), particle_spacing_shell);
            initializeSurfaceProperties(n2, wedge_thickness);
        }
    }
};
//----------------------------------------------------------------------
//	create wedge constrain shape
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (
            // std::abs(base_particles_.pos_[index_i][0]) < 0.7 * particle_spacing_shell * cos(angle) ||
            base_particles_.pos_[index_i][0] > wedge_length - particle_spacing_shell * cos(angle) ||
            base_particles_.pos_[index_i][0] < -wedge_length + particle_spacing_shell * cos(angle))
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class WedgeVelocity : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit WedgeVelocity(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part){};

  protected:
    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = Vec2d(0, -U_ref);
        angular_vel_[index_i] = Vec2d(0, 0);
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody wedge(sph_system, makeShared<DefaultShape>("Wedge"));
    wedge.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_shell);
    wedge.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wedge.generateParticles<WedgeParticleGenerator>();
    //----------------------------------------------------------------------
    //	Particle and body creation of wedge observer.
    //----------------------------------------------------------------------
    ObserverBody disp_observer(sph_system, "DisplacementObserver");
    disp_observer.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_shell);
    disp_observer.generateParticles<ParticleGeneratorObserver>(observation_location_disp);

    ObserverBody pressure_observer_A(sph_system, "PressureObserverA");
    disp_observer.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_shell);
    pressure_observer_A.generateParticles<ParticleGeneratorObserver>(observation_location_pressure_A);

    ObserverBody pressure_observer_D(sph_system, "PressureObserverD");
    disp_observer.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_shell);
    pressure_observer_D.generateParticles<ParticleGeneratorObserver>(observation_location_pressure_D);

    auto get_observer_closest_id = [&](const Vecd &observer_pos)
    {
        const auto &pos = wedge.getBaseParticles().pos_;
        Real distance2 = MaxReal;
        size_t id = 0;
        for (size_t index_i = 0; index_i < wedge.getBaseParticles().total_real_particles_; index_i++)
        {
            Vec2d disp = pos[index_i] - observer_pos;
            Real d_i2 = disp.x() * disp.x() + disp.y() * disp.y();

            if (d_i2 < distance2)
            {
                id = index_i;
                distance2 = d_i2;
            }
        }
        return id;
    };

    size_t observer_D_id = get_observer_closest_id(point_D);
    size_t observer_A_id = get_observer_closest_id(point_A);

    auto update_observer_pos = [&]()
    {
        pressure_observer_A.getBaseParticles().pos_[0] = wedge.getBaseParticles().pos_[observer_A_id] + Vecd(0, -particle_spacing_shell);
        pressure_observer_D.getBaseParticles().pos_[0] = wedge.getBaseParticles().pos_[observer_D_id] + Vecd(0, -particle_spacing_shell);
    };
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation wedge_inner(wedge);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelationToShell water_wedge_contact(water_block, {&wedge}, false);
    ContactRelationFromShell wedge_water_contact(wedge, {&water_block}, false);
    ContactRelation disp_observer_contact(disp_observer, {&wedge});
    ContactRelation pressure_observer_A_contact(pressure_observer_A, {&water_block});
    ContactRelation pressure_observer_D_contact(pressure_observer_D, {&water_block});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(wedge, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_wall_contact, &water_wedge_contact});
    //----------------------------------------------------------------------
    //	Define fluid methods which are used in this case.
    //----------------------------------------------------------------------
    /** Evaluation of fluid density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>, Contact<>>> update_fluid_density(water_block_inner, water_wall_contact, water_wedge_contact);
    /** Compute time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Compute time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(water_block_inner, water_wall_contact, water_wedge_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact, water_wedge_contact);
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> wedge_corrected_configuration(wedge_inner);
    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> wedge_computing_time_step_size(wedge);
    /** Stress relaxation stepping for the elastic wedge. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> wedge_stress_relaxation_first_half(wedge_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> wedge_stress_relaxation_second_half(wedge_inner);
    /**Constrain a solid body part.  */
    BoundaryGeometry wedge_constraint_part(wedge, "BcGeometry");
    SimpleDynamics<WedgeVelocity> wedge_velocity_constraint(wedge_constraint_part);
    /** Update the norm of elastic wedge. */
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> wedge_update_normal(wedge);
    /** Curvature calculation for elastic wedge. */
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> wedge_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Define fsi methods which are used in this case.
    //----------------------------------------------------------------------
    /** Compute the force exerted on elastic wedge due to fluid pressure. */
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid> fluid_pressure_force_on_wedge(wedge_water_contact);
    auto update_average_velocity = [&]()
    {
        auto avg_vel = *wedge.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
        auto avg_force = *wedge.getBaseParticles().getVariableByName<Vecd>("AverageForce");
        particle_for(
            par,
            wedge.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                avg_vel[index_i] = wedge.getBaseParticles().vel_[index_i];
                avg_force[index_i] = wedge.getBaseParticles().force_[index_i];
            });
    };
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    wedge.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    wedge.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    /** Output body states for visualization. */
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system.real_bodies_);
    /** Output the observed displacement of wedge center. */
    ObservedQuantityRecording<Vecd> write_displacement("Displacement", disp_observer_contact);
    ObservedQuantityRecording<Real> write_pressure_A("Pressure", pressure_observer_A_contact);
    ObservedQuantityRecording<Real> write_pressure_D("Pressure", pressure_observer_D_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing linear reproducing configuration for the insert body. */
    wedge_corrected_configuration.exec();
    wedge_curvature.exec();
    /** update fluid-shell contact*/
    water_wedge_contact.updateConfiguration();
    wedge_water_contact.updateConfiguration();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile(0);
    write_displacement.writeToFile(0);
    write_pressure_A.writeToFile(0);
    write_pressure_D.writeToFile(0);
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    const Real Dt_ref = get_fluid_advection_time_step_size.exec();
    const Real dt_ref = get_fluid_time_step_size.exec();
    const Real dt_s_ref = wedge_computing_time_step_size.exec();
    std::cout << "Dt_ref = " << Dt_ref << std::endl;
    std::cout << "dt_ref = " << dt_ref << std::endl;
    std::cout << "dt_s_ref = " << dt_s_ref << std::endl;

    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = initial_gap / U_ref + 2.5e-3; /**< End time. */
    Real output_interval = end_time / 200.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    // initial velocity
    {
        particle_for(
            par,
            wedge.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                wedge.getBaseParticles().vel_[index_i][1] = -U_ref;
            });
    }
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            Real dt_temp = get_fluid_time_step_size.exec();
            Real dt_s_temp = wedge_computing_time_step_size.exec();
            Real dt = std::min({dt_s_temp, dt_temp, Dt});

            update_fluid_density.exec();

            /** Fluid relaxation and force computation. */
            pressure_relaxation.exec(dt);
            fluid_pressure_force_on_wedge.exec();
            density_relaxation.exec(dt);
            /** Solid dynamics time stepping. */
            wedge_stress_relaxation_first_half.exec(dt);
            wedge_velocity_constraint.exec();
            wedge_stress_relaxation_second_half.exec(dt);
            update_average_velocity();

            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** Update normal direction on elastic body. */
            wedge_update_normal.exec();
            /** Update curvature. */
            wedge.updateCellLinkedList();
            shell_curvature_inner.updateConfiguration();
            wedge_curvature.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList(); // water particle motion is small

            /** one need update configuration after periodic condition. */
            water_block_complex.updateConfiguration();
            wedge_water_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();

        /** Output the observed data. */
        write_real_body_states_to_vtp.writeToFile();
        write_displacement.writeToFile(number_of_iterations);
        /** update pressure observer position*/
        update_observer_pos();
        pressure_observer_A_contact.updateConfiguration();
        pressure_observer_D_contact.updateConfiguration();
        write_pressure_A.writeToFile(number_of_iterations);
        write_pressure_D.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
