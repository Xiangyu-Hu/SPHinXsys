/**
 * @file 	mr_flow_around_cylinder.cpp
 * @author 	Xiangyu Hu
 */
#include "sphinxsys_ck.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;                               /**< Channel length. */
Real DH = 16.0;                               /**< Channel height. */
Real particle_spacing_ref = 0.4;              /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;         /**< Sponge region to impose injection. */
Vec2d circle_center(10.0, 0.5 * DH);          /**< Location of the cylinder center. */
Real radius = 1.0;                            /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                              /**< Density. */
Real U_f = 1.0;                                 /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                          /**< Speed of sound. */
Real Re = 100.0;                                /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometric parameters
//----------------------------------------------------------------------
Vecd water_body_halfsize = Vecd(0.5 * (DL_sponge + DL), 0.5 * DH);
Vecd water_body_translation = Vec2d(-DL_sponge, 0.0) + water_body_halfsize;
Vecd refinement_region_halfsize = Vec2d(0.5 * (DL_sponge + DL) + BW, 0.25 * DH);
Vec2d refinement_region_translation = Vec2d(-DL_sponge - BW, 0.25 * DH) + refinement_region_halfsize;
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(
            Transform(water_body_translation), water_body_halfsize, "OuterBoundary");
        subtract<GeometricShapeBall>(circle_center, radius);
    }
};
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    FreeStreamVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        target_velocity[0] = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineAdaptation<ParticleRefinementWithinShape>(1.3, 1.0, 1);
    water_body.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(sph_system);
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    GeometricShapeBox refinement_region(
        Transform(refinement_region_translation), refinement_region_halfsize, "RefinementRegion");
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_body.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_body.getName())
        : water_body.generateParticles<BaseParticles, Lattice, Adaptive>(refinement_region);

    SolidBody cylinder(sph_system, makeShared<GeometricShapeBall>(circle_center, radius, "Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 4.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        AdaptiveInnerRelation water_body_inner(water_body);
        AdaptiveContactRelation water_contact(water_body, {&cylinder});
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_body);
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_real_body_particle_reload_files({&water_body, &cylinder});
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(cylinder_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            DynamicsArgs(water_body_inner, std::string("OuterBoundary")), water_contact);
        SimpleDynamics<UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(water_body, refinement_region);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            update_smoothing_length_ratio.exec();
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finished !" << std::endl;
        /** Output results. */
        write_real_body_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Body partitioning for multi-resolution simulation.
    //----------------------------------------------------------------------
    BodyPartitionSpatial water_low_resolution_level(water_body, 0);
    BodyPartitionSpatial water_high_resolution_level(water_body, 1);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Relation<Inner<BodyPartitionSpatial, SmoothingLength<Continuous>>> water_low_resolution_inner(water_low_resolution_level);
    Relation<Inner<BodyPartitionSpatial, SmoothingLength<Continuous>>> water_high_resolution_inner(water_high_resolution_level);
    Relation<Contact<BodyPartitionSpatial, BodyPartitionSpatial, SmoothingLength<Continuous>>>
        water_increase_resolution_contact(water_low_resolution_level, {&water_high_resolution_level});
    Relation<Contact<BodyPartitionSpatial, BodyPartitionSpatial, SmoothingLength<Continuous>>>
        water_decrease_resolution_contact(water_high_resolution_level, {&water_low_resolution_level});
    Relation<Contact<BodyPartitionSpatial, RealBody, SmoothingLength<Continuous, SingleValued>>> water_cylinder_contact(water_high_resolution_level, {&cylinder});
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, BodyPartitionSpatial> water_low_resolution_cell_linked_list(water_low_resolution_level);
    UpdateCellLinkedList<MainExecutionPolicy, BodyPartitionSpatial> water_high_resolution_cell_linked_list(water_high_resolution_level);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> cylinder_cell_linked_list(cylinder);

    UpdateRelation<
        MainExecutionPolicy,
        Inner<BodyPartitionSpatial, SmoothingLength<Continuous>>,
        Contact<BodyPartitionSpatial, BodyPartitionSpatial, SmoothingLength<Continuous>>,
        Contact<BodyPartitionSpatial, RealBody, SmoothingLength<Continuous, SingleValued>>>
        water_high_resolution_update_complex_relation(
            water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);

    UpdateRelation<
        MainExecutionPolicy,
        Inner<BodyPartitionSpatial, SmoothingLength<Continuous>>,
        Contact<BodyPartitionSpatial, BodyPartitionSpatial, SmoothingLength<Continuous>>>
        water_low_resolution_update_complex_relation(
            water_low_resolution_inner, water_increase_resolution_contact);

    StateDynamics<execution::SequencedPolicy, AdaptLevelIndication<Refinement<Continuous, Fixed>>> water_adapt_level_indication(water_body);
    StartupAcceleration time_dependent_acceleration(Vec2d(U_f, 0.0), 2.0);
    StateDynamics<MainExecutionPolicy, GravityForceCK<StartupAcceleration>> constant_gravity(water_body, time_dependent_acceleration);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> cylinder_normal_direction(cylinder); // run on CPU
                                                                                                         /*
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::FreeSurfaceIndicationCK<
                                                                                                                     Inner<WithUpdate, Internal, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous, SingleValued>>>>
                                                                                                                 water_high_resolution_boundary_indicator(
                                                                                                                     water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::FreeSurfaceIndicationCK<
                                                                                                                     Inner<WithUpdate, Internal, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>>>
                                                                                                                 water_low_resolution_boundary_indicator(
                                                                                                                     water_low_resolution_inner, water_increase_resolution_contact);
                                                                                                     
                                                                                                             StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_body);
                                                                                                             StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_body);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 LinearCorrectionMatrix<
                                                                                                                     Inner<WithUpdate, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous, SingleValued>>>>
                                                                                                                 water_high_resolution_linear_correction_matrix(
                                                                                                                     water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 LinearCorrectionMatrix<
                                                                                                                     Inner<WithUpdate, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>>>
                                                                                                                 water_low_resolution_linear_correction_matrix(
                                                                                                                     water_low_resolution_inner, water_increase_resolution_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::DensityRegularization<
                                                                                                                     Inner<WithUpdate, Internal, AllParticles, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous, SingleValued>>>>
                                                                                                                 water_high_resolution_density_regularization(
                                                                                                                     water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::DensityRegularization<
                                                                                                                     Inner<WithUpdate, Internal, AllParticles, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<BodyPartitionSpatial, SmoothingLength<Continuous>>>>
                                                                                                                 water_low_resolution_density_regularization(
                                                                                                                     water_low_resolution_inner, water_increase_resolution_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::AcousticStep1stHalf<
                                                                                                                     Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous, SingleValued>>>>
                                                                                                                 water_high_resolution_acoustic_step_1st_half(
                                                                                                                     water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::AcousticStep2ndHalf<
                                                                                                                     Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous, SingleValued>>>>
                                                                                                                 water_high_resolution_acoustic_step_2nd_half(
                                                                                                                     water_high_resolution_inner, water_decrease_resolution_contact, water_cylinder_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::AcousticStep1stHalf<
                                                                                                                     Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>>>
                                                                                                                 water_low_resolution_acoustic_step_1st_half(
                                                                                                                     water_low_resolution_inner, water_increase_resolution_contact);
                                                                                                     
                                                                                                             InteractionDynamicsCK<
                                                                                                                 MainExecutionPolicy,
                                                                                                                 fluid_dynamics::AcousticStep2ndHalf<
                                                                                                                     Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>,
                                                                                                                     Contact<AcousticRiemannSolverCK, LinearCorrectionCK, BodyPartitionSpatial, SmoothingLength<Continuous>>>>
                                                                                                                 water_low_resolution_acoustic_step_2nd_half(
                                                                                                                     water_low_resolution_inner, water_increase_resolution_contact);
                                                                                                         */
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<int>(water_body, "AdaptLevel");
    body_states_recording.addToWrite<Real>(water_body, "SmoothingLength");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_adapt_level_indication.exec();
    water_low_resolution_cell_linked_list.exec();
    water_high_resolution_cell_linked_list.exec();
    cylinder_cell_linked_list.exec();
    water_low_resolution_update_complex_relation.exec();
    water_high_resolution_update_complex_relation.exec();

    body_states_recording.writeToFile(MainExecutionPolicy{});
    return 0;
}
