/**
 * @file fsi2.cpp
 * @brief This is the benchmark test of fluid-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 * The case can be found in Chi Zhang, Massoud Rezavand, Xiangyu Hu,
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * Journal of Computation Physics 404 (2020) 109135.
 * @author Chi Zhang and Xiangyu Hu
 */
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
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

void run_relaxation(RealBody &body)
{
    InnerRelation inner(body);
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

namespace SPH
{
class SlipWall;

namespace fluid_dynamics
{
template <class RiemannSolverType>
class Integration2ndHalf<Contact<SlipWall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit Integration2ndHalf(BaseContactRelation &wall_contact_relation)
        : BaseIntegrationWithWall(wall_contact_relation),
          riemann_solver_(this->fluid_, this->fluid_) {}
    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real density_change_rate = 0.0;
        Vecd p_dissipation = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const auto *vel_ave_k = wall_vel_ave_[k];
            const auto *n_k = wall_n_[k];
            const auto *wall_Vol_k = wall_Vol_[k];
            Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                size_t index_j = wall_neighborhood.j_[n];
                Vecd &e_ij = wall_neighborhood.e_ij_[n];
                Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];

                Vecd face_to_fluid_n = SGN(e_ij.dot(n_k[index_j])) * n_k[index_j];

                Vecd vel_j_n_in_wall = vel_[index_i].dot(face_to_fluid_n) * face_to_fluid_n;
                Vecd vel_j_t_in_wall = vel_[index_i] - vel_j_n_in_wall;
                Vecd vel_j_in_wall = vel_j_t_in_wall - vel_j_n_in_wall;

                density_change_rate += (vel_[index_i] - vel_j_in_wall).dot(e_ij) * dW_ijV_j;
                Real u_jump = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(face_to_fluid_n);
                p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * face_to_fluid_n;
            }
        }
        drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
        force_[index_i] += p_dissipation * this->Vol_[index_i];
    }

  private:
    RiemannSolverType riemann_solver_;
};
} // namespace fluid_dynamics
} // namespace SPH

void run_fsi2()
{
    constexpr Real cm_to_m = 0.01;
    constexpr Real g_to_kg = 0.001;

    // geometry
    const Real channel_length_1 = 5.5 * cm_to_m; // from inlet to the right edge of the square
    const Real channel_length_2 = 15.5 * cm_to_m;
    const Real channel_length = channel_length_1 + channel_length_2;
    const Real channel_height = 12.0 * cm_to_m;
    const Real square_length = 1.0 * cm_to_m;
    const Real beam_length = 4.0 * cm_to_m;
    const Real beam_thickness = 0.06 * cm_to_m;

    const Real dp_solid = square_length / 10.0;
    const Real dp_fluid = 2 * dp_solid;
    std::cout << "dp_solid/thickness = " << dp_solid / beam_thickness << std::endl;

    const Real buffer_length = dp_fluid * 3.0;
    const Real channel_total_length = channel_length + 2 * buffer_length;
    const Real wall_thickness = dp_fluid * 4.0;

    const Vec2d square_center(channel_length_1 - 0.5 * square_length, 0.5 * channel_height);

    // material
    // fluid
    const Real rho0_f = 1.18e-3 * g_to_kg / std::pow(cm_to_m, 3); /**< Density. */
    const Real U_f = 51.3 * cm_to_m;                              /**< Characteristic velocity. */
    const Real c_f = 20.0 * U_f;                                  /**< Speed of sound. */
    const Real mu_f = 1.82e-4 * g_to_kg / cm_to_m;                /**< Dynamics viscosity. */
    const Real Re = rho0_f * U_f * square_length / mu_f;
    std::cout << "The Reynolds number is " << Re << std::endl;

    // solid
    const Real rho0_s = 0.1 * g_to_kg / std::pow(cm_to_m, 3); /**< Reference density.*/
    const Real poisson = 0.35;                                /**< Poisson ratio.*/
    const Real Youngs_modulus = 2.5e6 * g_to_kg / cm_to_m;    /**< Youngs modulus.*/

    // create shape
    // fluid shape
    std::cout << "Creating fluid shape ... " << std::endl;
    MultiPolygon fluid_shape_poly;
    fluid_shape_poly.addABox(Transform(Vec2d(0.5 * channel_total_length - buffer_length, 0.5 * channel_height)),
                             0.5 * Vec2d(channel_total_length, channel_height),
                             ShapeBooleanOps::add);
    fluid_shape_poly.addABox(Transform(Vec2d(square_center)), 0.5 * square_length * Vec2d::Ones(), ShapeBooleanOps::sub);
    MultiPolygonShape fluid_shape(fluid_shape_poly, "Fluid");
    std::cout << "Creating fluid shape - done. " << std::endl;

    // wall shape
    std::cout << "Creating wall shape ... " << std::endl;
    MultiPolygon wall_shape_poly;
    wall_shape_poly.addABox(Transform(Vec2d(0.5 * channel_total_length - buffer_length, 0.5 * channel_height)),
                            Vec2d(0.5 * channel_total_length, 0.5 * channel_height + wall_thickness),
                            ShapeBooleanOps::add);
    wall_shape_poly.addABox(Transform(Vec2d(0.5 * channel_total_length - buffer_length, 0.5 * channel_height)),
                            Vec2d(0.5 * channel_total_length, 0.5 * channel_height),
                            ShapeBooleanOps::sub);
    MultiPolygonShape wall_shape(wall_shape_poly, "Wall");
    std::cout << "Creating wall shape - done. " << std::endl;

    // square
    std::cout << "Creating square shape ... " << std::endl;
    MultiPolygon square_shape_poly;
    square_shape_poly.addABox(Transform(Vec2d(square_center)), 0.5 * square_length * Vec2d::Ones(), ShapeBooleanOps::add);
    MultiPolygonShape square_shape(square_shape_poly, "Square");
    std::cout << "Creating square shape - done. " << std::endl;

    // beam
    StdVec<Vec2d> positions_beam;
    {
        Real x = square_center.x() + 0.5 * square_length + 0.5 * dp_solid;
        while (x < square_center.x() + 0.5 * square_length + beam_length)
        {
            positions_beam.emplace_back(x, square_center.y());
            x += dp_solid;
        }
    }
    StdVec<Vec2d> normals_beam(positions_beam.size(), Vec2d::UnitY());
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    const auto bbox = wall_shape.getBounds();
    SPHSystem sph_system(bbox, dp_fluid);
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, fluid_shape);
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody square_body(sph_system, square_shape);
    square_body.defineAdaptationRatios(1.15, dp_fluid / dp_solid);
    square_body.defineBodyLevelSetShape()->writeLevelSet();
    square_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? square_body.generateParticles<BaseParticles, Reload>(square_body.getName())
        : square_body.generateParticles<BaseParticles, Lattice>();

    SolidBody beam(sph_system, makeShared<DefaultShape>("Beam"));
    beam.defineAdaptationRatios(1.15, dp_fluid / dp_solid);
    beam.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam.generateParticles<SurfaceParticles, Shell>(positions_beam, normals_beam, dp_solid, beam_thickness);

    // reset mass
    SimpleDynamics<ShellFluidMixtureMass> reset_shell_mass(beam, rho0_f);
    // reset_shell_mass.exec();

    // remove fluid and square particles too close to the beam and run relaxation
    ContactRelation fluid_to_beam(water_block, {&beam});
    sph_system.initializeSystemCellLinkedLists();
    fluid_to_beam.updateConfiguration();
    delete_particles(fluid_to_beam);
    // run relaxation
    if (sph_system.RunParticleRelaxation())
        run_relaxation(square_body);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation beam_inner(beam);
    // Upper and lower walls are slip walls
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation water_square_contact(water_block, {&square_body});
    ContactRelationFSI2 water_beam_contact(water_block, {&beam});
    ContactRelationSFI2 beam_water_contact(beam, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_wall_contact, &water_square_contact, &water_beam_contact});
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
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(square_body);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> beam_corrected_configuration(beam_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> beam_stress_relaxation_first_half(beam_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> beam_stress_relaxation_second_half(beam_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> beam_computing_time_step_size(beam);

    auto clamp_id = [&]()
    {
        IndexVector ids;
        auto *pos = beam.getBaseParticles().getVariableDataByName<Vec2d>("Position");
        for (size_t i = 0; i < beam.getBaseParticles().TotalRealParticles(); ++i)
        {
            if (pos[i].x() < square_center.x() + 0.5 * square_length + 1.2 * dp_solid)
                ids.push_back(i);
        }
        return ids;
    }();
    BodyPartByParticle beam_base(beam);
    beam_base.body_part_particles_ = clamp_id;
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_beam_base(beam_base);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<ComplexInteraction<NablaWV<Inner<>, Contact<>, Contact<>, Contact<>>>> kernel_summation(water_block_inner, water_wall_contact, water_square_contact, water_beam_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>, Contact<>>>> boundary_indicator(water_block_inner, water_wall_contact, water_square_contact, water_beam_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(water_block_inner, water_wall_contact, water_square_contact, water_beam_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<SlipWall>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact, water_square_contact, water_beam_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SPHAdaptation, NoLimiter>, Contact<Boundary>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_wall_contact, water_square_contact, water_beam_contact);
    // No viscosity at the upper and lower wall boundaries
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity, NoKernelCorrection>> viscous_force(water_block_inner, water_square_contact, water_beam_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    //-----------Inflow and periodic condition --------//
    Vec2d bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_length, 1.2 * channel_height);
    Vec2d left_bidirectional_translation = 0.5 * Vec2d(-buffer_length, channel_height);
    Vec2d right_bidirectional_translation = Vec2d(channel_length + 0.5 * buffer_length, 0.5 * channel_height);
    AlignedBoxByCell left_emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    AlignedBoxByCell right_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationPressureComplex<Inner<>, Contact<>, Contact<>, Contact<>>> update_fluid_density(water_block_inner, water_wall_contact, water_square_contact, water_beam_contact);

    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(beam);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> beam_update_normal(beam);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(beam_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(beam_water_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vec2d>(water_block, "Velocity");
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    write_real_body_states.addToWrite<int>(water_block, "BufferIndicator");
    write_real_body_states.addDerivedVariableRecording<SimpleDynamics<Displacement>>(beam);
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
    /** computing surface normal direction for the insert body. */
    cylinder_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    beam_corrected_configuration.exec();
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

            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        std::cout << "Fluid relaxation finished !" << std::endl;
    };
    run_fluid_relaxation();
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 4.0;
    Real fsi_start_time = 0.5;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
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
            if (physical_time > fsi_start_time)
                viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            beam_update_normal.exec();

            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                /** FSI for pressure force. */
                if (physical_time > fsi_start_time)
                    pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                if (physical_time > fsi_start_time)
                {
                    inner_ite_dt_s = 0;
                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        Real dt_s = SMIN(beam_computing_time_step_size.exec(), dt - dt_s_sum);
                        beam_stress_relaxation_first_half.exec(dt_s);
                        constraint_beam_base.exec();
                        beam_stress_relaxation_second_half.exec(dt_s);
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

            // if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            // {
            //     particle_sorting.exec();
            // }
            water_block.updateCellLinkedList();
            if (physical_time > fsi_start_time)
                beam.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            if (physical_time > fsi_start_time)
                beam_water_contact.updateConfiguration();

            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        beam.setNewlyUpdated();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
}
