/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    Kernel *kernel_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_),
          kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
          inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
        {
            contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
        inner_relation.base_particles_.registerVariable(dW_ijV_je_ij_ttl, "TotalKernelGrad");
        inner_relation.base_particles_.registerVariable(number_of_inner_neighbor, "InnerNeighborNumber");
        inner_relation.base_particles_.registerVariable(number_of_contact_neighbor, "ContactNeighborNumber");
    }

    inline void exec()
    {
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                int N_inner_number = 0;
                int N_contact_number = 0;
                Real W_ijV_j_ttl_i = particles_->Vol_[index_i] * kernel_->W(0, ZeroVecd);
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }
};

class SmoothingNormal : public ParticleSmoothing<Vecd>
{
  public:
    explicit SmoothingNormal(BaseInnerRelation &inner_relation)
        : ParticleSmoothing<Vecd>(inner_relation, "NormalDirection"),
          n0_(*inner_relation.getSPHBody().getBaseParticles().getVariableByName<Vecd>("InitialNormalDirection")){};
    void update(size_t index_i, Real dt = 0.0)
    {
        ParticleSmoothing<Vecd>::update(index_i, dt);
        smoothed_[index_i] = temp_[index_i] / (temp_[index_i].norm() + TinyReal);
        // smoothed_[index_i] /= temp_[index_i].norm() + TinyReal;
        smoothed_[index_i] = smoothed_[index_i].normalized();
        n0_[index_i] = smoothed_[index_i];
    };

  private:
    StdLargeVec<Vecd> n0_;
};
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;                 /**< Reference length. */
Real DH = 3.0;                 /**< Reference and the height of main channel. */
Real DL1 = 0.7 * DL;           /**< The length of the main channel. */
Real DH2 = DL - DL1;           /**< Height of side channel. */
Real DL2 = 1.0 * DH;           /**< Length of side channel. */
Real angle = 7.5 / 180.0 * Pi; /**< Angle between side channel and horizontal line. */
Real cos_a = cos(angle);
Real sin_a = sin(angle);
Real DL3 = 0.5 * DH / sin_a + DL2 - DH2 * cos_a / sin_a; /**< Length of side channel on the other side. */

Real resolution_ref = 0.15;                  /**< Initial reference particle spacing. */
Real resolution_wall = resolution_ref * 0.5; /**< Initial reference particle spacing of wall. */
Real BW = resolution_ref * 4;                /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;        /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;    /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * DH2));
Real Re = 100.0;                    /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the water block in T shape polygon. */
Vec2d point_water_1(-DL_sponge, -0.5 * DH);
Vec2d point_water_2(-DL_sponge, 0.5 * DH);
Vec2d point_water_3(DL1, 0.5 * DH);
Vec2d point_water_4 = point_water_3 + Vec2d(DL2 * cos_a, DL2 *sin_a);
Vec2d point_water_5 = point_water_4 + Vec2d(DH2 * sin_a, -DH2 *cos_a);
Vec2d point_water_6(point_water_5.x() - (0.5 * DH / sin_a + DL2 - DH2 * 1 / tan(angle)) * cos_a, 0);
Vec2d point_water_7(point_water_5.x(), -point_water_5.y());
Vec2d point_water_8(point_water_4.x(), -point_water_4.y());
Vec2d point_water_9(point_water_3.x(), -point_water_3.y());
std::vector<Vecd> water_block_shape{
    point_water_1, point_water_2, point_water_3, point_water_4,
    point_water_5, point_water_6, point_water_7, point_water_8,
    point_water_9, point_water_1};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
/** Particle generator for the wall boundary. */
class WallBoundaryParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit WallBoundaryParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        const Vec2d n_main_channel_1(0, -1); // normal of main channel
        const Vec2d n_main_channel_2 = -n_main_channel_1;
        const Vec2d n_side_channel_1(sin_a, -cos_a); // normal of side channel
        const Vec2d n_side_channel_2(sin_a, cos_a);
        const Vec2d n_main_corner_1 = (n_main_channel_1 + n_side_channel_1).normalized(); //  normal of corner between main and side channel
        const Vec2d n_main_corner_2 = (n_main_channel_2 + n_side_channel_2).normalized();
        const Vec2d n_side_corner(-1, 0); // normal of corner between two side channels

        // main channel corner
        const Real x_main_corner = (point_water_3 + 0.5 * Vec2d(-resolution_wall * tan(0.5 * angle), resolution_wall)).x();
        const Real y_main_corner = 0.5 * DH + 0.5 * resolution_wall;

        initializePositionAndVolumetricMeasure(Vecd(x_main_corner, y_main_corner), resolution_wall);
        initializeSurfaceProperties(n_main_corner_1, resolution_wall);

        initializePositionAndVolumetricMeasure(Vecd(x_main_corner, -y_main_corner), resolution_wall);
        initializeSurfaceProperties(n_main_corner_2, resolution_wall);

        // main channel
        {
            int particle_number_mid_surface = int((DL1 + DL_sponge + BW) / resolution_wall);
            for (int i = 1; i < particle_number_mid_surface; i++)
            {
                Real x = x_main_corner - double(i) * resolution_wall;
                Real y = y_main_corner;
                initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_wall);
                initializeSurfaceProperties(n_main_channel_1, resolution_wall);

                initializePositionAndVolumetricMeasure(Vecd(x, -y), resolution_wall);
                initializeSurfaceProperties(n_main_channel_2, resolution_wall);
            }
        }

        // left side channel
        {
            int particle_number_mid_surface = int((DL2 + BW) / resolution_wall);
            for (int i = 1; i < particle_number_mid_surface; i++)
            {
                Real x = x_main_corner + double(i) * resolution_wall * cos_a;
                Real y = y_main_corner + double(i) * resolution_wall * sin_a;
                initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_wall);
                initializeSurfaceProperties(n_side_channel_1, resolution_wall);

                initializePositionAndVolumetricMeasure(Vecd(x, -y), resolution_wall);
                initializeSurfaceProperties(n_side_channel_2, resolution_wall);
            }
        }

        // side channel corner
        const Real x_side_corner = point_water_6.x() + 0.5 * resolution_wall / sin_a;
        const Real y_side_corner = 0.0;

        initializePositionAndVolumetricMeasure(Vecd(x_side_corner, y_side_corner), resolution_wall);
        initializeSurfaceProperties(n_side_corner, resolution_wall);

        // right side channel
        {
            int particle_number_mid_surface = int((DL3 + BW) / resolution_wall);
            for (int i = 1; i < particle_number_mid_surface; i++)
            {
                Real x = x_side_corner + double(i) * resolution_wall * cos_a;
                Real y = y_side_corner + double(i) * resolution_wall * sin_a;
                initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_wall);
                initializeSurfaceProperties(-n_side_channel_1, resolution_wall);

                initializePositionAndVolumetricMeasure(Vecd(x, -y), resolution_wall);
                initializeSurfaceProperties(-n_side_channel_2, resolution_wall);
            }
        }
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        return target_velocity;
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5 * DH - DL2),
                                     Vec2d(DL1 + DL3, 0.5 * DH + DL2));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_wall);
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0);
    wall_boundary.generateParticles<WallBoundaryParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_boundary_inner(wall_boundary);
    ShellInnerRelationWithContactKernel wall_curvature_inner(wall_boundary, water_block);
    ContactRelationToShell water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    water_block.addBodyStateForRecording<Real>("Pressure"); // output for debug
    water_block.addBodyStateForRecording<int>("Indicator"); // output for debug

    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> wall_corrected_configuration(wall_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(wall_curvature_inner);
    InteractionWithUpdate<SmoothingNormal> smoothing_normal(shell_boundary_inner);

    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d emitter_translation = Vec2d(-DL_sponge + 0.5 * BW, 0.0);
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);

    Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d inlet_buffer_translation = Vec2d(-0.5 * DL_sponge, 0.0);
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    Vec2d disposer_up_halfsize = Vec2d(0.5 * BW, 0.3 * DH2);
    Vec2d disposer_up_translation = 0.5 * (point_water_4 + point_water_5);
    BodyAlignedBoxByCell disposer_up(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(angle), Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up, xAxis);

    Vec2d disposer_down_halfsize = disposer_up_halfsize;
    Vec2d disposer_down_translation = 0.5 * (point_water_7 + point_water_8);
    BodyAlignedBoxByCell disposer_down(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(-angle), Vec2d(disposer_down_translation)), disposer_down_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_down_outflow_deletion(disposer_down, xAxis);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    wall_boundary.addBodyStateForRecording<Real>("TotalMeanCurvature");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    // for (int i = 0; i < 2; i++)
    //     smoothing_normal.exec();
    // wall_corrected_configuration.exec();
    shell_curvature.compute_initial_curvature();
    water_block_complex.updateConfiguration();

    // Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness(water_block_inner,
                                                      water_wall_contact);
    check_kernel_completeness.exec();
    water_block.addBodyStateForRecording<Real>("TotalKernel");
    water_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 100.0;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                pressure_relaxation.exec(dt);
                emitter_buffer_inflow_condition.exec();
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_up_outflow_deletion.exec();
            disposer_down_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        check_kernel_completeness.exec();
        write_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
