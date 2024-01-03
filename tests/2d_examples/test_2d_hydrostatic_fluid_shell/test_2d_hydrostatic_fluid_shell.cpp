/**
 * @file 	hydrostatic_fsi.cpp
 * @brief 	structure deformation due to hydrostatic pressure under gravity.
 * @details This is the one of the basic test cases
 * for understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Yujie Zhu, Chi Zhang and Xiangyu Hu
 * @version 0.1
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_), inner_configuration_(&inner_relation.inner_configuration_)
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
                Real W_ijV_j_ttl_i = 0;
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

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real DL = 1.0;                                     /**< Tank length. */
const Real DH = 2.1;                                     /**< Tank height. */
const Real Dam_L = 1.0;                                  /**< Water block width. */
const Real Dam_H = 2.0;                                  /**< Water block height. */
const Real Gate_thickness = 0.05;                        /**< Width of the gate. */
const Real particle_spacing_gate = Gate_thickness / 4.0; /**< Initial reference particle spacing*/
const Real particle_spacing_ref = particle_spacing_gate; /**< Initial reference particle spacing*/
const Real BW = particle_spacing_ref * 4.0;
const BoundingBox system_domain_bounds(Vec2d(-BW, -std::max(particle_spacing_gate, Gate_thickness)), Vec2d(DL + BW, DH + Gate_thickness));
//----------------------------------------------------------------------
//	Define the corner point of water block geometry.
//----------------------------------------------------------------------
const Vec2d DamP_lb(0.0, 0.0);     /**< Left bottom. */
const Vec2d DamP_lt(0.0, Dam_H);   /**< Left top. */
const Vec2d DamP_rt(Dam_L, Dam_H); /**< Right top. */
const Vec2d DamP_rb(Dam_L, 0.0);   /**< Right bottom. */
//----------------------------------------------------------------------
//	Define the geometry for gate constrain.
//----------------------------------------------------------------------
Vec2d ConstrainLP_lb(-BW, -particle_spacing_ref);
Vec2d ConstrainLP_lt(-BW, 0.0);
Vec2d ConstrainLP_rt(0.0, 0.0);
Vec2d ConstrainLP_rb(0.0, -particle_spacing_ref);
Vec2d ConstrainRP_lb(Dam_L, -particle_spacing_ref);
Vec2d ConstrainRP_lt(Dam_L, 0.0);
Vec2d ConstrainRP_rt(Dam_L + BW, 0.0);
Vec2d ConstrainRP_rb(Dam_L + BW, -particle_spacing_ref);
// observer location
const StdVec<Vecd> observation_location = {Vecd(0.5 * Dam_L, -0.5 * particle_spacing_gate)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
const Real rho0_f = 1000.0;                       /**< Reference density of fluid. */
const Real gravity_g = 9.81;                      /**< Value of gravity. */
const Real U_ref = 2.0 * sqrt(Dam_H * gravity_g); /**< Characteristic velocity. */
const Real c_f = 10.0 * U_ref;                    /**< Reference sound speed. */
const Real Re = 0.1;                              /**< Reynolds number. */
const Real mu_f = rho0_f * U_ref * DL / Re;       /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Material properties of the elastic gate.
//----------------------------------------------------------------------
const Real rho0_s = 2700.0; /**< Reference solid density. */
const Real poisson = 0.49;  /**< Poisson ratio. */
const Real Ae = 6.75e10;    /**< Normalized Youngs Modulus. */
const Real Youngs_modulus = Ae;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * Gate_thickness;
//----------------------------------------------------------------------
//	Geometry definition.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(DamP_lb);
    water_block_shape.push_back(DamP_lt);
    water_block_shape.push_back(DamP_rt);
    water_block_shape.push_back(DamP_rb);
    water_block_shape.push_back(DamP_lb);

    return water_block_shape;
}
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	wall body shape definition.
//----------------------------------------------------------------------
const int particle_number_wall = int(DH / particle_spacing_gate);
class WallBoundaryParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit WallBoundaryParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_wall; i++)
        {
            Real x1 = -0.5 * particle_spacing_gate;
            Real x2 = DL + 0.5 * particle_spacing_gate;
            Real y = 0.5 * particle_spacing_gate + Real(i) * particle_spacing_gate;
            const Vec2d normal_direction_1(1.0, 0.0);
            const Vec2d normal_direction_2(-1.0, 0.0);
            initializePositionAndVolumetricMeasure(Vecd(x1, y), particle_spacing_gate);
            initializeSurfaceProperties(normal_direction_1, Gate_thickness);
            initializePositionAndVolumetricMeasure(Vecd(x2, y), particle_spacing_gate);
            initializeSurfaceProperties(normal_direction_2, Gate_thickness);
        }
    }
};
//----------------------------------------------------------------------
//	wall body shape definition.
//----------------------------------------------------------------------
const int particle_number_gate = int((DL + 2 * BW) / particle_spacing_gate);
class GateParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit GateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        // generate particles for the elastic gate
        for (int i = 0; i < particle_number_gate; i++)
        {
            Real x = -BW + 0.5 * particle_spacing_gate + Real(i) * particle_spacing_gate;
            Real y = -0.5 * particle_spacing_gate;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_gate);
            const Vec2d normal_direction(0, 1.0);
            initializeSurfaceProperties(normal_direction, Gate_thickness);
        }
    }
};
//----------------------------------------------------------------------
//	create Gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape_left;
    gate_constraint_shape_left.push_back(ConstrainLP_lb);
    gate_constraint_shape_left.push_back(ConstrainLP_lt);
    gate_constraint_shape_left.push_back(ConstrainLP_rt);
    gate_constraint_shape_left.push_back(ConstrainLP_rb);
    gate_constraint_shape_left.push_back(ConstrainLP_lb);

    std::vector<Vecd> gate_constraint_shape_right;
    gate_constraint_shape_right.push_back(ConstrainRP_lb);
    gate_constraint_shape_right.push_back(ConstrainRP_lt);
    gate_constraint_shape_right.push_back(ConstrainRP_rt);
    gate_constraint_shape_right.push_back(ConstrainRP_rb);
    gate_constraint_shape_right.push_back(ConstrainRP_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(gate_constraint_shape_left, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(gate_constraint_shape_right, ShapeBooleanOps::add);
    return multi_polygon;
}
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
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_gate);
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1, 1, 0);
    wall_boundary.generateParticles<WallBoundaryParticleGenerator>();

    SolidBody gate(sph_system, makeShared<DefaultShape>("Gate"));
    gate.defineAdaptation<SPHAdaptation>(1.15, particle_spacing_ref / particle_spacing_gate);
    gate.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<GateParticleGenerator>();
    //----------------------------------------------------------------------
    //	Particle and body creation of gate observer.
    //----------------------------------------------------------------------
    ObserverBody gate_observer(sph_system, "Observer");
    gate_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation gate_inner(gate);
    ContactRelationToShell water_block_contact(water_block, {&wall_boundary, &gate});
    ContactRelationFromShell gate_contact(gate, {&water_block});
    ContactRelation gate_observer_contact(gate_observer, {&gate});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(gate, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define fluid methods which are used in this case.
    //----------------------------------------------------------------------
    /** Initialize particle acceleration. */ // TODO: this is missing for solid body.
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    /** Evaluation of fluid density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_fluid_density(water_block_inner, water_block_contact);
    /** Compute time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_ref);
    /** Compute time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec2d, DampingPairwiseInner>>>
        fluid_damping(0.2, water_block_inner, water_block_contact, "Velocity", mu_f);
    //----------------------------------------------------------------------
    //	Define solid methods which are used in this case.
    //----------------------------------------------------------------------
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> gate_corrected_configuration(gate_inner);
    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> gate_computing_time_step_size(gate);
    /** Stress relaxation stepping for the elastic gate. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> gate_stress_relaxation_first_half(gate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> gate_stress_relaxation_second_half(gate_inner);
    /**Constrain a solid body part.  */
    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> gate_constraint(gate_constraint_part);
    /** Update the norm of elastic gate. */
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> gate_update_normal(gate);
    /** Curvature calculation for elastic gate. */
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> gate_curvature(shell_curvature_inner);
    /**Damping.  */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        gate_position_damping(0.2, gate_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        gate_rotation_damping(0.2, gate_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define fsi methods which are used in this case.
    //----------------------------------------------------------------------
    /** Compute the average velocity of gate. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
    /** Compute the force exerted on elastic gate due to fluid pressure. */
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_pressure_force_on_gate(gate_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    gate.addBodyStateForRecording<Real>("AverageTotalMeanCurvature");
    gate.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    gate.addBodyStateForRecording<Vecd>("PriorForce");
    /** Output body states for visualization. */
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(io_environment, sph_system.real_bodies_);
    /** Output the observed displacement of gate center. */
    ObservedQuantityRecording<Vecd>
        write_beam_tip_displacement("Displacement", io_environment, gate_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing linear reproducing configuration for the insert body. */
    gate_corrected_configuration.exec();
    /** calculate initial curvature after corrected configuration*/
    gate_curvature.exec();
    /** update fluid-shell contact*/
    water_block_contact.updateConfiguration();

    // Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness(water_block_inner,
                                                      water_block_contact);
    check_kernel_completeness.exec();
    water_block.addBodyStateForRecording<Real>("TotalKernel");
    water_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.5; /**< End time. */
    Real output_interval = end_time / 50.0;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force and gravity. */
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            /** Update normal direction on elastic body. */
            gate_update_normal.exec();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                fluid_damping.exec(dt);
                /** Fluid relaxation and force computation. */
                pressure_relaxation.exec(dt);
                fluid_pressure_force_on_gate.exec();
                density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    gate_stress_relaxation_first_half.exec(dt_s);
                    gate_constraint.exec();
                    gate_position_damping.exec(dt_s);
                    gate_rotation_damping.exec(dt_s);
                    gate_constraint.exec();
                    gate_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    dt_s = gate_computing_time_step_size.exec();
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update curvature. */
            gate_curvature.exec();
            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList(); // water particle motion is small
            gate.updateCellLinkedList();
            /** one need update configuration after periodic condition. */
            water_block_complex.updateConfiguration();
            gate_contact.updateConfiguration();

            /** Output the observed data. */
            write_beam_tip_displacement.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        write_real_body_states_to_vtp.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    // Compare with analytical solution
    const Real p = rho0_f * gravity_g * Dam_H;
    const Real I = 1 / 12.0 * 1.0 * std::pow(Gate_thickness, 3);
    const Real max_disp_analytical = p * std::pow(Dam_L, 4) / 384.0 / Youngs_modulus / I;
    const Real max_disp = std::abs((*write_beam_tip_displacement.interpolated_quantities_)[0][1]);
    const Real error = std::abs((max_disp_analytical - max_disp) / max_disp_analytical) * 100.0;

    std::cout << "Analytical displacement: " << max_disp_analytical
              << "\t Displacement: " << max_disp
              << "\t Error: " << error << "%"
              << std::endl;

    // gtest
    EXPECT_NEAR(max_disp_analytical, max_disp, max_disp_analytical * 20e-2); // it's below 5% but 20% for CI

    return 0;
}
