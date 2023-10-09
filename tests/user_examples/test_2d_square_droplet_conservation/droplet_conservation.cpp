/**
 * @file 	droplet.cpp
 * @brief 	A square droplet deforms to circle due to surface tension.
 * @details A momentum-conservative formulation for surface twnsion is used here 
 *          to reach a long-rerm stable simulation.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.

using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;                         /**< Tank length. */
Real DH = 2.0;                         /**< Tank height. */
Real LL = 1.0;                         /**< Liquid column length. */
Real LH = 1.0;                         /**< Liquid column height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of water. */
Real rho0_a = 0.001;     /**< Reference density of air. */
Real U_ref = 2.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
Real mu_f = 0.2;         /**< Water viscosity. */
Real mu_a = 0.002;       /**< Air viscosity. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
Vecd air_halfsize = inner_wall_halfsize; // local center at origin
Vecd air_translation = inner_wall_translation;   // translation to global coordinates

Vecd droplet_center(DL / 2, DH / 2);
Real drolet_radius = LL / 2;
Vecd droplet_halfsize = Vec2d(drolet_radius, drolet_radius); // local center at origin
Vecd droplet_translation = droplet_center;   // translation to global coordinates
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(droplet_translation), droplet_halfsize);
    }
};
//----------------------------------------------------------------------
// Air body shape definition.
//----------------------------------------------------------------------/**
class AirBlock : public ComplexShape
{
public:
    explicit AirBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(air_translation), air_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(droplet_translation), droplet_halfsize);
    }
};
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};

typedef DataDelegateContact<BaseParticles, BaseParticles> BaseDataContact;
class SurfaceTensionStress : public LocalDynamics, public BaseDataContact
{
public:
    explicit SurfaceTensionStress(BaseContactRelation& conact_relation, StdVec<Real> contact_surface_tension)
        : LocalDynamics(conact_relation.getSPHBody()), BaseDataContact(conact_relation)
    {
        particles_->registerVariable(color_gradient_, "ColorGradient");
        particles_->registerSortableVariable<Vecd>("ColorGradient");
        particles_->addVariableToWrite<Vecd>("ColorGradient");
        particles_->registerVariable(surface_tension_stress_, "SurfaceTensionStress");
        particles_->registerSortableVariable<Matd>("SurfaceTensionStress");
        particles_->addVariableToWrite<Matd>("SurfaceTensionStress");
        Real rho0 = getSPHBody().base_material_->ReferenceDensity();
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_surface_tension_.push_back(contact_surface_tension[k]);
            Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
            contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
        }
    };
    virtual ~SurfaceTensionStress() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        color_gradient_[index_i] = ZeroData<Vecd>::value;
        surface_tension_stress_[index_i] = ZeroData<Matd>::value;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd weighted_color_gradeint = ZeroData<Vecd>::value;
            Real contact_fraction_k = contact_fraction_[k];
            Real surface_tension_k = contact_surface_tension_[k];
            const Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                weighted_color_gradeint -= contact_fraction_k *
                    contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
            color_gradient_[index_i] = weighted_color_gradeint;
            Real norm = weighted_color_gradeint.norm();
            surface_tension_stress_[index_i] += surface_tension_k / (norm + Eps) *
                (norm * norm * Matd::Identity() -
                    weighted_color_gradeint * weighted_color_gradeint.transpose());
        }
    };

protected:
    StdLargeVec<Vecd> color_gradient_;
    StdLargeVec<Matd> surface_tension_stress_;
    StdVec<Real> contact_surface_tension_, contact_fraction_;
};

using fluid_dynamics::FluidDataInner;
class SurfaceStressAccelerationInner : public LocalDynamics, public FluidDataInner
{
public:
    SurfaceStressAccelerationInner(BaseInnerRelation& inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
        rho_(particles_->rho_), acc_prior_(particles_->acc_prior_),
        color_gradient_(*particles_->getVariableByName<Vecd>("ColorGradient")),
        surface_tension_stress_(*particles_->getVariableByName<Matd>("SurfaceTensionStress")) {};
    virtual ~SurfaceStressAccelerationInner() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd summation = ZeroData<Vecd>::value;
        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            summation += inner_neighborhood.dW_ijV_j_[n] *
                (surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
                inner_neighborhood.e_ij_[n];
        }
        acc_prior_[index_i] += summation / rho_[index_i];
    };

protected:
    StdLargeVec<Real>& rho_;
    StdLargeVec<Vecd>& acc_prior_;
    StdLargeVec<Vecd>& color_gradient_;
    StdLargeVec<Matd>& surface_tension_stress_;
};

class SurfaceStressAccelerationContact : public LocalDynamics, public BaseDataContact
{
public:
    explicit SurfaceStressAccelerationContact(BaseContactRelation& conact_relation)
        : LocalDynamics(conact_relation.getSPHBody()), BaseDataContact(conact_relation),
        rho_(particles_->rho_), acc_prior_(particles_->acc_prior_),
        color_gradient_(*particles_->getVariableByName<Vecd>("ColorGradient")),
        surface_tension_stress_(*particles_->getVariableByName<Matd>("SurfaceTensionStress"))
    {
        Real rho0 = getSPHBody().base_material_->ReferenceDensity();
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
            contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
            contact_color_gradient_.push_back(
                contact_particles_[k]->getVariableByName<Vecd>("ColorGradient"));
            contact_surface_tension_stress_.push_back(
                contact_particles_[k]->getVariableByName<Matd>("SurfaceTensionStress"));
        }
    };
    virtual ~SurfaceStressAccelerationContact() {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd summation = ZeroData<Vecd>::value;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real contact_fraction_k = contact_fraction_[k];
            StdLargeVec<Vecd>& contact_color_gradient_k = *(contact_color_gradient_[k]);
            StdLargeVec<Matd>& contact_surface_tension_stress_k = *(contact_surface_tension_stress_[k]);
            const Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real r_ij = contact_neighborhood.r_ij_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Real mismatch = 1.0 - 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]).dot(e_ij) * r_ij;
                summation += contact_neighborhood.dW_ijV_j_[n] *
                    (-0.1 * mismatch * Matd::Identity() +
                        (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] +
                        contact_surface_tension_stress_k[index_j] * contact_fraction_k) *
                    contact_neighborhood.e_ij_[n];
            }
        }
        acc_prior_[index_i] += summation / rho_[index_i];
    };

protected:
    StdLargeVec<Real>& rho_;
    StdLargeVec<Vecd>& acc_prior_;
    StdLargeVec<Vecd>& color_gradient_;
    StdLargeVec<Matd>& surface_tension_stress_;
    StdVec<StdLargeVec<Vecd>*> contact_color_gradient_;
    StdVec<StdLargeVec<Matd>*> contact_surface_tension_stress_;
    StdVec<Real> contact_surface_tension_, contact_fraction_;
};

template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class InteractionWithInitialization : public InteractionDynamics<LocalDynamicsType, ExecutionPolicy>
{
public:
    template <typename... Args>
    InteractionWithInitialization(Args &&...args)
        : InteractionDynamics<LocalDynamicsType, ExecutionPolicy>(false, std::forward<Args>(args)...)
    {
        static_assert(!has_update<LocalDynamicsType>::value,
            "LocalDynamicsType does not fulfill InteractionWithInitialization requirements");
    }
    virtual ~InteractionWithInitialization() {};

    virtual void exec(Real dt = 0.0) override
    {
        particle_for(ExecutionPolicy(),
            this->identifier_.LoopRange(),
            [&](size_t i)
        { this->initialization(i, dt); });
        InteractionDynamics<LocalDynamicsType, ExecutionPolicy>::exec(dt);
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<Real>("Pressure");

    FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
    air_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
    air_block.generateParticles<ParticleGeneratorLattice>();
    air_block.addBodyStateForRecording<Real>("Density");
    air_block.addBodyStateForRecording<Real>("Pressure");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_block, { &air_block });
    ContactRelation water_wall_contact(water_block, { &wall_boundary });
    ComplexRelation air_water_complex(air_block, { &water_block });
    ContactRelation air_wall_contact(air_block, { &wall_boundary });
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block);
    SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
        update_air_density_by_summation(air_wall_contact, air_water_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
        update_water_density_by_summation(water_wall_contact, water_air_complex);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_wall_contact, air_water_complex, 0.02);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
        water_transport_correction(water_air_complex, 0.02);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_ref);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);
    /** Pressure relaxation for water by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
        water_pressure_relaxation(water_wall_contact, water_air_complex);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
        water_density_relaxation(water_wall_contact, water_air_complex);
    /** Extend Pressure relaxation is used for air. */
    Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
        air_pressure_relaxation(air_wall_contact, air_water_complex, 1.0);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
        air_density_relaxation(air_wall_contact, air_water_complex);
    /** Viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall>
        air_viscous_acceleration(air_wall_contact, air_water_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall>
        water_viscous_acceleration(water_wall_contact, water_air_complex);
    /** Surface tension. */
    InteractionDynamics<SurfaceTensionStress>
        water_surface_tension_stress(water_air_complex.getContactRelation(), StdVec<Real>{Real(1.0)});
    InteractionDynamics<SurfaceTensionStress>
        air_surface_tension_stress(air_water_complex.getContactRelation(), StdVec<Real>{Real(1.0e-3)});
    InteractionDynamics<ComplexInteraction<SurfaceStressAccelerationInner, SurfaceStressAccelerationContact>>
        water_surface_tension_acceleration(water_air_complex.getInnerRelation(), water_air_complex.getContactRelation());
    water_block.addBodyStateForRecording<Matd>("SurfaceTensionStress");
    InteractionDynamics<ComplexInteraction<SurfaceStressAccelerationInner, SurfaceStressAccelerationContact>>
        air_surface_tension_acceleration(air_water_complex.getInnerRelation(), air_water_complex.getContactRelation());
    air_block.addBodyStateForRecording<Vecd>("PriorAcceleration");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_water_mechanical_energy(io_environment, water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    Real end_time = 10.0;
    Real output_interval = end_time / 50; /**< Time stamps for output of body states. */
    Real dt = 0.0;               /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force and gravity. */
            time_instance = TickCount::now();
            initialize_a_water_step.exec();
            initialize_a_air_step.exec();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_f, Dt_a);

            update_air_density_by_summation.exec();
            update_water_density_by_summation.exec();
            air_transport_correction.exec();
            water_transport_correction.exec();

            air_viscous_acceleration.exec();
            water_viscous_acceleration.exec();

            water_surface_tension_stress.exec();
            air_surface_tension_stress.exec();
            water_surface_tension_acceleration.exec();
            air_surface_tension_acceleration.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                    << GlobalStaticVariables::physical_time_
                    << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            water_block.updateCellLinkedList();
            water_air_complex.updateConfiguration();
            water_wall_contact.updateConfiguration();

            air_block.updateCellLinkedList();
            air_water_complex.updateConfiguration();
            air_wall_contact.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
        }
        write_water_mechanical_energy.writeToFile(number_of_iterations);
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
        << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
        << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
        << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
        << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_water_mechanical_energy.testResult();
    }

    return 0;
}