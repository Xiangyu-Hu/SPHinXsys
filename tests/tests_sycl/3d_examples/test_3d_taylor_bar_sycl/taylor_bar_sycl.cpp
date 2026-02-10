/**
 * @file taylor_bar_sycl.cpp
 * @brief This test offloads the levelset computation to the GPU device.
 *        All other components remain unchanged from the original Taylor bar test.
 *        The regression test is also adapted to validate the results.
 * @author Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 */
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 0.0032; /**< X-direction domain. */
Real PW = 0.0324; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 5.0;
/** YOU can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real SL = particle_spacing_ref * 4.0; /**< Length of the holder is one layer particle. */
Real column_radius = PL;
Vec3d translation_column(0.0, 0.0, 0.6 * PW);
Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
Vecd translation_holder(0.0, 0.0, -0.5 * SL);
Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; /**< Reference density. */
Real poisson = 0.35;  /**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;
class InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][2] = -227.0;
    }
};

/**
 * @class DynamicContactForceWithWall
 * @brief Computing the contact force with a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class DynamicContactForceWithWall : public LocalDynamics,
                                    public DataDelegateContact
{
  public:
    explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0)
        : LocalDynamics(solid_body_contact_relation.getSPHBody()),
          DataDelegateContact(solid_body_contact_relation),
          solid_(DynamicCast<Solid>(this, sph_body_->getBaseMaterial())),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          penalty_strength_(penalty_strength)
    {
        impedance_ = sqrt(solid_.ReferenceDensity() * solid_.ContactStiffness());
        reference_pressure_ = solid_.ContactStiffness();
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_n_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        }
    };
    virtual ~DynamicContactForceWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->getSPHAdaptation().ReferenceSpacing();
            Real particle_spacing_ratio2 =
                1.0 / (getSPHAdaptation().ReferenceSpacing() * particle_spacing_j1);
            particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

            Vecd *n_k = contact_n_[k];
            Vecd *vel_n_k = contact_vel_[k];
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Vecd n_k_j = n_k[index_j];

                Real impedance_p = 0.5 * impedance_ * (vel_[index_i] - vel_n_k[index_j]).dot(-n_k_j);
                Real overlap = contact_neighborhood.r_ij_[n] * n_k_j.dot(e_ij);
                Real delta = 2.0 * overlap * particle_spacing_j1;
                Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
                Real penalty_p = penalty_strength_ * beta * fabs(overlap) * reference_pressure_;

                // force due to pressure
                force -= 2.0 * (impedance_p + penalty_p) * e_ij.dot(n_k_j) *
                         n_k_j * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            }
        }

        force_prior_[index_i] += force * Vol_[index_i];
    };

  protected:
    Solid &solid_;
    Real *Vol_;
    Vecd *vel_, *force_prior_; // note that prior force directly used here
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_n_;
    Real penalty_strength_;
    Real impedance_, reference_pressure_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Setup geometry first.
    //----------------------------------------------------------------------
    auto &column_shape = sph_system.addShape<TriangleMeshShapeCylinder>(
        Vec3d(0, 0, 1.0), column_radius, 0.5 * PW, resolution, translation_column, "Column");
    auto &wall_shape = sph_system.addShape<TriangleMeshShapeBrick>(
        halfsize_holder, resolution, translation_holder, "Wall");
    //----------------------------------------------------------------------
    //	Run particle relaxation first if needed.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        // setup a sub-system for particle relaxation and delete it after particle relaxation.
        SPHSystem relaxation_system(system_domain_bounds, particle_spacing_ref);
        auto &column = relaxation_system.addBody<RealBody>(column_shape);
        auto &wall = relaxation_system.addBody<SolidBody>(wall_shape);

        LevelSetShape *level_set_shape = column.defineBodyLevelSetShape(par_ck, 2.0)->writeLevelSet();
        column.generateParticles<BaseParticles, Lattice>();
        wall.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_body_surface(column);
        Inner<> column_inner(column);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(relaxation_system);
        auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
        auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

        auto &input_body_cell_linked_list = main_methods.addCellLinkedListDynamics(column);
        auto &input_body_update_inner_relation = main_methods.addRelationDynamics(column_inner);
        auto &random_input_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(column);
        auto &relaxation_residual =
            main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(column_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(column, *level_set_shape);
        auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(column);
        auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(column);
        auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
        //----------------------------------------------------------------------
        //	Run on CPU after relaxation finished and output results.
        //----------------------------------------------------------------------
        auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(relaxation_system);
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&column, &wall});
        write_particle_reload_files.addToReload<Vecd>(wall, "NormalDirection");
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        random_input_body_particles.exec();

        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        body_state_recorder.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            input_body_cell_linked_list.exec();
            input_body_update_inner_relation.exec();

            relaxation_residual.exec();
            Real relaxation_step = relaxation_scaling.exec();
            update_particle_position.exec(relaxation_step);
            level_set_bounding.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;
        wall_boundary_normal_direction.exec();
        write_particle_reload_files.writeToFile();

        if (!sph_system.ReloadParticles())
        {
            return 0;
        }
        else
        {
            std::cout << "To reload particles and start the main simulation." << std::endl;
        }
    }
    //----------------------------------------------------------------------
    //	Simulation setup start here.
    //----------------------------------------------------------------------
    auto &column = sph_system.addBody<RealBody>(column_shape);
    column.defineMaterial<HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
    column.generateParticles<BaseParticles, Reload>(column.getName());

    auto &wall = sph_system.addBody<RealBody>(wall_shape);
    wall.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall.generateParticles<BaseParticles, Reload>(wall.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");

    /** Define Observer. */
    ObserverBody my_observer(sph_system, "MyObserver");
    StdVec<Vecd> observation_location = {Vecd(0.0, 0.0, PW)};
    my_observer.generateParticles<ObserverParticles>(observation_location);
    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    SurfaceContactRelation column_wall_contact(column, {&wall});
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(sph_system);
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<::InitialCondition> initial_condition(column);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(column_inner);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf> stress_relaxation_first_half(column_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(column_inner);
    InteractionDynamics<DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(column, 0.2);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", my_observer_contact);

    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    corrected_configuration.exec();
    initial_condition.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 1.0e-4;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real output_period = 1.0e-6; // anyway 50 write_states files in total
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % screen_output_interval == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";

                if (ite != 0 && ite % observation_sample_interval == 0)
                {
                    write_displacement.writeToFile(ite);
                }
            }
            column_wall_contact_force.exec(dt);
            stress_relaxation_first_half.exec(dt);
            stress_relaxation_second_half.exec(dt);

            column.updateCellLinkedList();
            column_wall_contact.updateConfiguration();

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_displacement.generateDataBase(0.1);
    }
    else
    {
        write_displacement.testResult();
    }

    return 0;
}
