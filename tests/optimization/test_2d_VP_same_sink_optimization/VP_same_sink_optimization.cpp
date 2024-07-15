/**
 * @file 	VP_same_sink_optimized.cpp
 * @brief 	This is the optimized test for the same sink (2/10) temperature.
 * @author 	Bo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 50.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 0.0;
Real high_temperature = 300.0;
Real low_temperature = 300.0;
Real heat_source = 1000.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
    std::vector<Vecd> thermalDomainShape;
    thermalDomainShape.push_back(Vecd(0.0, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, H));
    thermalDomainShape.push_back(Vecd(L, H));
    thermalDomainShape.push_back(Vecd(L, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, 0.0));

    return thermalDomainShape;
};

std::vector<Vecd> createBoundaryDomain()
{
    std::vector<Vecd> boundaryDomain;
    boundaryDomain.push_back(Vecd(-BW, -BW));
    boundaryDomain.push_back(Vecd(-BW, H + BW));
    boundaryDomain.push_back(Vecd(L + BW, H + BW));
    boundaryDomain.push_back(Vecd(L + BW, -BW));
    boundaryDomain.push_back(Vecd(-BW, -BW));

    return boundaryDomain;
};
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createBoundaryDomain(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DiffusionBodyInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit DiffusionBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          pos_(*particles_->getVariableDataByName<Vecd>("Position")),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
          heat_source_(*(particles_->registerSharedVariable<Real>("HeatSource"))){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = 650;
        heat_source_[index_i] = heat_source;
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &phi_, &heat_source_;
};

class ThermalConductivityRandomInitialization : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermalConductivityRandomInitialization(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          thermal_conductivity(*(particles_->getVariableDataByName<Real>("ThermalConductivity"))){};
    void update(size_t index_i, Real dt)
    {
        thermal_conductivity[index_i] = 0.5 + rand_uniform(0.0, 1.0);
    };

  protected:
    StdLargeVec<Real> &thermal_conductivity;
};

class WallBoundaryInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit WallBoundaryInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          pos_(*particles_->getVariableDataByName<Vecd>("Position")),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;
        if (pos_[index_i][1] < 0.0 && pos_[index_i][00] > 0.4 * L && pos_[index_i][0] < 0.6 * L)
        {
            phi_[index_i] = low_temperature;
        }
        if (pos_[index_i][1] > 1.0 && pos_[index_i][0] > 0.4 * L && pos_[index_i][0] < 0.6 * L)
        {
            phi_[index_i] = high_temperature;
        }
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &phi_;
};
//----------------------------------------------------------------------
//  Impose constraints on the objective function
//----------------------------------------------------------------------
class ImposeObjectiveFunction : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ImposeObjectiveFunction(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
          species_modified_(*particles_->getVariableDataByName<Real>("SpeciesModified")),
          species_recovery_(*particles_->getVariableDataByName<Real>("SpeciesRecovery")){};

    void update(size_t index_i, Real learning_rate)
    {
        species_recovery_[index_i] = phi_[index_i];
        species_modified_[index_i] = phi_[index_i] - learning_rate;
    };

  protected:
    StdLargeVec<Real> &phi_, &species_modified_, &species_recovery_;
};

class StoreGlobalPDEResidual : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit StoreGlobalPDEResidual(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          residual_T_local_(*particles_->getVariableDataByName<Real>("ResidualTLocal")),
          residual_T_global_(*particles_->getVariableDataByName<Real>("ResidualTGlobal")){};

    void update(size_t index_i, Real learning_rate)
    {
        residual_T_local_[index_i] = residual_T_global_[index_i];
    };

  protected:
    StdLargeVec<Real> &residual_T_local_, &residual_T_global_;
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
TEST(test_optimization, test_problem1_optimized)
{
    /* Here is the initialization of the random generator. */
    srand((double)clock());
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    LocalIsotropicDiffusion *local_isotropic_diffusion =
        diffusion_body.defineMaterial<LocalIsotropicDiffusion>("Phi", "Phi", diffusion_coeff);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation diffusion_body_inner(diffusion_body);
    ContactRelation diffusion_body_contact(diffusion_body, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation diffusion_body_complex(diffusion_body_inner, diffusion_body_contact);
    //----------------------------------------------------------------------
    // Obtain the time step size.
    //----------------------------------------------------------------------
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, *local_isotropic_diffusion);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RestartIO restart_io(sph_system);
    //----------------------------------------------------------------------
    //	Setup parameter for optimization control
    //----------------------------------------------------------------------
    int ite = 0;                  /* define the loop of all operations for optimization. */
    int ite_T = 0;                /* define loop index for temperature splitting iteration. */
    int ite_k = 0;                /* define loop index for parameter splitting iteration. */
    int ite_rg = 0;               /* define loop index for parameter regularization. */
    int ite_T_total = 1;          /* define the total iteration for temperature splitting. */
    int ite_k_total = 1;          /* define the total iteration for parameter splitting. */
    int ite_loop = 0;             /* define loop index for optimization cycle. */
    int ite_T_comparison_opt = 0; /* define the real step for splitting temperature by solving PDE. */
    int ite_output = 50;          /* define the interval for state output. */
    int ite_restart = 50;         /* define the interval for restart output. */
    int dt_ratio_k = 1;           /* ratio for adjusting the time step for parameter evolution. */
    int dt_ratio_rg = 1;          /* ratio for adjusting the time step for regularization. */

    Real dt = get_time_step_size.exec();
    Real averaged_residual_T_last_global(10.0);
    Real averaged_variation_last_global(10.0);
    Real averaged_residual_T_current_global(0.0);
    Real averaged_variation_current_global(0.0);
    Real maximum_variation_current_global(10.0);
    Real opt_averaged_temperature = 0.0;
    Real nonopt_averaged_temperature = MaxReal;
    Real averaged_k_parameter = 0.0;
    Real initial_eta_regularization = 0.4;
    Real current_eta_regularization = initial_eta_regularization;
    Real relative_temperature_difference = 2.0;
    Real last_averaged_temperature = 0.0;
    Real current_averaged_temperature = 0.0;
    Real relative_average_variation_difference = 1.0;

    /* Gradient descent parameter for objective function.*/
    Real initial_learning_rate = 0.2;
    Real learning_rate_alpha = initial_learning_rate;
    //----------------------------  ------------------------------------------
    //	Define the main numerical methods used for optimization
    //  Note that there may be data dependence on the constructors of tested methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> diffusion_body_normal_direction(diffusion_body);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    InteractionSplit<TemperatureSplittingByPDEWithBoundary<Real>>
        temperature_splitting_pde_complex(diffusion_body_inner, diffusion_body_contact, "Phi");
    InteractionSplit<UpdateTemperaturePDEResidual<TemperatureSplittingByPDEWithBoundary<Real>>>
        update_temperature_pde_residual(diffusion_body_inner, diffusion_body_contact, "Phi");

    SimpleDynamics<StoreGlobalPDEResidual> store_global_PDE_residual(diffusion_body);
    SimpleDynamics<ImposeObjectiveFunction> impose_objective_function(diffusion_body);
    InteractionSplit<ParameterSplittingByPDEWithBoundary<Real>>
        parameter_splitting_pde_complex(diffusion_body_inner, diffusion_body_contact, "ThermalConductivity");
    InteractionSplit<RegularizationByDiffusionAnalogy<Real>>
        thermal_diffusivity_regularization(diffusion_body_inner, "ThermalConductivity",
                                           initial_eta_regularization, maximum_variation_current_global);
    InteractionSplit<UpdateRegularizationVariation<Real>>
        update_regularization_global_variation(diffusion_body_inner, "ThermalConductivity");
    ReduceDynamics<Average<ComputeTotalErrorOrPositiveParameter<SPHBody>>>
        total_averaged_thermal_diffusivity(diffusion_body, "ThermalConductivity");
    SimpleDynamics<ThermalConductivityConstraint<SPHBody>>
        thermal_diffusivity_constrain(diffusion_body, "ThermalConductivity");
    ReduceDynamics<Average<ComputeTotalErrorOrPositiveParameter<SPHBody>>>
        calculate_temperature_global_residual(diffusion_body, "ResidualTGlobal");
    ReduceDynamics<Average<ComputeTotalErrorOrPositiveParameter<SPHBody>>>
        calculate_regularization_global_variation(diffusion_body, "VariationGlobal");
    ReduceDynamics<ComputeMaximumError<SPHBody>>
        calculate_maximum_variation(diffusion_body, "VariationGlobal");
    ReduceDynamics<Average<QuantitySummation<Real, SPHBody>>>
        calculate_averaged_opt_temperature(diffusion_body, "Phi");
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<DiffusionBodyInitialCondition> setup_diffusion_initial_condition(diffusion_body);
    SimpleDynamics<WallBoundaryInitialCondition> setup_diffusion_boundary_condition(wall_boundary);
    SimpleDynamics<ThermalConductivityRandomInitialization> thermal_diffusivity_random_initialization(diffusion_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    diffusion_body_normal_direction.exec();
    wall_boundary_normal_direction.exec();
    setup_diffusion_initial_condition.exec();
    setup_diffusion_boundary_condition.exec();
    thermal_diffusivity_random_initialization.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        diffusion_body.updateCellLinkedList();
        diffusion_body_complex.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    // record the temperature with modified design variable.
    std::string filefullpath_opt_temperature =
        sph_system.getIOEnvironment().output_folder_ + "/" + "opt_temperature.dat";
    std::ofstream out_file_opt_temperature(filefullpath_opt_temperature.c_str(), std::ios::app);

    // record the temperature without modified design variable.
    std::string filefullpath_nonopt_temperature =
        sph_system.getIOEnvironment().output_folder_ + "/" + "nonopt_temperature.dat";
    std::ofstream out_file_nonopt_temperature(filefullpath_nonopt_temperature.c_str(), std::ios::app);
    //----------------------------------------------------------------------
    //	Initial States update.
    //----------------------------------------------------------------------
    write_states.writeToFile(ite); // output the initial states.

    update_regularization_global_variation.exec(dt_ratio_rg * dt);
    averaged_variation_current_global = calculate_regularization_global_variation.exec(dt);
    maximum_variation_current_global = calculate_maximum_variation.exec(dt);

    update_temperature_pde_residual.exec(dt);
    averaged_residual_T_current_global = calculate_temperature_global_residual.exec(dt);
    averaged_residual_T_last_global = averaged_residual_T_current_global;

    current_averaged_temperature = calculate_averaged_opt_temperature.exec();
    out_file_nonopt_temperature << std::fixed << std::setprecision(12) << ite << "   " << current_averaged_temperature << "\n";
    out_file_opt_temperature << std::fixed << std::setprecision(12) << ite_T_comparison_opt << "   " << current_averaged_temperature << "\n";

    /** the converged criterion contains three parts respect to target function, PDE constrain, and maximum step. */
    while ((relative_temperature_difference > 0.00001 || averaged_residual_T_current_global > 0.000005 ||
            relative_average_variation_difference > 0.0001) &&
           ite_loop < 10000)
    {
        std::cout << "This is the beginning of the " << ite_loop << " iteration loop." << std::endl;

        //----------------------------------------------------------------------
        //	Impose Objective Function.
        //----------------------------------------------------------------------
        ite++;

        /* Store the global PDE residual to provide the reference for design variable splitting based on PDE. */
        store_global_PDE_residual.exec();

        /* Impose objective function and PDE residual may increase. */
        impose_objective_function.exec(learning_rate_alpha);

        // if (ite_loop % ite_output == 0) { write_states.writeToFile(ite); }
        std::cout << "N=" << ite << " and the objective function has been imposed. "
                  << "\n";

        //----------------------------------------------------------------------
        //	Parameter (design variable) splitting.
        //----------------------------------------------------------------------
        /* Parameter splitting should recovery the increased residual by imposing objective function. */
        while (ite_k < ite_k_total)
        {
            //----------------------------------------------------------------------
            //	Parameter splitting by PDE.
            //----------------------------------------------------------------------
            ite++;
            ite_k++;
            parameter_splitting_pde_complex.exec(dt_ratio_k * dt);
            //----------------------------------------------------------------------
            //	Constraint of the summation of parameter.
            //----------------------------------------------------------------------
            if (ite_k % 1 == 0 || ite_k == ite_k_total)
            {
                ite++;
                averaged_k_parameter = total_averaged_thermal_diffusivity.exec(dt);
                thermal_diffusivity_constrain.UpdateAverageParameter(averaged_k_parameter);
                thermal_diffusivity_constrain.exec(dt);
            }

            //----------------------------------------------------------------------
            //	Regularization.
            //----------------------------------------------------------------------
            if (ite_k % 1 == 0 || ite_k == ite_k_total)
            {
                ite++;
                ite_rg++;
                thermal_diffusivity_regularization.UpdateCurrentEta(current_eta_regularization);
                thermal_diffusivity_regularization.UpdateMaximumVariation(maximum_variation_current_global);
                thermal_diffusivity_regularization.UpdateAverageVariation(averaged_variation_current_global);
                thermal_diffusivity_regularization.exec(dt_ratio_rg * dt);

                update_temperature_pde_residual.exec(dt);
                averaged_residual_T_current_global = calculate_temperature_global_residual.exec(dt);

                update_regularization_global_variation.exec(dt);
                averaged_variation_current_global = calculate_regularization_global_variation.exec(dt);
                maximum_variation_current_global = calculate_maximum_variation.exec(dt);
            }
        }
        ite_k = 0;
        ite_rg = 0;
        if (ite_loop % ite_output == 0)
        {
            write_states.writeToFile(ite);
        }
        std::cout << "N=" << ite << " and the k splitting is finished."
                  << "\n";

        //----------------------------------------------------------------------
        //	Temperature splitting.
        //----------------------------------------------------------------------
        std::cout << "averaged_residual_T_last_global is " << averaged_residual_T_last_global << std::endl;
        while (((averaged_residual_T_current_global > 0.9 * averaged_residual_T_last_global) &&
                averaged_residual_T_current_global > 0.000005) ||
               ite_T < ite_T_total)
        {
            ite++;
            ite_T++;
            ite_T_comparison_opt++;
            temperature_splitting_pde_complex.exec(dt);

            update_temperature_pde_residual.exec(dt);
            averaged_residual_T_current_global = calculate_temperature_global_residual.exec(dt);
        }

        opt_averaged_temperature = calculate_averaged_opt_temperature.exec(dt);
        out_file_opt_temperature << std::fixed << std::setprecision(12) << ite_T_comparison_opt << "   " << opt_averaged_temperature << "\n";
        out_file_nonopt_temperature << std::fixed << std::setprecision(12) << ite << "   " << opt_averaged_temperature << "\n";

        /* Decay learning rate by learning process. */
        if (nonopt_averaged_temperature > opt_averaged_temperature)
        {
            learning_rate_alpha = 1.05 * learning_rate_alpha;
            current_eta_regularization = 1.05 * current_eta_regularization;
            std::cout << "The learning rate is fixed!" << std::endl;
        }
        else
        {
            learning_rate_alpha = 0.8 * learning_rate_alpha;
            current_eta_regularization = 0.8 * current_eta_regularization;
            std::cout << "The learning rate is decreased by optimization process!" << std::endl;
        }

        nonopt_averaged_temperature = opt_averaged_temperature;
        averaged_residual_T_last_global = averaged_residual_T_current_global;

        std::cout << "averaged_residual_T_current_global is " << averaged_residual_T_current_global << std::endl;
        ite_T = 0;
        write_states.writeToFile(ite);
        std::cout << "N=" << ite << " and the temperature splitting is finished."
                  << "\n";

        //----------------------------------------------------------------------
        //	Decision Making.
        //----------------------------------------------------------------------
        last_averaged_temperature = current_averaged_temperature;
        current_averaged_temperature = calculate_averaged_opt_temperature.exec();

        ite_loop++;
        std::cout << "This is the " << ite_loop << " iteration loop and the averaged temperature is " << opt_averaged_temperature
                  << " and the learning rate is " << learning_rate_alpha
                  << " and the regularization is " << current_eta_regularization << std::endl;
        relative_temperature_difference = abs(current_averaged_temperature - last_averaged_temperature) / last_averaged_temperature;
        relative_average_variation_difference = abs(averaged_variation_current_global - averaged_variation_last_global) / abs(averaged_variation_last_global);
        averaged_variation_last_global = averaged_variation_current_global;
        if (ite_loop % ite_restart == 0)
        {
            restart_io.writeToFile(ite_loop);
        }
    }
    out_file_opt_temperature.close();
    out_file_nonopt_temperature.close();

    TickCount t2 = TickCount::now();
    TickCount::interval_t tt;
    tt = t2 - t1;
    std::cout << "Total time for optimization: " << tt.seconds() << " seconds." << std::endl;

    EXPECT_GT(502.0, calculate_averaged_opt_temperature.exec());
};

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}