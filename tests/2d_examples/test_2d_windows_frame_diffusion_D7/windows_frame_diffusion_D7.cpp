/*
 * @file 	windows_frame_diffusion_D7.cpp
 * @brief 	Heat diffusion of windows frame based on international standard ISO 10077-2:2012,
 *          Application 7: PVC frame section and insulation panel.
 * @details https://www.comsol.de/model/thermal-performances-of-windows-16077
 * @author 	Haotian Ji, Dong Wu, Chi Zhang and Xiangyu Hu
 */
#include "windows_frame_diffusion_D7.h"
#include "sphinxsys.h"

using namespace SPH; // Namespace cite here
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<MultiPolygonShape>(createOverallStructureBody(), "DiffusionBody"));
    LocalIsotropicDiffusion *frame_diffusion = diffusion_body.defineMaterial<LocalIsotropicDiffusion>("Phi", "Phi", pvc_cond);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody boundary_Robin_in(sph_system, makeShared<MultiPolygonShape>(createInternalAirBody(), "InternalConvectionBoundary"));
    boundary_Robin_in.generateParticles<Lattice>();

    SolidBody boundary_Robin_ex(sph_system, makeShared<MultiPolygonShape>(createExternalAirBody(), "ExternalConvectionBoundary"));
    boundary_Robin_ex.generateParticles<Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of temperature observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles(ParticleGeneratorTemperatureObserver(temperature_observer));
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation inner_relation(diffusion_body);
    ContactRelation contact_Robin_in(diffusion_body, {&boundary_Robin_in});
    ContactRelation contact_Robin_ex(diffusion_body, {&boundary_Robin_ex});
    ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    DiffusionBodyRelaxation temperature_relaxation(inner_relation, contact_Robin_in, contact_Robin_ex);

    RobinFluxCalculation heat_flux_calculation_in(contact_Robin_in);
    RobinFluxCalculation heat_flux_calculation_ex(contact_Robin_ex);

    GetDiffusionTimeStepSize<WallParticles> get_time_step_size(boundary_Robin_ex);

    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_in(boundary_Robin_in);
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_ex(boundary_Robin_ex);

    SimpleDynamics<InternalHeatTransferInitialCondition> setup_heat_transfer_initial_condition_in(boundary_Robin_in);
    SimpleDynamics<ExternalHeatTransferInitialCondition> setup_heat_transfer_initial_condition_ex(boundary_Robin_ex);

    // Overall Body definition
    BodyRegionByParticle epdm_body(diffusion_body, makeShared<MultiPolygonShape>(createEPDMBody()));
    BodyRegionByParticle panel_body(diffusion_body, makeShared<MultiPolygonShape>(createPanelBody()));
    BodyRegionByParticle polyamide_body(diffusion_body, makeShared<MultiPolygonShape>(createPolyamideBody()));
    BodyRegionByParticle ac_body1(diffusion_body, makeShared<MultiPolygonShape>(createACBody1()));
    BodyRegionByParticle ac_body2(diffusion_body, makeShared<MultiPolygonShape>(createACBody2()));
    BodyRegionByParticle ac_body3(diffusion_body, makeShared<MultiPolygonShape>(createACBody3()));
    BodyRegionByParticle ac_body4(diffusion_body, makeShared<MultiPolygonShape>(createACBody4()));
    BodyRegionByParticle ac_body5(diffusion_body, makeShared<MultiPolygonShape>(createACBody5()));
    BodyRegionByParticle ac_body6(diffusion_body, makeShared<MultiPolygonShape>(createACBody6()));
    BodyRegionByParticle ac_body7(diffusion_body, makeShared<MultiPolygonShape>(createACBody7()));
    BodyRegionByParticle ac_open_body1(diffusion_body, makeShared<MultiPolygonShape>(createACOpenBody1()));
    BodyRegionByParticle decreased_convection_body(boundary_Robin_in, makeShared<MultiPolygonShape>(createDecreasedInternalConvectionBody()));

    // Diffusion coefficient initialization
    SimpleDynamics<LocalDiffusivityDefinition> epdm_diffusivity(epdm_body, epdm_cond);
    SimpleDynamics<LocalDiffusivityDefinition> panel_diffusivity(panel_body, pane_cond);
    SimpleDynamics<LocalDiffusivityDefinition> polyamide_diffusivity(polyamide_body, poly_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_diffusivity(ac_body1, ac1_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac2_diffusivity(ac_body2, ac2_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac3_diffusivity(ac_body3, ac3_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac4_diffusivity(ac_body4, ac4_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac5_diffusivity(ac_body5, ac5_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac6_diffusivity(ac_body6, ac6_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac7_diffusivity(ac_body7, ac7_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_open_diffusivity(ac_open_body1, ac1_open_cond);

    SimpleDynamics<LocalConvectionDefinition> decreased_convection_initial_condition(decreased_convection_body, convection_i_decreased);
    SimpleDynamics<LocalHeatTransferConvection> decreased_ht_convection(decreased_convection_body, convection_i_decreased);

    SimpleDynamics<NormalDirectionFromBodyShape> diffusion_body_normal_direction(diffusion_body);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_in(boundary_Robin_in);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_ex(boundary_Robin_ex);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system.real_bodies_);
    ReducedQuantityRecording<QuantitySummation<Real>> write_heat_flux(diffusion_body, "HT_Flux");
    ReducedQuantityRecording<QuantitySummation<Real>> write_diffusion_changing_rate(diffusion_body, "PhiChangeRate");
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
        write_solid_temperature("Phi", temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    setup_diffusion_initial_condition.exec();
    robin_boundary_condition_in.exec();
    robin_boundary_condition_ex.exec();

    setup_heat_transfer_initial_condition_in.exec();
    setup_heat_transfer_initial_condition_ex.exec();

    // thermal conductivity initialization
    epdm_diffusivity.exec();
    panel_diffusivity.exec();
    polyamide_diffusivity.exec();
    ac1_diffusivity.exec();
    ac2_diffusivity.exec();
    ac3_diffusivity.exec();
    ac4_diffusivity.exec();
    ac5_diffusivity.exec();
    ac6_diffusivity.exec();
    ac7_diffusivity.exec();
    ac1_open_diffusivity.exec();

    decreased_convection_initial_condition.exec();
    decreased_ht_convection.exec();

    diffusion_body_normal_direction.exec();
    Robin_normal_direction_in.exec();
    Robin_normal_direction_ex.exec();

    diffusion_body.addBodyStateForRecording<Real>("HT_Flux");
    diffusion_body.addBodyStateForRecording<Real>("ThermalConductivity");

    boundary_Robin_in.addBodyStateForRecording<Real>("Convection");
    boundary_Robin_in.addBodyStateForRecording<Real>("HT_Convection");
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 0.02; // for steady state
    Real End_Time = T0;
    Real Observe_time = 0.01 * End_Time;
    Real Output_Time = 0.1 * End_Time;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_solid_temperature.writeToFile();
    write_diffusion_changing_rate.writeToFile();
    write_heat_flux.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
                if (ite % 500 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }

                temperature_relaxation.exec(dt);
                /*heat_flux_calculation_ex.exec(dt);*/
                heat_flux_calculation_in.exec(dt); // Absolute value of Heat flux _in and _ex are identical

                ite++;
                dt = get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        write_solid_temperature.writeToFile(ite);
        write_diffusion_changing_rate.writeToFile(ite);
        write_heat_flux.writeToFile(ite);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << GlobalStaticVariables::physical_time_ << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_solid_temperature.generateDataBase(1.0e-3, 1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_solid_temperature.testResult();
    }

    return 0;
}
