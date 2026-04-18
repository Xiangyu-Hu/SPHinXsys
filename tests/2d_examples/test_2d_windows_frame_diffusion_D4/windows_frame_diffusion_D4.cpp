/*
 * @file 	windows_frame_diffusion_D4.cpp
 * @brief 	Heat diffusion of windows frame based on international standard ISO 10077-2:2012,
 *          Application 4: Wood frame section and insulation panel.
 * @details https://www.comsol.de/model/glazing-influence-on-thermal-performances-of-a-window-16075
 * @author 	Haotian Ji, Dong Wu, Chi Zhang and Xiangyu Hu
 */
#include "windows_frame_diffusion_D4.h"
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
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<MultiPolygonShape>(createOverallStructureBody(), "DiffusionBody"));
    diffusion_body.defineClosure<Solid, LocalIsotropicDiffusion>(
        Solid(), ConstructArgs(diffusion_species_name, wood_cond, epdm_cond));
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody boundary_Robin_in(sph_system, makeShared<MultiPolygonShape>(createInternalAirBody(), "InternalConvectionBoundary"));
    boundary_Robin_in.generateParticles<BaseParticles, Lattice>();

    SolidBody boundary_Robin_ex(sph_system, makeShared<MultiPolygonShape>(createExternalAirBody(), "ExternalConvectionBoundary"));
    boundary_Robin_ex.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of temperature observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation inner_relation(diffusion_body);
    ContactRelation contact_Robin(diffusion_body, {&boundary_Robin_in, &boundary_Robin_ex});
    ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> diffusion_body_normal_direction(diffusion_body);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_in(boundary_Robin_in);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_ex(boundary_Robin_ex);

    // Define body regions
    BodyRegionByParticle epdm_body(diffusion_body, makeShared<MultiPolygonShape>(createEPDMBody()));
    BodyRegionByParticle panel_body(diffusion_body, makeShared<MultiPolygonShape>(createPanelBody()));
    BodyRegionByParticle ac_body1(diffusion_body, makeShared<MultiPolygonShape>(createACBody1()));
    BodyRegionByParticle ac_body2(diffusion_body, makeShared<MultiPolygonShape>(createACBody2()));
    BodyRegionByParticle ac_open_body1(diffusion_body, makeShared<MultiPolygonShape>(createACOpenBody1()));
    // Define diffusion coefficient
    SimpleDynamics<LocalDiffusivityDefinition> epdm_diffusivity(epdm_body, epdm_cond);
    SimpleDynamics<LocalDiffusivityDefinition> panel_diffusivity(panel_body, pane_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_diffusivity(ac_body1, ac1_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac2_diffusivity(ac_body2, ac2_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_open_diffusivity(ac_open_body1, ac1_open_cond);

    DiffusionBodyRelaxation temperature_relaxation(inner_relation, contact_Robin);
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body);

    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_in(boundary_Robin_in);
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_ex(boundary_Robin_ex);
    BodyRegionByParticle decreased_convection_body(boundary_Robin_in, makeShared<MultiPolygonShape>(createDecreasedInternalConvectionBody()));
    SimpleDynamics<LocalConvectionDefinition> decreased_convection_initial_condition(decreased_convection_body, convection_i_decreased);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    ReducedQuantityRecording<Average<QuantitySummation<Real, SPHBody>>>
        write_heat_transfer_from_internal_boundary(diffusion_body, diffusion_species_name + "TransferFromInternalConvectionBoundary");
    ReducedQuantityRecording<Average<QuantitySummation<Real, SPHBody>>>
        write_heat_transfer_from_external_boundary(diffusion_body, diffusion_species_name + "TransferFromExternalConvectionBoundary");
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
        write_solid_temperature(diffusion_species_name, temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    setup_diffusion_initial_condition.exec();
    robin_boundary_condition_in.exec();
    robin_boundary_condition_ex.exec();
    decreased_convection_initial_condition.exec();

    // thermal conductivity initialization
    epdm_diffusivity.exec();
    panel_diffusivity.exec();
    ac1_diffusivity.exec();
    ac2_diffusivity.exec();
    ac1_open_diffusivity.exec();

    diffusion_body_normal_direction.exec();
    Robin_normal_direction_in.exec();
    Robin_normal_direction_ex.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < End_Time)
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
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                temperature_relaxation.exec(dt);

                ite++;
                dt = get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        write_solid_temperature.writeToFile(ite);
        write_heat_transfer_from_internal_boundary.writeToFile(ite);
        write_heat_transfer_from_external_boundary.writeToFile(ite);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << physical_time << " seconds." << std::endl;

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
