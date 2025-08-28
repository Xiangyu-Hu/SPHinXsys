It’s now time to look at our first example. 
Here we’ll introduce features as we go. 
In the next chapter we’ll step back and talk more about the SPHinXsys applications in general.

=======================
Example 1: 2D dam break
=======================

The following program creates a system for the dam break problem: 
a water block is released and moves under 
the action of gravity within a water tank.
It simulates the behavior of this system over a time interval 
in which the water wave impacts the tank wall and produces splashes.

.. code-block:: cpp

	/**
	 * @file 	Dambreak.cpp
	 * @brief 	2D dambreak example.
	 * @details This is the one of the basic test cases, also the first case for
	 * 			understanding SPH method for fluid simulation.
	 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
	 */
	#include "sphinxsys.h" //SPHinXsys Library.
	using namespace SPH;   //Namespace cite here.
	//----------------------------------------------------------------------
	//	Basic geometry parameters and numerical setup.
	//----------------------------------------------------------------------
	Real DL = 5.366;					/**< Tank length. */
	Real DH = 5.366;					/**< Tank height. */
	Real LL = 2.0;						/**< Liquid column length. */
	Real LH = 1.0;						/**< Liquid column height. */
	Real particle_spacing_ref = 0.025;	/**< Initial reference particle spacing. */
	Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
	//----------------------------------------------------------------------
	//	Material properties of the fluid.
	//----------------------------------------------------------------------
	Real rho0_f = 1.0;						 /**< Reference density of fluid. */
	Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
	Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
	Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
	//----------------------------------------------------------------------
	//	Geometric shapes used in this case.
	//----------------------------------------------------------------------
	std::vector<Vecd> water_block_shape{
		Vecd(0.0, 0.0), Vecd(0.0, LH), Vecd(LL, LH), Vecd(LL, 0.0), Vecd(0.0, 0.0)};
	std::vector<Vecd> outer_wall_shape{
		Vecd(-BW, -BW), Vecd(-BW, DH + BW), Vecd(DL + BW, DH + BW), Vecd(DL + BW, -BW), Vecd(-BW, -BW)};
	std::vector<Vecd> inner_wall_shape{
		Vecd(0.0, 0.0), Vecd(0.0, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(0.0, 0.0)};
	//----------------------------------------------------------------------
	//	Fluid body with cases-dependent geometries (ComplexShape).
	//----------------------------------------------------------------------
	class WaterBlock : public FluidBody
	{
	public:
		WaterBlock(SPHSystem &sph_system, const std::string &body_name)
			: FluidBody(sph_system, body_name)
		{
			MultiPolygon multi_polygon;
			multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
			body_shape_.add<MultiPolygonShape>(multi_polygon);
		}
	};
	//----------------------------------------------------------------------
	//	Wall boundary body with cases-dependent geometries.
	//----------------------------------------------------------------------
	class WallBoundary : public SolidBody
	{
	public:
		WallBoundary(SPHSystem &sph_system, const std::string &body_name)
			: SolidBody(sph_system, body_name)
		{
			MultiPolygon multi_polygon;
			multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
			multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

			body_shape_.add<MultiPolygonShape>(multi_polygon);
		}
	};
	//----------------------------------------------------------------------
	//	Observer particle generator.
	//----------------------------------------------------------------------
	class FluidObserverParticleGenerator : public ParticleGeneratorDirect
	{
	public:
		FluidObserverParticleGenerator() : ParticleGeneratorDirect()
		{
			positions_volumes_.push_back(std::make_pair(Vecd(DL, 0.2), 0.0));
		}
	};
	//----------------------------------------------------------------------
	//	Main program starts here.
	//----------------------------------------------------------------------
	int main(int ac, char *av[])
	{
		//----------------------------------------------------------------------
		//	Build up the environment of a SPHSystem.
		//----------------------------------------------------------------------
		SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
		sph_system.handleCommandlineOptions(ac, av);
		/** Tag for computation from restart files. 0: not from restart files. */
		system.restart_step_ = 0;
		/** I/O environment. */
		In_Output in_output(sph_system);
		//----------------------------------------------------------------------
		//	Creating body, materials and particles.
		//----------------------------------------------------------------------
		WaterBlock water_block(sph_system, "WaterBody");
		FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f));

		WallBoundary wall_boundary(sph_system, "Wall");
		SolidParticles wall_particles(wall_boundary);

		ProbeBody fluid_observer(sph_system, "Fluidobserver");
		ObserverParticles observer_particles(fluid_observer, makeShared<FluidObserverParticleGenerator>());
		//----------------------------------------------------------------------
		//	Define body relation map.
		//	The contact map gives the topological connections between the bodies.
		//	Basically the the range of bodies to build neighbor particle lists.
		//----------------------------------------------------------------------
		ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
		BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
		//----------------------------------------------------------------------
		//	Define the main numerical methods used in the simulation.
		//	Note that there may be data dependence on the constructors of these methods.
		//----------------------------------------------------------------------
		Gravity gravity(Vecd(0.0, -gravity_g));
		TimeStepInitialization fluid_step_initialization(water_block, gravity);
		fluid_dynamics::DensitySummationFreeSurfaceComplex fluid_density_by_summation(water_block_complex);
		fluid_dynamics::AdvectionTimeStepSize fluid_advection_time_step(water_block, U_max);
		fluid_dynamics::AcousticTimeStepSize fluid_acoustic_time_step(water_block);
		fluid_dynamics::PressureRelaxationRiemannWithWall fluid_pressure_relaxation(water_block_complex);
		fluid_dynamics::DensityRelaxationRiemannWithWall fluid_density_relaxation(water_block_complex);
		//----------------------------------------------------------------------
		//	Define the methods for I/O operations, observations
		//	and regression tests of the simulation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
		RestartIO restart_io(in_output, sph_system.real_bodies_);
		RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
			write_water_mechanical_energy(in_output, water_block, gravity);
		RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
			write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
		//----------------------------------------------------------------------
		//	Prepare the simulation with cell linked list, configuration
		//	and case specified initial condition if necessary.
		//----------------------------------------------------------------------
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		wall_particles.initializeNormalDirectionFromBodyShape();
		//----------------------------------------------------------------------
		//	Load restart file if necessary.
		//----------------------------------------------------------------------
		if (sph_system.restart_step_ != 0)
		{
			GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
		}
		//----------------------------------------------------------------------
		//	Setup for time-stepping control
		//----------------------------------------------------------------------
		size_t number_of_iterations = sph_system.restart_step_;
		int screen_output_interval = 100;
		int observation_sample_interval = screen_output_interval * 2;
		int restart_output_interval = screen_output_interval * 10;
		Real End_Time = 20.0; /**< End time. */
		Real D_Time = 0.1;	  /**< Time stamps for output of body states. */
		Real dt = 0.0;		  /**< Default acoustic time step sizes. */
		//----------------------------------------------------------------------
		//	Statistics for CPU time
		//----------------------------------------------------------------------
		tick_count t1 = tick_count::now();
		tick_count::interval_t interval;
		tick_count::interval_t interval_computing_time_step;
		tick_count::interval_t interval_computing_fluid_pressure_relaxation;
		tick_count::interval_t interval_updating_configuration;
		tick_count time_instance;
		//----------------------------------------------------------------------
		//	First output before the main loop.
		//----------------------------------------------------------------------
		body_states_recording.writeToFile();
		write_water_mechanical_energy.writeToFile(number_of_iterations);
		write_recorded_water_pressure.writeToFile(number_of_iterations);
		//----------------------------------------------------------------------
		//	Main loop starts here.
		//----------------------------------------------------------------------
		while (GlobalStaticVariables::physical_time_ < End_Time)
		{
			Real integration_time = 0.0;
			/** Integrate time (loop) until the next output time. */
			while (integration_time < D_Time)
			{
				/** outer loop for dual-time criteria time-stepping. */
				time_instance = tick_count::now();
				fluid_step_initialization.parallel_exec();
				Real Dt = fluid_advection_time_step.parallel_exec();
				fluid_density_by_summation.parallel_exec();
				interval_computing_time_step += tick_count::now() - time_instance;

				time_instance = tick_count::now();
				Real relaxation_time = 0.0;
				while (relaxation_time < Dt)
				{
					/** inner loop for dual-time criteria time-stepping.  */
					fluid_pressure_relaxation.parallel_exec(dt);
					fluid_density_relaxation.parallel_exec(dt);
					dt = fluid_acoustic_time_step.parallel_exec();
					relaxation_time += dt;
					integration_time += dt;
					GlobalStaticVariables::physical_time_ += dt;
				}
				interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;

				/** screen output, write body reduced values and restart files  */
				if (number_of_iterations % screen_output_interval == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
							  << GlobalStaticVariables::physical_time_
							  << "	Dt = " << Dt << "	dt = " << dt << "\n";

					if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.restart_step_)
					{
						write_water_mechanical_energy.writeToFile(number_of_iterations);
						write_recorded_water_pressure.writeToFile(number_of_iterations);
					}
					if (number_of_iterations % restart_output_interval == 0)
						restart_io.writeToFile(number_of_iterations);
				}
				number_of_iterations++;

				/** Update cell linked list and configuration. */
				time_instance = tick_count::now();
				water_block.updateCellLinkedList();
				water_block_complex.updateConfiguration();
				fluid_observer_contact.updateConfiguration();
				interval_updating_configuration += tick_count::now() - time_instance;
			}

			tick_count t2 = tick_count::now();
			body_states_recording.writeToFile();
			tick_count t3 = tick_count::now();
			interval += t3 - t2;
		}
		tick_count t4 = tick_count::now();

		tick_count::interval_t tt;
		tt = t4 - t1 - interval;
		std::cout << "Total wall time for computation: " << tt.seconds()
				  << " seconds." << std::endl;
		std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
				  << interval_computing_time_step.seconds() << "\n";
		std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
				  << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
		std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
				  << interval_updating_configuration.seconds() << "\n";

		write_water_mechanical_energy.newResultTest();
		write_recorded_water_pressure.newResultTest();

		return 0;
	};


If you run the test_2d_dambreak correctly, 
you should see a new folder :code:`output` is created.
In the :code:`output` folder, you can see particle state files start with :code:`SPHBody`, :code:`Fluidobserver_Pressure_0.dat`, 
and :code:`WaterBody_TotalMechanicalEnergy_0.dat` which is the global information file.
In the visualization software Paraview you can produces the particle distribution as shown in the following figure. 

.. figure:: ../figures/dambreak.png
   :width: 500 px
   :align: center

   An snapshot of the particle distribution in the dam break problem


Let’s go through the program line by line and see how it works. 
It begins with the include statement:

.. code-block:: cpp

	/**
	 * @file 	Dambreak.cpp
	 * @brief 	2D dambreak example.
	 * @details This is the one of the basic test cases, also the first case for
	 * 			understanding SPH method for fluid simulation.
	 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
	 */
	#include "sphinxsys.h" //SPHinXsys Library.


That gets us all the declarations we need to write a SPHinXsys-using application.

Next we import the :code:`SPH` namespace, 
which includes nearly all of the symbols used by SPHinXsys:

.. code-block:: cpp

	using namespace SPH;   // Namespace cite here.


Now, we provide the parameters for geometric modeling.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Basic geometry parameters and numerical setup.
	//----------------------------------------------------------------------
	Real DL = 5.366;					/**< Tank length. */
	Real DH = 5.366;					/**< Tank height. */
	Real LL = 2.0;						/**< Liquid column length. */
	Real LH = 1.0;						/**< Liquid column height. */
	Real particle_spacing_ref = 0.025;	/**< Initial reference particle spacing. */
	Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));


Here, :code:`particle_spacing_ref` gives the reference initial particle spacing. 
:code:`BW` is the size (thickness) of a wall boundary, which is usually 4 times of particle spacing. 
We give the the coordinates of lower and upper bounds of the domain 
in :code:`system_domain_bounds` 
which will be used as the bounds for a mesh used for building cell linked lists.

We also provide parameters for physical modeling, 
such as material properties of the fluid and physical parameters of the dam break problem.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Material properties of the fluid.
	//----------------------------------------------------------------------
	Real rho0_f = 1.0;						 /**< Reference density of fluid. */
	Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
	Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
	Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */


As we are using a weakly compressible model for imposing incompressibility, 
the maximum speed in the flow and artificial speed of sound are estimated.

Then, we define the realization of :code:`SPHBody` s.
First, the geometric shapes, 
water_block_shape, outer_wall_shape, and inner_wall_shape, 
are defined form the coordinates based on the geometric parameters.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Geometric shapes used in this case.
	//----------------------------------------------------------------------
	std::vector<Vecd> water_block_shape{
		Vecd(0.0, 0.0), Vecd(0.0, LH), Vecd(LL, LH), Vecd(LL, 0.0), Vecd(0.0, 0.0)};
	std::vector<Vecd> outer_wall_shape{
		Vecd(-BW, -BW), Vecd(-BW, DH + BW), Vecd(DL + BW, DH + BW), Vecd(DL + BW, -BW), Vecd(-BW, -BW)};
	std::vector<Vecd> inner_wall_shape{
		Vecd(0.0, 0.0), Vecd(0.0, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(0.0, 0.0)};
	//----------------------------------------------------------------------
	//	Fluid body with cases-dependent geometries (ComplexShape).
	//----------------------------------------------------------------------
	class WaterBlock : public FluidBody
	{
	public:
		WaterBlock(SPHSystem &sph_system, const std::string &body_name)
			: FluidBody(sph_system, body_name)
		{
			MultiPolygon multi_polygon;
			multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
			body_shape_.add<MultiPolygonShape>(multi_polygon);
		}
	};
	//----------------------------------------------------------------------
	//	Wall boundary body with cases-dependent geometries.
	//----------------------------------------------------------------------
	class WallBoundary : public SolidBody
	{
	public:
		WallBoundary(SPHSystem &sph_system, const std::string &body_name)
			: SolidBody(sph_system, body_name)
		{
			MultiPolygon multi_polygon;
			multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
			multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

			body_shape_.add<MultiPolygonShape>(multi_polygon);
		}
	};
	//----------------------------------------------------------------------
	//	Observer particle generator.
	//----------------------------------------------------------------------
	class FluidObserverParticleGenerator : public ParticleGeneratorDirect
	{
	public:
		FluidObserverParticleGenerator() : ParticleGeneratorDirect()
		{
			positions_volumes_.push_back(std::make_pair(Vecd(DL, 0.2), 0.0));
		}
	};


The :code:`WaterBlock` and  :code:`WallBoundary`, 
which are the derived class of :code:`FluidBody` and :code:`SolidBody` respectively, 
are difined with boolean operation, 
such as :code:`add` and :code:`sub`.
The :code:`FluidObserverParticleGenerator` defines the observation body 
through adding the observation point :code:`Vecd(DL, 0.2)`.
The observation body obtains data from the body it is observing at.

After all :code:`SPHBody` s are defined, here comes to the :code:`int main()` function,
in which the application is defined.
In the first part of :code:`main` function, 
an object of :code:`SPHSystem` is created, 
whether the computation begin from restart files is checked, 
and input/output environment is initialized.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	/** I/O environment. */
	InOutput in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock water_block(sph_system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f));

	WallBoundary wall_boundary(sph_system, "Wall");
	SolidParticles wall_particles(wall_boundary);

	ProbeBody fluid_observer(sph_system, "Fluidobserver");
	ObserverParticles observer_particles(fluid_observer, makeShared<FluidObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});


The material, particles and bodies are also created for water block, wall and observer. 
Then, the collection of topological relations,
which specifies for each body the possible interacting bodies, 
are defined. 

After this, the physical dynamics of system is defined 
as method classes in the form of particle discretization.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	TimeStepInitialization fluid_step_initialization(water_block, gravity);
	fluid_dynamics::DensitySummationFreeSurfaceComplex fluid_density_by_summation(water_block_complex);
	fluid_dynamics::AdvectionTimeStepSize fluid_advection_time_step(water_block, U_max);
	fluid_dynamics::AcousticTimeStepSize fluid_acoustic_time_step(water_block);
	fluid_dynamics::PressureRelaxationRiemannWithWall fluid_pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall fluid_density_relaxation(water_block_complex);


First, the external force is defined.
Then, the methods that will used for multiple times are defined.
They are the SPH algorithms for the fluid dynamics and the time step criteria.

After the dynamics, we also define the outputs, 
including the particle states, restart files, global values and observations.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
		write_water_mechanical_energy(in_output, water_block, gravity);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);


The :code:`Vtp` files can be read directly by the open-source visualization code ParaView.
You also have the option to save the files in Tecplot format.
The global information and observation data are written in simple data format 
and are used to check the accuracy of the simulation results in the regression tests. 
The restart files are in :code:`XML` data format. 

Before the computation, 
we need to prepare the simulation with the cell linked list, configuration and the wall normal direction.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();


Finally, the time stepping will almost start. 
However, if the computation begin from restart files. 
The system will be reset.  

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
		fluid_observer_contact.updateConfiguration();
	}


Note that, because the particles have been moved in the previous simulation, 
one need to update the cell-linked list and particle configuration.

The basic control parameter for the simulation is defined,
such as the restart file, output frequency, total simulation time, 
interval for writing output files, etc. 

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0; /**< End time. */
	Real D_Time = 0.1;	  /**< Time stamps for output of body states. */
	Real dt = 0.0;		  /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_fluid_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	write_water_mechanical_energy.writeToFile(number_of_iterations);
	write_recorded_water_pressure.writeToFile(number_of_iterations);


Also the statistic for computation time is initialized and the initial body states and data are outputed.

Here comes the time-stepping loops. 
The computation is carried out with a dual-criteria time-stepping scheme,
as discussed in SPHinXsys's theory section.

.. code-block:: cpp

	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			time_instance = tick_count::now();
			fluid_step_initialization.parallel_exec();
			Real Dt = fluid_advection_time_step.parallel_exec();
			fluid_density_by_summation.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** inner loop for dual-time criteria time-stepping.  */
				fluid_pressure_relaxation.parallel_exec(dt);
				fluid_density_relaxation.parallel_exec(dt);
				dt = fluid_acoustic_time_step.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;

			/** screen output, write body reduced values and restart files  */
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.restart_step_)
				{
					write_water_mechanical_energy.writeToFile(number_of_iterations);
					write_recorded_water_pressure.writeToFile(number_of_iterations);
				}
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			  << interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
			  << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	write_water_mechanical_energy.newResultTest();
	write_recorded_water_pressure.newResultTest();

	return 0;


During the looping outputs are scheduled.
On screen output will be the number of time steps, 
the current physical time, and the advection and acoustic time-step sizes.
After the simulation is terminated, the statistics of computation time 
and the accuracy of the global information and observation data are output on the screen.