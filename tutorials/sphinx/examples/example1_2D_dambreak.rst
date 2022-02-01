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
in which the water wave impacts the tank wall and producing splashes.

.. code-block:: cpp

	/**
	* @file 	Dambreak.cpp
	* @brief 	2D dambreak exaple.
	* @details This is the one of the basic test cases, also the first case for
	* 			understanding SPH method for fluid similation.
	* @author 	Luhui Han, Chi Zhang and Xiangyu Hu
	* @version 0.1
	*/
	/**
	* @brief SPHinXsys Library.
	*/
	#include "sphinxsys.h"
	/**
	* @brief Namespace cite here.
	*/
	using namespace SPH;
	/**
	* @brief Basic geometry parameters and numerical setup.
	*/
	Real DL = 5.366; 						/**< Tank length. */
	Real DH = 5.366; 						/**< Tank height. */
	Real LL = 2.0; 							/**< Liquid colume length. */
	Real LH = 1.0; 							/**< Liquid colume height. */
	Real particle_spacing_ref = 0.025; 		/**< Initial reference particle spacing. */
	Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
	/**
	* @brief Material properties of the fluid.
	*/
	Real rho0_f = 1.0;						/**< Reference density of fluid. */
	Real gravity_g = 1.0;					/**< Gravity force of fluid. */
	Real U_f = 2.0*sqrt(gravity_g*LH);		/**< Characteristic velocity. */
	Real c_f = 10.0*U_f;					/**< Reference sound speed. */

	Real initial_pressure = 0.0;			/**< Initial pressure field. */
	Vec2d intial_velocity(0.0, 0.0);		/**< Initial velocity field. */
	/**
	* @brief 	Fluid body definition.
	*/
	class WaterBlock : public WeaklyCompressibleFluidBody
	{
	public:
		WaterBlock(SPHSystem &system, string body_name,
		WeaklyCompressibleFluid* material,
		WeaklyCompressibleFluidParticles
		&weakly_compressible_fluid_particles, int refinement_level, ParticlesGeneratorOps op)
		: WeaklyCompressibleFluidBody(system, body_name, material,
			weakly_compressible_fluid_particles, refinement_level, op)
		{
			/** Geomerty definition. */
			std::vector<Point> water_block_shape;
			water_block_shape.push_back(Point(0.0, 0.0));
			water_block_shape.push_back(Point(0.0, LH));
			water_block_shape.push_back(Point(LL, LH));
			water_block_shape.push_back(Point(LL, 0.0));
			water_block_shape.push_back(Point(0.0, 0.0));
			Geometry *water_block_geometry = new Geometry(water_block_shape);
			body_region_.add_geometry(water_block_geometry, RegionBooleanOps::add);

			body_region_.done_modeling();
		}
		/** Initialize every fluid particle data. */
		void InitialCondition()
		{
			for (int i = 0; i < number_of_particles_; ++i) {
				BaseParticleData  &base_particle_data_i = weakly_compressible_fluid_particles_.base_particle_data_[i];
				WeaklyCompressibleFluidParticleData 
					&fluid_data_i
					= weakly_compressible_fluid_particles_.fluid_data_[i];

				fluid_data_i.p_ = initial_pressure;
				base_particle_data_i.vel_n_ = intial_velocity;
				base_particle_data_i.dvel_dt_(0);
				fluid_data_i.rho_0_
					= material_->ReinitializeRho(initial_pressure);
				fluid_data_i.rho_n_ = fluid_data_i.rho_0_;
				fluid_data_i.mass_
					= fluid_data_i.rho_0_*base_particle_data_i.Vol_;
			}
		}
	};
	/**
	* @brief 	Wall boundary body definition.
	*/
	class WallBoundary : public SolidBody
	{
	public:
		WallBoundary(SPHSystem &system, string body_name,
			SolidBodyParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op)
			: SolidBody(system, body_name, solid_particles, refinement_level, op)
		{
			/** Geomerty definition. */
			std::vector<Point> outer_wall_shape;
			outer_wall_shape.push_back(Point(-BW, -BW));
			outer_wall_shape.push_back(Point(-BW, DH + BW));
			outer_wall_shape.push_back(Point(DL + BW, DH + BW));
			outer_wall_shape.push_back(Point(DL + BW, -BW));
			outer_wall_shape.push_back(Point(-BW, -BW));
			Geometry *outer_wall_geometry = new Geometry(outer_wall_shape);
			body_region_.add_geometry(outer_wall_geometry, RegionBooleanOps::add);

			std::vector<Point> inner_wall_shape;
			inner_wall_shape.push_back(Point(0.0, 0.0));
			inner_wall_shape.push_back(Point(0.0, DH));
			inner_wall_shape.push_back(Point(DL, DH));
			inner_wall_shape.push_back(Point(DL, 0.0));
			inner_wall_shape.push_back(Point(0.0, 0.0));
			Geometry *inner_wall_geometry = new Geometry(inner_wall_shape);
			body_region_.add_geometry(inner_wall_geometry, RegionBooleanOps::sub);

			body_region_.done_modeling();
		}
		/** Initialize every wallboundary particle data. */
		void InitialCondition()
		{
			for (int i = 0; i < solid_particles_.number_of_particles_; ++i) {
				BaseParticleData &base_particle_data_i
					= solid_particles_.base_particle_data_[i];
				SolidBodyParticleData &solid_body_data_i
					= solid_particles_.solid_body_data_[i];
		
				base_particle_data_i.vel_n_ = intial_velocity;
				Vec2d zero(0);
				base_particle_data_i.dvel_dt_ = zero;
				solid_body_data_i.vel_ave_ = zero;
				solid_body_data_i.dvel_dt_ave_ = zero;
			}
		}
	};
	/**
	* @brief 	Fluid observer body definition.
	*/
	class FluidObserver : public ObserverBody
	{
	public:
		FluidObserver(SPHSystem &system, string body_name,
			ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op)
			: ObserverBody(system, body_name, observer_particles, refinement_level, op)
		{
			body_input_points_volumes_.push_back(make_pair(Point(DL, 0.2), 0.0));
		}
	};
	/**
	* @brief 	Main program starts here.
	*/
	int main()
	{
		/**
		* @brief Build up -- a SPHSystem --
		*/
		SPHSystem system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
		/** Set the starting time. */
		GlobalStaticVariables::physical_time_ = 0.0;
		/** Tag for computation from restart files. 0: not from restart files. */
		system.restart_step_ = 0;
		/**
		* @brief Material property, partilces and body creation of fluid.
		*/
		WeaklyCompressibleFluid 			fluid("Water", rho0_f, c_f, mu_f, k_f);
		WeaklyCompressibleFluidParticles 	fluid_particles("WaterBody");
		WaterBlock *water_block = new WaterBlock(system, "WaterBody", &fluid,
			fluid_particles, 0, ParticlesGeneratorOps::lattice);
		/**
		* @brief 	Particle and body creation of wall boundary.
		*/
		SolidBodyParticles 					solid_particles("Wall");
		WallBoundary *wall_boundary = new WallBoundary(system, "Wall",
			solid_particles, 0, ParticlesGeneratorOps::lattice);
		/**
		* @brief 	Particle and body creation of fluid observer.
		*/
		ObserverParticles 					observer_particles("Fluidobserver");
		FluidObserver *fluid_observer = new FluidObserver(system, "Fluidobserver",
			observer_particles, 0, ParticlesGeneratorOps::direct);
		/**
		* @brief 	Body contact map.
		* @details The contact map gives the data conntections between the bodies.
		* 			Basically the the rang of bidies to build neighbor particle lists.
		*/
		SPHBodyTopology 	body_topology = { { water_block, { wall_boundary } },
												{ wall_boundary, {} },{ fluid_observer,{ water_block} } };
		system.SetBodyTopology(&body_topology);
		/**
		* @brief 	Simulation set up.
		*/
		system.SetupSPHSimulation();
		/**
		* @brief 	Define all numerical methods which are used in this case.
		*/
		/** Define external force. */
		Gravity 							gravity(Vecd(0.0, -gravity_g));
		/**
		* @brief 	Methods used only once.
		*/
		/** Initialize normal direction of the wall boundary. */
		solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
		get_wall_normal.exec();
		/** Obtain the initial number density of fluid. */
		fluid_dynamics::InitialNumberDensity 		
			fluid_initial_number_density(water_block, { wall_boundary });
		fluid_initial_number_density.exec();
		/**
		* @brief 	Methods used for time stepping.
		*/
		/** Initialize particle acceleration. */
		InitializeOtherAccelerations 	initialize_fluid_acceleration(water_block);
		/** Add particle acceleration due to gravity force. */
		AddGravityAcceleration 			add_fluid_gravity(water_block, &gravity);
		/**
		* @brief 	Algorithms of fluid dynamics.
		*/
		/** Wvaluation of density by summation approach. */
		fluid_dynamics::DensityBySummationFreeSurface 		
			update_fluid_desnity(water_block, { wall_boundary });
		/** Time step size without considering sound wave speed. */
		fluid_dynamics::FluidAdvectionTimeStepSize 			
			get_fluid_adevction_time_step_size(water_block, U_f);
		/** Time step size with considering sound wave speed. */
		fluid_dynamics::WeaklyCompressibleFluidTimeStepSize get_fluid_time_step_size(water_block);
		/** Pressure relaxation algorithm by using verlet time stepping. */
		fluid_dynamics::PressureRelaxationVerletFreeSurface 
			pressure_relaxation(water_block, { wall_boundary }, &gravity);
		/**
		* @brief 	Methods used for updating data structure.
		*/
		/** Update the cell linked list of bodies when neccessary. */
		ParticleDynamicsCellLinkedList			update_cell_linked_list(water_block);
		/** Update the configuration of bodies when neccessary. */
		ParticleDynamicsConfiguration 			update_particle_configuration(water_block);
		/** Update the interact configuration of bodies when neccessary. */
		ParticleDynamicsInteractionConfiguration 	
			update_observer_interact_configuration(fluid_observer, { water_block });
		/**
		* @brief Output.
		*/
		Output output(system);
		/** Output the body states. */
		WriteBodyStatesToVtu 		write_body_states(output, system.real_bodies_);
		/** Output the body states for restart simulation. */
		WriteRestartFileToXml		write_restart_body_states(output, system.real_bodies_);
		/** Output the mechanical energy of fluid body. */
		WriteWaterMechanicalEnergy 	write_water_mechanical_energy(output, water_block, &gravity);
		/** output the observed data from fluid body. */
		WriteObservedFluidPressure	write_recorded_water_pressure(output, fluid_observer, { water_block });
		/**
		* @brief The time stepping starts here.
		*/
		/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
		if (system.restart_step_ != 0)
		{
			system.ResetSPHSimulationFromRestart();
			update_cell_linked_list.parallel_exec();
			update_particle_configuration.parallel_exec();
		}
		/** Output the start states of bodies. */
		write_body_states
			.WriteToFile(GlobalStaticVariables::physical_time_);
		/** Output the Hydrostatic mechanical energy of fluid. */
		write_water_mechanical_energy
			.WriteToFile(GlobalStaticVariables::physical_time_);
		/**
		* @brief 	Basic parameters.
		*/
		int ite = system.restart_step_;
		int rst_out = 1000;
		Real End_Time = 20.0; 	/**< End time. */
		Real D_Time = 0.1;		/**< Time stamps for output of body states. */
		Real Dt = 0.0;			/**< Default advection time step sizes. */
		Real dt = 0.0; 			/**< Default accoustic time step sizes. */
		/** statistics for computing CPU time. */
		tick_count t1 = tick_count::now();
		tick_count::interval_t interval;
		/** Output global basic parameters. */
		output.WriteCaseSetup(End_Time, D_Time, GlobalStaticVariables::physical_time_);
		/**
		* @brief 	Main loop starts here.
		*/
		while (GlobalStaticVariables::physical_time_ < End_Time)
		{
				Real integeral_time = 0.0;
				/** Integrate time (loop) until the next output time. */
				while (integeral_time < D_Time)
				{
					/** Acceleration due to viscous force and gravity. */
					initialize_fluid_acceleration.parallel_exec();
					add_fluid_gravity.parallel_exec();
					Dt = get_fluid_adevction_time_step_size.parallel_exec();
					update_fluid_desnity.parallel_exec();
					
					/** Dynamics including pressure relaxation. */
					Real relaxation_time = 0.0;
					while (relaxation_time < Dt)
					{
						if (ite % 100 == 0)
						{
							cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_
							<< "	dt: " << dt << "\n";
							if (ite % rst_out == 0)
							write_restart_body_states.WriteToFile(Real(ite));
						}
						pressure_relaxation.parallel_exec(dt);
						
						ite++;
						dt = get_fluid_time_step_size.parallel_exec();
						relaxation_time += dt;
						integeral_time += dt;
						GlobalStaticVariables::physical_time_ += dt;
					
					}
					/** Update cell linked list and configuration. */
					update_cell_linked_list.parallel_exec();
					update_particle_configuration.parallel_exec();
					update_observer_interact_configuration.parallel_exec();
				}
				
				tick_count t2 = tick_count::now();
				write_water_mechanical_energy
					.WriteToFile(GlobalStaticVariables::physical_time_);
				write_body_states
					.WriteToFile(GlobalStaticVariables::physical_time_);
				write_recorded_water_pressure
					.WriteToFile(GlobalStaticVariables::physical_time_);
				tick_count t3 = tick_count::now();
				interval += t3 - t2;
			}
		tick_count t4 = tick_count::now();
		
		tick_count::interval_t tt;
		tt = t4 - t1 - interval;
		cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;
		
		return 0;
	}


Before you can compile and run this program, 
you need to have SPHinXsys installed. 
The installation directory has subdirectories :code:`include`, :code:`lib`, and :code:`bin`. 
Make sure the :code:`include` directory is part of your compiler’s include path, and the :code:`lib` directory is available to the linker. At runtime the shared library directory (lib for Mac and Linux, bin for Windows) must be on the appropriate  :code:`path` environment variable. Exactly how you do this will depend on the compiler and operating system you are using.

If everything is working correctly, 
you should see a new folder :code:`output` is created and particle state files start with :code:`SPHBody`
and :code:`Fluidobserver_fluid_pressure.dat` and :code:`WaterBody_water_mechnical_energy.dat`, 
which are global information files.
In the visualization software Paraview you can produces the particle distribution as shown in the following figure. 

.. figure:: ../figures/dambreak.png
   :width: 500 px
   :align: center

   An snapshot of the particle distribution in the dam break problem


Let’s go through the program line by line and see how it works. 
It begins with the include statements:

.. code-block:: cpp

	/**
	* @file 	Dambreak.cpp
	* @brief 	2D dambreak exaple.
	* @details This is the one of the basic test cases, also the first case for
	* 			understanding SPH method for fluid similation.
	* @author 	Luhui Han, Chi Zhang and Xiangyu Hu
	* @version 0.1
	*/
	/**
	* @brief SPHinXsys Library.
	*/
	#include "sphinxsys.h"


That gets us all the declarations we need to write a SPHinXsys-using application.

Next we import the :code:`SPH` namespace, 
which includes nearly all of the symbols used by SPHinXsys:

.. code-block:: cpp

	/**
	* @brief Namespace cite here.
	*/
	using namespace SPH;


Now, we provide the parameters for geometric modeling.

.. code-block:: cpp

	/**
	* @brief Basic geometry parameters and numerical setup.
	*/
	Real DL = 5.366; 						/**< Tank length. */
	Real DH = 5.366; 						/**< Tank height. */
	Real LL = 2.0; 							/**< Liquid colume length. */
	Real LH = 1.0; 							/**< Liquid colume height. */
	Real particle_spacing_ref = 0.025; 		/**< Initial reference particle spacing. */
	Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */


Here, :code:`particle_spacing_ref` gives 
the reference initial particle spacing for multi-resolution modeling, e.g. for refinement level 0. 
:code:`BW` is the size (thickness) of a wall boundary, which is usually 4 times of particle spacing.

We also provide parameters for physical modeling, 
such as material properties of the fluid and physical parameters of the dam break problem.

.. code-block:: cpp

	/**
	* @brief Material properties of the fluid.
	*/
	Real rho0_f = 1.0;						/**< Reference density of fluid. */
	Real gravity_g = 1.0;					/**< Gravity force of fluid. */
	Real U_f = 2.0*sqrt(gravity_g*LH);		/**< Characteristic velocity. */
	Real c_f = 10.0*U_f;					/**< Reference sound speed. */

	Real initial_pressure = 0.0;			/**< Initial pressure field. */
	Vec2d intial_velocity(0.0, 0.0);		/**< Initial velocity field. */


As we are using a weakly compressible model for imposing incompressibility, 
the maximum speed in the flow and artificial speed of sound are estimated.

Then, we define the realization of :code:`SPHBody` s. 
First, a :code:`WaterBlock`, which is a derived class of :code:`WeaklyCompressibleFluidBody`, 
is defined with constructor parameters, such as material, particles, refinement level 
and the option for particle generator.

.. code-block:: cpp

	/**
	* @brief 	Fluid body definition.
	*/
	class WaterBlock : public WeaklyCompressibleFluidBody
	{
		public:
		WaterBlock(SPHSystem &system, string body_name,
			WeaklyCompressibleFluid* material,
			WeaklyCompressibleFluidParticles
			&weakly_compressible_fluid_particles, int refinement_level, ParticlesGeneratorOps op)
			: WeaklyCompressibleFluidBody(system, body_name, material,
			weakly_compressible_fluid_particles, refinement_level, op)
		{
			/** Geomerty definition. */
			std::vector<Point> water_block_shape;
			water_block_shape.push_back(Point(0.0, 0.0));
			water_block_shape.push_back(Point(0.0, LH));
			water_block_shape.push_back(Point(LL, LH));
			water_block_shape.push_back(Point(LL, 0.0));
			water_block_shape.push_back(Point(0.0, 0.0));
			Geometry *water_block_geometry = new Geometry(water_block_shape);
			body_region_.add_geometry(water_block_geometry, RegionBooleanOps::add);
			
			body_region_.done_modeling();
		}
		/** Initialize every fluid particle data. */
		void InitialCondition()
		{
			for (int i = 0; i < number_of_particles_; ++i) {
				BaseParticleData &base_particle_data_i
					= weakly_compressible_fluid_particles_.base_particle_data_[i];
				WeaklyCompressibleFluidParticleData &fluid_data_i
					= weakly_compressible_fluid_particles_.fluid_data_[i];
				
				fluid_data_i.p_ = initial_pressure;
				base_particle_data_i.vel_n_ = intial_velocity;
				base_particle_data_i.dvel_dt_(0);
				fluid_data_i.rho_0_
					= material_->ReinitializeRho(initial_pressure);
				fluid_data_i.rho_n_ = fluid_data_i.rho_0_;
				fluid_data_i.mass_
					= fluid_data_i.rho_0_*base_particle_data_i.Vol_;
			}
		}
	};


Here, the body geometry is defined from the coordinates 
based on the geometric parameters and binary operations, 
such as :code:`add` and :code:`sub`.
Note that, the initial condition of the :code:`WaterBlock` 
is also given in a member function :code:`void InitialCondition()`.
Similarly, we define the :code:`WallBoundary` and :code:`FluidObserver`.
Note that there is no initial condition for the observation body
as it usually only obtain data from the body it is observing at.

After all :code:`SPHBody` s are defined, here comes to the :code:`int main()` function,
which the application is defined.
In the first part of :code:`main` function, 
an object of :code:`SPHSystem` is created, global physical time initialized,
and whether the computation begin from restart files is checked.

.. code-block:: cpp

	/**
	* @brief Build up -- a SPHSystem --
	*/
	SPHSystem system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	/**
	* @brief Material property, partilces and body creation of fluid.
	*/
	WeaklyCompressibleFluid 			fluid("Water", rho0_f, c_f, mu_f, k_f);
	WeaklyCompressibleFluidParticles 	fluid_particles("WaterBody");
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 
		&fluid, fluid_particles, 0, ParticlesGeneratorOps::lattice);
	/**
	* @brief 	Particle and body creation of wall boundary.
	*/
	SolidBodyParticles 					solid_particles("Wall");
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall",
	solid_particles, 0, ParticlesGeneratorOps::lattice);
	/**
	* @brief 	Particle and body creation of fluid observer.
	*/
	ObserverParticles 					observer_particles("Fluidobserver");
	FluidObserver *fluid_observer = new FluidObserver(system, "Fluidobserver",
	observer_particles, 0, ParticlesGeneratorOps::direct);
	/**
	* @brief 	Body contact map.
	* @details The contact map gives the data conntections between the bodies.
	* 			Basically the the rang of bidies to build neighbor particle lists.
	*/
	SPHBodyTopology 	body_topology = { { water_block, { wall_boundary } },
	{ wall_boundary, {} },{ fluid_observer,{ water_block} } };
	system.SetBodyTopology(&body_topology);
	/**
	* @brief 	Simulation set up.
	*/
	system.SetupSPHSimulation();


Note that the constructor of :code:`SPHSystem` requires the coordinates of 
lower and upper bounds of the domain, which will be used as the bounds 
for a mesh used for building cell linked lists.
The material, particles and bodies are also created for water block, wall and observer. 
Then, the collection of topological relations,
which specifies for each body the possible interacting bodies, 
are defined. The function :code:`SetupSPHSimulation()` creates SPH particles,
builds particle configurations and set initial condition if necessary.

After this, the physical dynamics of system is defined 
as method classes in the form of particle discretization.

.. code-block:: cpp

		/**
	* @brief 	Define all numerical methods which are used in this case.
	*/
	/** Define external force. */
	Gravity 							gravity(Vecd(0.0, -gravity_g));
	/**
	* @brief 	Methods used only once.
	*/
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
	get_wall_normal.exec();
	/** Obtain the initial number density of fluid. */
	fluid_dynamics::InitialNumberDensity 		
	fluid_initial_number_density(water_block, { wall_boundary });
	fluid_initial_number_density.exec();
	/**
	* @brief 	Methods used for time stepping.
	*/
	/** Initialize particle acceleration. */
	InitializeOtherAccelerations 	initialize_fluid_acceleration(water_block);
	/** Add particle acceleration due to gravity force. */
	AddGravityAcceleration 			add_fluid_gravity(water_block, &gravity);
	/**
	* @brief 	Algorithms of fluid dynamics.
	*/
	/** Wvaluation of density by summation approach. */
	fluid_dynamics::DensityBySummationFreeSurface 		
	update_fluid_desnity(water_block, { wall_boundary });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::FluidAdvectionTimeStepSize 			
	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::WeaklyCompressibleFluidTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	fluid_dynamics::PressureRelaxationVerletFreeSurface 
	pressure_relaxation(water_block, { wall_boundary }, &gravity);


First, the external force is defined.
Then comes the methods that will be used only once,
such as computing normal direction of the static wall surface, 
and  the initial particle number density. 
Then, the methods that will used for multiple times are defined.
They are the SPH algorithms for fluid dynamics, time step criteria.

The methods for updating particle configurations will be realized in the following,
including update cell linked list, inner (within the body) 
and contact (with the interacting bodies) neighboring particles.

.. code-block:: cpp

	/**
	* @brief 	Methods used for updating data structure.
	*/
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList			update_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 			update_particle_configuration(water_block);
	/** Update the interact configuration of bodies when neccessary. */
	ParticleDynamicsInteractionConfiguration 	
		update_observer_interact_configuration(fluid_observer, { water_block });


Note that such updating can be specified for a given body for its inner and/or 
all contact configuration or cell-linked list, 
or given pair of bodies for the interact configuration.

Before the computation, we also define the outputs, 
including the particle states, restart files, global values and observations.

.. code-block:: cpp

		/**
	* @brief Output.
	*/
	Output output(system);
	/** Output the body states. */
	WriteBodyStatesToVtu 		write_body_states(output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	WriteRestartFileToXml		write_restart_body_states(output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	WriteWaterMechanicalEnergy 	write_water_mechanical_energy(output, water_block, &gravity);
	/** output the observed data from fluid body. */
	WriteObservedFluidPressure	write_recorded_water_pressure(output, fluid_observer, { water_block });


The :code:`Vtu` files can be read directly by the open-source visualization code ParaView.
You also have the option to save the files in Tecplot format.
The global information and observation data are written in simple data format. 
The restart files are in :code:`XML` data format.

Finally, the time stepping will almost start. 
However, if the computation begin from restart files. 
The system will be reset.  

.. code-block:: cpp

		/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		system.ResetSPHSimulationFromRestart();
		update_cell_linked_list.parallel_exec();
		update_particle_configuration.parallel_exec();
	}
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);


Note that, because the particles have been moved in the previous simulation, 
one need to update the cell-linked list and particle configuration. 
After that, the states from the starting time step will be outputted. 

The basic control parameter for the simulation is defined.
Such as the restart file output frequency, total simulation time 
and interval for writing output files. 

.. code-block:: cpp

		/**
	* @brief 	Basic parameters.
	*/
	int ite = system.restart_step_;
	int rst_out = 1000;
	Real End_Time = 20.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default accoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Output global basic parameters. */
	output.WriteCaseSetup(End_Time, D_Time, GlobalStaticVariables::physical_time_);


Also the statistic for computation time is initialized.
A case setup file will be written as a summary of the case. This file goes together with other output data for later reference.

Here comes the time-stepping loops. 
The computation is carried out with a dual-criteria time-stepping scheme,
as discussed in SPHinXsys's theory section.

.. code-block:: cpp

		/**
	* @brief 	Main loop starts here.
	*/
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_fluid_acceleration.parallel_exec();
			add_fluid_gravity.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				if (ite % 100 == 0)
				{
					cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_
					<< "	dt: " << dt << "\n";
					if (ite % rst_out == 0)
					write_restart_body_states.WriteToFile(Real(ite));
			}
			pressure_relaxation.parallel_exec(dt);
			
			ite++;
			dt = get_fluid_time_step_size.parallel_exec();
			relaxation_time += dt;
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			
			}
			/** Update cell linked list and configuration. */
			update_cell_linked_list.parallel_exec();
			update_particle_configuration.parallel_exec();
			update_observer_interact_configuration.parallel_exec();
		}
		
		tick_count t2 = tick_count::now();
		write_water_mechanical_energy
		.WriteToFile(GlobalStaticVariables::physical_time_);
		write_body_states
		.WriteToFile(GlobalStaticVariables::physical_time_);
		write_recorded_water_pressure
		.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
	<< " seconds." << endl;

	return 0;


During the looping outputs are scheduled.
On screen output will be the number of time steps, 
the current physical time and acoustic time-step size.
After the simulation is terminated, the statistics of computation time are output to the screen.
Note that the total computation time has excluded the time for writing files.


