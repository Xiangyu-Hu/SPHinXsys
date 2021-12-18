By introducing the Uflyand-Mindlin shell formulation, 
we proposed a simple total Lagrangian SPH method for plate/shell structures,
which is an enhancement of the classical stabilized SPH method commonly used for 3D continua. 
This method allows the modeling of moderately thin structure using only one layer of particles on the plate/shell mid-surface. 
The proposed SPH-based plate/shell modeling is of high-efficiency compared to the classical continuum SPH modeling.

==================================================================================
Example 10: Shell cases: a 2D thin plate
==================================================================================

Here, we consider the even-distributed upward body force applys on a 2D plate, 
when the displacements of the two ends of this plate are fixed. 
Note that the quasi-static analysis is performed here.

First, we provide the parameters for geometric modeling and numerical setup.

.. code-block:: cpp

	/**
	* @file 	2d_plate.cpp
	* @brief 	This is the benchmark test of the 2D shell.
	* @details  We consider the body force applys on a 2D plate.
	* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
	* @version  0.1
	 */
	#include "sphinxsys.h"

	using namespace SPH;

	/**
	 * @brief Basic geometry parameters and numerical setup.
	 */
	Real PL = 10.0;                                          /** Length of the square plate. */
	Real PT = 1.0;                                           /** Thickness of the square plate. */
	Vec2d n_0 = Vec2d(0.0, 1.0);                             /** Initial pseudo-normal. */
	int particle_number = 20;								 /** Particle number in the direction of the length */
	Real resolution_ref = PL / (Real)particle_number;        /** Initial reference particle spacing. */
	int BWD = 1;
	Real BW = resolution_ref * (Real)BWD;                    /** Boundary width, determined by specific layer of boundary particles. */
	/** Domain bounds of the system. */
	BoundingBox system_domain_bounds(Vec2d(-BW, -0.5 * resolution_ref), 
		Vec2d(PL + BW, 0.5 * resolution_ref));

There is one more angular DOF (Degree of Freedom) for 2D thin structure dynamics than that for 2D solid dynamics. 
And we get the pseudo normal direction by rotating the normal direction according to this angular DOF.
So we should define the initial normal direction, :code:`n_0`.
:code:`resolution_ref` gives the reference of initial particle spacing, 
and :code:`BW` gives the boundary width.
:code:`system_domain_bounds` defines the domain of this case.

Then we give the material properties and load applied.

.. code-block:: cpp

	/** For material properties of the solid. */
	Real rho0_s = 1.0; 			                             /** Normalized density. */
	Real Youngs_modulus = 1.3024653e6;	                     /** Normalized Youngs Modulus. */
	Real poisson = 0.3; 			                         /** Poisson ratio. */
	Real physical_viscosity = 200.0;                         /** physical damping used for quasi-static analysis. */
	Real q = 100.0 * Youngs_modulus * 1.0e-4;                /** Total distributed load. */
	Real time_to_full_external_force = 0.1;
	Real gravitational_acceleration = 0.0;

Here we give the physical viscosity for getting quasi-static numerical results. 
:code:`time_to_full_external_force` means the external load increases linearly 
from 0 to the specified value :code:`q` in the time of :code:`time_to_full_external_force` 
and then remains constant.
We ignore the gravity here, so :code:`gravitational_acceleration` is set to 0. 

We can generate particles directly or by Lattice. 
For the piece of code below, 
we define to generate particles directly by giving each particle position and volume. 

.. code-block:: cpp

	class Plate : public ThinStructure
	{
	public:
		Plate(SPHSystem &system, std::string body_name, ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
			: ThinStructure(system, body_name, particle_adaptation, particle_generator)
		{
			// the plate and boundary
			for (int i = 0; i < (particle_number + 2 * BWD); i++)
			{
				Real x = resolution_ref * i - BW + resolution_ref * 0.5;
				body_input_points_volumes_.push_back(std::make_pair(Vecd(x, 0.0), resolution_ref));
			}
		}
	};

Here, we push back each particle's position and volume. 
Note that the volume is :code:`resolution_ref` for 2D shell particles.
And we define initial condition, boundary geometry in following code piece.

.. code-block:: cpp

	/**
	 * application dependent initial condition
	 */
	class PlateDynamicsInitialCondition
		: public thin_structure_dynamics::ShellDynamicsInitialCondition
	{
	public:
		PlateDynamicsInitialCondition(SolidBody *plate)
			: thin_structure_dynamics::ShellDynamicsInitialCondition(plate) {};
	protected:
		void Update(size_t index_i, Real dt) override {
			/** initial pseudo-normal. */
			n_0_[index_i] = n_0;
			n_[index_i] = n_0;
			pseudo_n_[index_i] = n_0_[index_i];
		};
	};

	/** Define the boundary geometry. */
	class BoundaryGeometry : public BodyPartByParticle
	{
	protected:
		virtual void tagBodyPart() override
		{
			BaseParticles* base_particles = body_->base_particles_;
			for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
			{
				if (base_particles->pos_n_[i][0] < 0.0 || base_particles->pos_n_[i][0] > PL)
				{
					tagAParticle(i);
				}
			}
		};
	public:
		BoundaryGeometry(SPHBody *body, std::string body_part_name)
			: BodyPartByParticle(body, body_part_name) {
			tagBodyPart();
		};
		virtual ~BoundaryGeometry() {};
	};

We define the initial pseudo-normal direction :code:`pseudo_n_` of the shell same as the initial normal direction.
Then we use the command :code:`tagAParticle` to tag the boundary particles.
And we define external force, observer body and material properties in following code piece.

.. code-block:: cpp

	/**
	 * define time dependent external force
	 */
	class TimeDependentExternalForce : public Gravity
	{
	public:
		TimeDependentExternalForce(Vecd external_force)
			: Gravity(external_force) {}
		virtual Vecd InducedAcceleration(Vecd& position) override
		{
			Real current_time = GlobalStaticVariables::physical_time_;
			return current_time < time_to_full_external_force ?
				current_time * global_acceleration_ / time_to_full_external_force : global_acceleration_;
		}
	};

	/** Define an observer body. */
	class PlateObserver : public FictitiousBody
	{
	public:
		PlateObserver(SPHSystem &system, std::string body_name)
			: FictitiousBody(system, body_name)
		{
			/** the measuring particle with zero volume */
			body_input_points_volumes_.push_back(std::make_pair(Vecd(0.5 * PL, 0.0), 0.0));
		}
	};

	class PlateMaterial : public LinearElasticSolid
	{
	public:
		PlateMaterial(): LinearElasticSolid()
		{
			rho0_ = rho0_s;
			youngs_modulus_ = Youngs_modulus;
			poisson_ratio_ = poisson;

			assignDerivedMaterialParameters();
		}
	};

The observer body includes only one point, located at the middle of the thin plate.

Here we come to the :code:`int main()` function. 
In the first part of :code:`main` function, 
an object of :code:`SPHSystem` is created, and external force is defined.

.. code-block:: cpp

	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** Define the external force. */
	TimeDependentExternalForce external_force(Vec2d(0.0, q / (PT * rho0_s) - gravitational_acceleration));

Note that the external force applied gives each particle the same acceleration since the load is equally distributed.
The bodies, material and particles are also created in following code piece.

.. code-block:: cpp

	/** Creat a plate body. */
	Plate *plate_body = new Plate(system, "PlateBody", new ParticleAdaptation(1.15, 0), new ParticleGeneratorDirect());
	/** elastic soild material properties */
	PlateMaterial *plate_material = new PlateMaterial();
	/** Creat particles for the elastic body. */
	ShellParticles plate_body_particles(plate_body, plate_material, PT);

When defining :code:`plate_body`, four parameters are inputed.
In :code:`ParticleAdaptation(1.15, 0)`, 1.15 is the smooth length ratio, 
which means the cutoff radius for searching neighbor particls is 2.3 * :code:`resolution_ref`.
And 0 is global refinement level, which means the particle spacing is still :code:`resolution_ref`.
If 0 is changed to 1, the particle spacing will be half :code:`resolution_ref`.
And then the observer body and contact map are defined.

.. code-block:: cpp

	/** Define Observer. */
	PlateObserver *plate_observer = new PlateObserver(system, "PlateObserver");
	BaseParticles observer_particles(plate_observer);

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerBodyRelation* plate_body_inner = new InnerBodyRelation(plate_body);
	ContactBodyRelation* plate_observer_contact = new ContactBodyRelation(plate_observer, { plate_body });

Using class :code:`InnerBodyRelation` means :code:`plate_body_inner` defines the inner data connections.
And using class :code:`ContactBodyRelation` means :code:`plate_observer_contact` 
defines the :code:`palte_observer` has data connections with :code:`plate_body`,
e.g. the :code:`palte_observer` gets data from :code:`plate_body`.
After this, all the physical dynamics are defined in the form of particle discretization.

.. code-block:: cpp

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** Common particle dynamics. */
	InitializeATimeStep 	initialize_external_force(plate_body, &external_force);
	 /** initial condition */
	PlateDynamicsInitialCondition plate_initial_normal(plate_body);
	 /** Corrected strong configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration_in_strong_form(plate_body_inner);
	/** Time step size caclutation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	/** active-pative stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(plate_body, new BoundaryGeometry(plate_body, "BoundaryGeometry"));
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		plate_position_damping(plate_body_inner, 0.5, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		plate_rotation_damping(plate_body_inner, 0.5, "AngularVelocity", physical_viscosity);

First, the applying external force is defined. 
Then comes to the methods, initial condition and correted configuration, that will be executed only once.
Initial condition defines the initial normal and pseudo-normal direction, 
and configuration is corrected to ensure the first-order consistency.
Then, the methods that will used for multiple times are defined. 
They are the SPH algorithms for time step criteria, thin structure dynamics, boundary condition and physical damping.
Note that the time step is dependent on plate thickness and material properties,
and physical damping is applied for quasi-steady analysis.

Before the computation, we also define the outputs, 
including the particle states and obervations.

.. code-block:: cpp

	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToPlt write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<indexVector, Vecd> write_plate_max_displacement("Position", in_output, plate_observer_contact);

The :code:`Plt` files can be read directly by the Tecplot.
You can also save the files in ParaView format by changing :code:`WriteBodyStatesToPlt` to :code:`WriteBodyStatesToVtu`.

The initial conditions, including the cell-linked list and particle configuration, are executed once before the main loop.

.. code-block:: cpp

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	plate_initial_normal.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();

For solid dynamics, we do not change the cell-linked list and particle configuration. 
So they are calculated only once before the simulation.
The basic control parameter for the simulation is defined in the following, 
such as total simulation time 
and interval for writing output files. 

.. code-block:: cpp

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.WriteToFile(0);
	write_plate_max_displacement.WriteToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.5;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

Here, the initial particle states and obervations is written. 
Then we come to the time-stepping loop.

.. code-block:: cpp

	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < output_period)
		{
			if (ite % 1000 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			plate_position_damping.parallel_exec(dt);
			plate_rotation_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

		}
		write_plate_max_displacement.WriteToFile(ite);
		tick_count t2 = tick_count::now();
		write_states.WriteToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;

During the looping outputs are scheduled.
On screen output will be the number of time steps, 
the current physical time and acoustic time-step size.
After the simulation is terminated, the statistics of computation time are outputed to the screen.
Note that the total computation time has excluded the time for writing files.

After the simulation process, one can use the Tecplot to read the result files.
The following figure shows the von Mises stresses of this 2D thin plate.

.. figure:: ../figures/2D_thin_plate.png
   :width: 600 px
   :align: center

   The von Mises stresses of 2D thin plate.



