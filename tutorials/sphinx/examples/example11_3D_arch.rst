In Example 10, we explained how to build up a 2D shell case. 
Here we will give you a 3D arch case. 
Different form the previous case, the initial psudo-normal is not same between each particle, 
and dynamic responses are solved.
Here we will emphasise the method differences for the 3D arch compared to the 2D plate.

===================================
Example 11: Shell cases: a 3D arch
===================================

The example of a deep circular arch shell structure is shown in the figure, 
of which two ends are fully clamped and the even-distributed downward body force is applied.

.. figure:: ../figures/3D_arch.png
   :width: 600 px
   :align: center

   The dynamic responses of a deep circular arch shell structure.

First, we provide the parameters for geometric modeling, numerical setup, material properties, etc.

.. code-block:: cpp

	/**
	* @file 	3d_arch.cpp
	* @brief 	This is the benchmark test of the 3D shell.
	* @details  We consider large deformation of a circular arch.
	* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
	* @version  0.1
	 */
	#include "sphinxsys.h"

	using namespace SPH;

	/**
	 * @brief Basic geometry parameters and numerical setup.
	 */
	Real radius = 0.0975;                                    /** Radius of the inner boudary of the cylinder. */
	Real height = 0.02;                                      /** Height of the cylinder. */
	Real thickness = 0.005;                                  /** Thickness of the cylinder. */
	Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
	int particle_number = 10;								 /** Particle number in the height direction. */
	Real particle_spacing_ref = height / (Real)particle_number;
	int particle_number_mid_surface = 2 * radius_mid_surface * Pi * 215.0 / 360.0 / particle_spacing_ref;
	int BWD = 4;
	Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */
	/** Domain bounds of the system. */
	BoundingBox system_domain_bounds(Vec3d(-radius - thickness, 0.0, -radius - thickness),
		Vec3d(radius + thickness, height, radius + thickness));

	/** For material properties of the solid. */
	Real rho0_s = 7.800; 			                         /** Normalized density. */
	Real Youngs_modulus = 210e6;	                         /** Normalized Youngs Modulus. */
	Real poisson = 0.3; 			                         /** Poisson ratio. */

	Real time_to_full_external_force = 0.0;
	Real body_acceleration = -400000;

There is two more angular DOFs (Degrees of Freedom) for 3D thin structure dynamics than that for 3D solid dynamics. 
And we get the pseudo normal direction by rotating the normal direction according to these two angular DOFs.
We do not define the physical viscosity any more because of dynamic analysis.
:code:`time_to_full_external_force`  is equal to 0.0, 
which means the external force remains constant , that is to say, stepping load is applied.

Then we define to generate particles directly by giving each particle position and volume, just like 2D case. 

.. code-block:: cpp

	class Cylinder : public ThinStructure
	{
	public:
		Cylinder(SPHSystem &system, std::string body_name, ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
			: ThinStructure(system, body_name, particle_adaptation, particle_generator)
		{
			// the cylinder and boundary
			for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
			{
				for (int j = 0; j < particle_number; j++)
				{
					Real x = radius_mid_surface * cos(-17.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
					Real y = particle_spacing_ref * j + particle_spacing_ref * 0.5;
					Real z = radius_mid_surface * sin(-17.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
					body_input_points_volumes_.push_back(std::make_pair(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref));
				}
			}
		}
	};

And we define initial condition, boundary geometry, external force, observer body and material properties 
in following code piece.

.. code-block:: cpp

	/**
	 * application dependent initial condition
	 */
	class CylinderDynamicsInitialCondition
		: public thin_structure_dynamics::ShellDynamicsInitialCondition
	{
	public:
		CylinderDynamicsInitialCondition(SolidBody *plate)
			: thin_structure_dynamics::ShellDynamicsInitialCondition(plate) {};
	protected:
		void Update(size_t index_i, Real dt) override {
			/** initial pseudo-normal. */
			n_0_[index_i] = Vec3d(pos_0_[index_i][0] / radius_mid_surface, 0.0, pos_0_[index_i][2] / radius_mid_surface);
			n_[index_i] = n_0_[index_i];
			pseudo_n_[index_i] = n_0_[index_i];
			transformation_matrix_[index_i] = getTransformationMatrix(n_0_[index_i]);
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
				if (base_particles->pos_n_[i][2] < radius_mid_surface * sin(-17.5 / 180.0 * Pi))
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
	class CylinderObserver : public FictitiousBody
	{
	public:
		CylinderObserver(SPHSystem &system, std::string body_name)
			: FictitiousBody(system, body_name)
		{
			/** the measuring particle with zero volume */
			body_input_points_volumes_.push_back(std::make_pair(Vecd(0.0, height / 2.0, radius_mid_surface), 0.0));
		}
	};

	class CylinderMaterial : public LinearElasticSolid
	{
	public:
		CylinderMaterial(): LinearElasticSolid()
		{
			rho0_ = rho0_s;
			youngs_modulus_ = Youngs_modulus;
			poisson_ratio_ = poisson;

			assignDerivedMaterialParameters();
		}
	};

Note that we should also initialize the transformation matrices, 
which will be used for coordinate tranformations from global coordinates to local coordiantes.
The observer body includes only one point, located at the middle of the arch.

Here we come to the :code:`int main()` function. 
In the first part of :code:`main` function, 
an object of :code:`SPHSystem` is created, and external force is defined.

.. code-block:: cpp

	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Define the external force. */
	TimeDependentExternalForce external_force(Vec3d(0.0, 0.0, body_acceleration));

Note that the external force applied gives each particle the same acceleration since the load is equally distributed.
The bodies, material and particles are also created in following code piece.

.. code-block:: cpp

	/** Creat a Cylinder body. */
	Cylinder *cylinder_body = new Cylinder(system, "CylinderBody", new ParticleAdaptation(1.15, 0), new ParticleGeneratorDirect());
	/** elastic soild material properties */
	CylinderMaterial *cylinder_material = new CylinderMaterial();
	/** Creat particles for the elastic body. */
	ShellParticles cylinder_body_particles(cylinder_body, cylinder_material, thickness);

And then the observer body and contact map are defined.

.. code-block:: cpp

	/** Define Observer. */
	CylinderObserver *cylinder_observer = new CylinderObserver(system, "CylinderObserver");
	BaseParticles observer_particles(cylinder_observer);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerBodyRelation* cylinder_body_inner = new InnerBodyRelation(cylinder_body);
	ContactBodyRelation* cylinder_observer_contact = new ContactBodyRelation(cylinder_observer, { cylinder_body });

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
	InitializeATimeStep 	initialize_external_force(cylinder_body, &external_force);
	 /** initial condition */
	CylinderDynamicsInitialCondition cylinder_initial_pseudo_normal(cylinder_body);
	 /** Corrected strong configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration_in_strong_form(cylinder_body_inner);
	/** Time step size caclutation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(cylinder_body);
	/** active-pative stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(cylinder_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(cylinder_body_inner);
	/** Constrain the Boundary. */
	thin_structure_dynamics::ClampConstrainShellBodyRegion
		constrain_holder(cylinder_body_inner, new BoundaryGeometry(cylinder_body, "BoundaryGeometry"));

First, the applying external force is defined. 
Then comes to the methods, initial condition and correted configuration, that will be executed only once.
Initial condition defines the initial normal and pseudo-normal direction, and the transformation matrix, 
and configuration is corrected to ensure the first-order consistency.
Then, the methods that will used for multiple times are defined.

Before the computation, we also define the outputs, 
including the particle states and obervations.

.. code-block:: cpp

	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToPlt write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<indexVector, Vecd>
		write_cylinder_max_displacement("Position", in_output, cylinder_observer_contact);

The initial conditions, including the cell-linked list and particle configuration, are executed once before the main loop.

.. code-block:: cpp

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	cylinder_initial_pseudo_normal.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();

For solid dynamics, we do not change the cell-linked list and particle configuration. 
So they are calculated only once before the simulation.
The basic control parameter for the simulation is defined in the following. 

.. code-block:: cpp

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.WriteToFile(0);
	write_cylinder_max_displacement.WriteToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.005;
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
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			Real check_time = GlobalStaticVariables::physical_time_;

		}
		write_cylinder_max_displacement.WriteToFile(ite);
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

The main function is almost same with that of 2D case, 
except that the physical damping is not included.
