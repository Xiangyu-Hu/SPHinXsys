/**
 * @file 	test_3d_thin_plate.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider the body force applied on a quasi-static square plate.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.ijnonlinmec.2014.04.009, doi.org/10.1201/9780849384165
 */
#include "sphinxsys.h"
#include <pybind11/pybind11.h>
#include <string>
namespace py = pybind11;
using namespace SPH;   // Namespace cite here.

/**
 * @brief Basic geometry parameters and numerical setup.
 */
class Parameter
{
protected:
	Real PL = 10.0;									  /** Length of the square plate. */
	Real PH = 10.0;									  /** Width of the square plate. */
	Real PT = 1.0;									  /** Thickness of the square plate. */
	Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);				  /** Pseudo-normal. */
	int particle_number = 40;						  /** Particle number in the direction of the length */
	Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
	int BWD = 1;									  /** Width of the boundary layer measured by number of particles. */
	Real BW = resolution_ref * (Real)BWD;			  /** Boundary width, determined by specific layer of boundary particles. */

	/** For material properties of the solid. */
	Real rho0_s = 1.0;				   /** Normalized density. */
	Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
	Real poisson = 0.3;				   /** Poisson ratio. */
	Real physical_viscosity = 200.0;   /** physical damping, here we choose the same value as numerical viscosity. */

	//Real q = 100.0 * Youngs_modulus * 1.0e-4; /** Total distributed load. */
	Real q = Youngs_modulus * 1.0e-4; /** Total distributed load. */
	Real time_to_full_external_force = 0.1;

	Real gravitational_acceleration = 0.009646;
};

/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public SurfaceParticleGenerator, public Parameter
{
public:
	explicit PlateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			for (int j = 0; j < (particle_number + 2 * BWD); j++)
			{
				Real x = resolution_ref * i - BW + resolution_ref * 0.5;
				Real y = resolution_ref * j - BW + resolution_ref * 0.5;
				initializePositionAndVolumetricMeasure(Vecd(x, y, 0.0), resolution_ref * resolution_ref);
				initializeSurfaceProperties(n_0, PT);
			}
		}
	}
};
/** Define the boundary geometry. */
class BoundaryGeometryParallelToXAxis : public BodyPartByParticle, public Parameter
{
public:
	BoundaryGeometryParallelToXAxis(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometryParallelToXAxis::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometryParallelToXAxis() {};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][1] < 0.0 || base_particles_.pos_[index_i][1] > PH)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
class BoundaryGeometryParallelToYAxis : public BodyPartByParticle, public Parameter
{
public:
	BoundaryGeometryParallelToYAxis(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometryParallelToYAxis::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometryParallelToYAxis() {};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][0] < 0.0 || base_particles_.pos_[index_i][0] > PL)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity, public Parameter
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force
				   ? current_time * global_acceleration_ / time_to_full_external_force
				   : global_acceleration_;
	}
};
//----------------------------------------------------------------------
//  Define system, geometry, material, particles and all other things.
//----------------------------------------------------------------------
class PreSettingCase : public Parameter
{
protected:
	/** Domain bounds of the system. */
	BoundingBox system_domain_bounds;
	// Observer location
	StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.5 * PH, 0.0)};
	/** Setup the system. */
	SPHSystem system;

	/** create a plate body. */
	SolidBody plate_body;
	
	/** Define Observer. */
	ObserverBody plate_observer;
	
public:
	PreSettingCase():
		system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
									Vec3d(PL + BW, PH + BW, 0.5 * resolution_ref)),
		system(system_domain_bounds, resolution_ref),
		plate_body(system, makeShared<DefaultShape>("PlateBody")),
		plate_observer(system, "PlateObserver")
	{
		//----------------------------------------------------------------------
		//	Creating bodies with corresponding materials and particles.
		//----------------------------------------------------------------------
		plate_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
		plate_body.generateParticles<PlateParticleGenerator>();

		plate_observer.defineParticlesAndMaterial();
		plate_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	}
};
//----------------------------------------------------------------------
//  Define environment.
//----------------------------------------------------------------------
class Environment : public PreSettingCase
{
protected:
	Real loading_factor;
	/** Set body contact map
		 *  The contact map gives the data connections between the bodies
		 *  basically the the range of bodies to build neighbor particle lists
		 */
	InnerRelation plate_body_inner;
	ContactRelation plate_observer_contact;

	/** Common particle dynamics. */
	SimpleDynamics<TimeStepInitialization> initialize_external_force;
	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration;
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size;
	/** active-passive stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half;
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half;
	/** Constrain the Boundary. */
	BoundaryGeometryParallelToXAxis boundary_geometry_x;
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegionAlongAxis>
		constrain_holder_x;
	BoundaryGeometryParallelToYAxis boundary_geometry_y;
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegionAlongAxis>
		constrain_holder_y;
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
		plate_position_damping;
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
		plate_rotation_damping;
	/** Output */
	IOEnvironment io_environment;
	BodyStatesRecordingToVtp write_states;
	ObservedQuantityRecording<Vecd> write_plate_max_displacement;
	
	
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;

public:
	explicit Environment(Real loading_factor_new) :
		PreSettingCase(),
		loading_factor(loading_factor_new),
		plate_body_inner(plate_body),
		plate_observer_contact(plate_observer, {&plate_body}),
		initialize_external_force(plate_body, 
			makeShared<TimeDependentExternalForce>(Vec3d(0.0, 0.0, q*loading_factor / (PT * rho0_s) - gravitational_acceleration))),
		corrected_configuration(plate_body_inner),
		computing_time_step_size(plate_body),
		stress_relaxation_first_half(plate_body_inner),
		stress_relaxation_second_half(plate_body_inner),
		boundary_geometry_x(plate_body, "BoundaryGeometryParallelToXAxis"),
		constrain_holder_x(boundary_geometry_x, 0),
		boundary_geometry_y(plate_body, "BoundaryGeometryParallelToYAxis"),
		constrain_holder_y(boundary_geometry_y, 1),
		plate_position_damping(0.5, plate_body_inner, "Velocity", physical_viscosity),
		plate_rotation_damping(0.5, plate_body_inner, "AngularVelocity", physical_viscosity),
		io_environment(system),
		write_states(io_environment, system.real_bodies_),
		write_plate_max_displacement("Position", io_environment, plate_observer_contact)
	{	
		std::cout<<"Running simulation for loading factor = " << loading_factor <<"\n";
		/** Apply initial condition. */
		system.initializeSystemCellLinkedLists();
		system.initializeSystemConfigurations();
		corrected_configuration.exec();
		
		write_states.writeToFile(0);
		write_plate_max_displacement.writeToFile(0);
	}
	
	virtual ~Environment() {};	
	//----------------------------------------------------------------------
	//	For ctest.
	//----------------------------------------------------------------------
	int cmakeTest()
	{
		return 1;
	}

	/**
	 *  The main program
	 */
	void runCase()
	{

		/** Set the starting time.
	 	* From here the time stepping begins.
	 	*/
		GlobalStaticVariables::physical_time_ = 0.0;
		/** Setup physical parameters. */
		int ite = 0;
		Real end_time = 0.8;
		Real output_period = end_time / 100.0;
		Real dt = 0.0;

		/**
		 * Main loop
		 */
		while (GlobalStaticVariables::physical_time_ < end_time)
		{
			Real integral_time = 0.0;
			while (integral_time < output_period)
			{
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt: "
							<< dt << "\n";
				}
				initialize_external_force.exec(dt);
				stress_relaxation_first_half.exec(dt);
				constrain_holder_x.exec(dt);
				constrain_holder_y.exec(dt);
				plate_position_damping.exec(dt);
				plate_rotation_damping.exec(dt);
				constrain_holder_x.exec(dt);
				constrain_holder_y.exec(dt);
				stress_relaxation_second_half.exec(dt);

				ite++;
				dt = computing_time_step_size.exec();
				integral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_plate_max_displacement.writeToFile(ite);

			TickCount t2 = TickCount::now();
			TickCount t3 = TickCount::now();
			interval += t3 - t2;
		}
		TickCount t4 = TickCount::now();

		write_states.writeToFile(1);

		TimeInterval tt;
		tt = t4 - t1 - interval;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	}
};

PYBIND11_MODULE(test_3d_thin_plate_python, m)
{
	py::class_<Environment>(m, "thin_plate_from_sph_cpp")
		.def(py::init<const float&>())
		.def("CmakeTest", &Environment::cmakeTest)
		.def("RunCase", &Environment::runCase);
}
