/**
 * @file 	2d_plate.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider point force force apply on a 2D plate.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 */

#include "sphinxsys.h"

using namespace SPH;

//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 10.0;									  /** Length of the square plate. */
Vec2d n_0 = Vec2d(0.0, 1.0);					  /** Pseudo-normal. */
Real thickness = 1.0;							  /** Thickness of the square plate. */
int particle_number = 10;						  /** Particle number in the direction of the length */
Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;									  /** number of boundary particles layers . */
Real BW = resolution_ref * (Real)BWD;			  /** Boundary width, determined by specific layer of boundary particles. */
BoundingBox system_domain_bounds(Vec2d(-BW, -0.5 * resolution_ref), Vec2d(PL + BW, 0.5 * resolution_ref));
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.0)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.0;				   /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;				   /** Poisson ratio. */
//----------------------------------------------------------------------
//	Derived classes used in the case
//----------------------------------------------------------------------
class PlateParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit PlateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			Real position = resolution_ref * i - BW + resolution_ref * 0.5 - 0.5 * PL;
			initializePositionAndVolumetricMeasure(Vecd(position, 0.0), resolution_ref);
			initializeSurfaceProperties(n_0, thickness);
		}
	};
};

class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		body_part_particles_.push_back(index_i);
	};
};

/** Define the controled rotation. */
class ControlRotation : public thin_structure_dynamics::ConstrainShellBodyRegion
{
public:
	ControlRotation(BodyPartByParticle& body_part)
		: ConstrainShellBodyRegion(body_part),
		vel_(particles_->vel_), angular_vel_(particles_->angular_vel_), pos_(particles_->pos_), pos0_(particles_->pos0_) {};
	virtual ~ControlRotation() {};

protected:
	StdLargeVec<Vecd>& vel_, & angular_vel_, & pos_, & pos0_;
	Real ratation_v = Pi;
	void update(size_t index_i, Real dt = 0.0)
	{
		vel_[index_i] = Vecd(-ratation_v * pos_[index_i][1], ratation_v * pos_[index_i][0]);
		angular_vel_[index_i] = Vecd(-ratation_v, 0.0);
	};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody plate_body(system, makeShared<DefaultShape>("PlateBody"));
	plate_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	plate_body.generateParticles<PlateParticleGenerator>();
	plate_body.addBodyStateForRecording<Vecd>("PseudoNormal");
	plate_body.addBodyStateForRecording<Vecd>("Acceleration");
	plate_body.addBodyStateForRecording<Vecd>("AngularAcceleration");
	plate_body.addBodyStateForRecording<Vecd>("AngularVelocity");
	plate_body.addBodyStateForRecording<Matd>("DeformationGradient");
	plate_body.addBodyStateForRecording<Matd>("BendingDeformationGradient");
	plate_body.addBodyStateForRecording<Matd>("DeformationRate");
	plate_body.addBodyStateForRecording<Matd>("BendingDeformationGradientChangeRate");
	plate_body.addBodyStateForRecording<Matd>("GlobalStress");
	plate_body.addBodyStateForRecording<Matd>("GlobalMoment");
	plate_body.addBodyStateForRecording<Vecd>("GlobalShearStress");

	ObserverBody plate_observer(system, "PlateObserver");
	plate_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation plate_body_inner(plate_body);
	ContactRelation plate_observer_contact(plate_observer, {&plate_body});
	//----------------------------------------------------------------------
	//	Define all numerical methods which are used in this case.
	//----------------------------------------------------------------------
	/** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(plate_body_inner);
	/** Time step size. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
	/** stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(plate_body_inner);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	SimpleDynamics<ControlRotation> constrain_holder(boundary_geometry);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_plate_max_displacement("Position", io_environment, plate_observer_contact); //TODO: using ensemble better
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.exec();
	constrain_holder.exec();
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile(0);
	write_plate_max_displacement.writeToFile(0);
	//----------------------------------------------------------------------
	//	Basic control parameters for time stepping.
	//----------------------------------------------------------------------
	int ite = 0;
	Real end_time = 1.0;
	Real output_period = end_time / 1000.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	Main loop of time stepping starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 1000 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			stress_relaxation_first_half.exec(dt);
			//if (GlobalStaticVariables::physical_time_ < 0.5)
			//{
			//	constrain_holder.exec(dt);
			//}
			stress_relaxation_second_half.exec(dt);

			ite++;
			//dt = computing_time_step_size.exec();
			dt = 1.0e-4;

			//if (dt < 1.0e-9) {
			//	write_states.writeToFile(ite);
			//}

			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_plate_max_displacement.writeToFile(ite);
		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
