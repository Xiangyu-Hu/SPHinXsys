/* ---------------------------------------------------------------------------*
*            SPHinXsys: 3D cantilever           							  *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases, also the  case for            	  *
* understanding SPH method for 3D solid similation under large deformation.   *
* ----------------------------------------------------------------------------*/
/**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------

//for geometry
Real PL = 0.1; 	//cantilever lenght
Real PH = 0.01; //cantilever thickness for thick plate;
Real PW = 0.01;	//cantilever width				
Real SL = 0.03; //depth of the insert

			   //particle spacing, at least three particles
Real particle_spacing_ref = PH / 5.0;
Real BW = particle_spacing_ref * 4; 	//boundary width

										//SimTK geometric modeling resolution
int resolution(20);

//for material properties of the solid
Real rho0_s = 7800.0; //reference density
Real poisson = 0.3; //Poisson ratio
Real Youngs_modulus = 2.1e11;

//------------------------------------------------------------------------------
//define geometry and initial conditions of SPH bodies
//------------------------------------------------------------------------------
/**
* @brief define the cantilever geometry
*/
Geometry *CreateCantilever()
{
	Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_cantilever(0.5 * (PL - SL), 0.0, 0.0);
	Geometry *geometry_cantilever = new Geometry(halfsize_cantilever, resolution,
		translation_cantilever);

	return geometry_cantilever;
}
/**
* @brief define the cantilever holder geometry
*/
Geometry *CreateCantileverHolder()
{
	Vecd halfsize_cantilever_holder(0.5*(SL + BW), 0.5*(PH + BW), 0.5*(PW + BW));
	Vecd translation_cantilever_holder(-0.5*(SL + BW), 0.0, 0.0);
	Geometry *geometry_cantilever_holder = new Geometry(halfsize_cantilever_holder,
		resolution, translation_cantilever_holder);

	return geometry_cantilever_holder;
}
/**
* @brief define the spot for force imposing
*/
Geometry *CreateForceSpot()
{
	Vecd halfsize_cantilever_force_impose(particle_spacing_ref, 0.5*PH, 0.5*PW);
	Vecd translation_cantilever_force_impose(PL, 0.0, 0.0);
	Geometry *geometry_cantilever_force_impose = new Geometry(halfsize_cantilever_force_impose,
		resolution, translation_cantilever_force_impose);

	return geometry_cantilever_force_impose;
}
//define the cantilever body
class Cantilever : public SolidBody
{
public:
	Cantilever(SPHSystem &system, string body_name, ElasticSolid &material,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, material, refinement_level, op)
	{
		//geometry
		body_region_.add_geometry(CreateCantilever(), RegionBooleanOps::add);
		body_region_.add_geometry(CreateCantileverHolder(),	RegionBooleanOps::add);
		//finish the region modeling
		body_region_.done_modeling();
	}
};

/**
* @brief define the Cantilever Holder base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class CantileverHolder : public SolidBodyPart
{
public:
	CantileverHolder(SolidBody *solid_body, string constrianed_region_name)
		: SolidBodyPart(solid_body, constrianed_region_name)
	{
		//geometry
		soild_body_part_region_.add_geometry(CreateCantileverHolder(), RegionBooleanOps::add);
		//finish the region modeling
		soild_body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};

/**
* @brief define the spot where an external force is imposed on
*/
class ForceSpot : public SolidBodyPartForSimbody
{
public:
	ForceSpot(SolidBody *solid_body, string constrianed_region_name, Real solid_body_density)
		: SolidBodyPartForSimbody(solid_body, constrianed_region_name, solid_body_density)
	{
		//geometry
		soild_body_part_region_.add_geometry(CreateForceSpot(), RegionBooleanOps::add);
		//finish the region modeling
		soild_body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};
/**
* @brief imposing a uniform force at the spot
*/
class ImposingUnifromForceOnSpot
	: public solid_dynamics::ImposeExternalForce
{
	Vec3d uniform_accelaeration_;
	Real max_acceleration_;
protected:

	Vec3d GetAcceleration(Vec3d &pos) override
	{
		return uniform_accelaeration_;
	}

	void PrepareConstraint() override
	{
		Real imposing_time_ = 0.5;
		Real time = GlobalStaticVariables::physical_time_;

		uniform_accelaeration_[1]
			= time < 0.5 ? max_acceleration_ * time / imposing_time_
			: max_acceleration_;
	}
public:
	ImposingUnifromForceOnSpot(SolidBody *body, SolidBodyPartForSimbody *body_part)
		: ImposeExternalForce(body, body_part)
	{
		uniform_accelaeration_ = Vec3d(0);
		max_acceleration_ = -1.0e5 / body_part_->body_part_mass_properties_->getMass();
	}
};

//define an observer body
class CantileverObserver : public ObserverLagrangianBody
{
public:
	CantileverObserver(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: ObserverLagrangianBody(system, body_name, refinement_level, op)
	{
		//add observation point
		body_input_points_volumes_.push_back(make_pair(Point(PL, 0.5 * PH, 0.5 * PW), 0.0));
	}
};

//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------

int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(Vecd(-SL - BW, -0.5 * (PL + BW), -0.5 * (SL + BW)),
		Vecd(PL, 0.5 * (PL + BW), 0.5 * (SL + BW)), particle_spacing_ref);

	//Configuration of soild materials
	ElasticSolid solid_material("ElasticSolid", rho0_s, Youngs_modulus, poisson, 0.0);

	//the water block
	Cantilever *cantilever_body =
		new Cantilever(system, "CantileverBody", solid_material, 0, ParticlesGeneratorOps::lattice);
	//creat particles the elastic body
	ElasticSolidParticles cantilever_particles(cantilever_body);

	CantileverObserver *cantilever_observer
		= new CantileverObserver(system, "CantileverObserver", 0, ParticlesGeneratorOps::direct);
	//create observer particles 
	ObserverParticles observer_particles(cantilever_observer);
	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the range of bidies to build neighbor particle lists
	SPHBodyTopology body_topology
		= { { cantilever_body,{} },{ cantilever_observer,{ cantilever_body } } };
	system.SetBodyTopology(&body_topology);

	//setting up the simulation
	system.SetupSPHSimulation();
	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	/** initial condition for the elastic solid bodies */
	solid_dynamics::ElasticSolidDynamicsInitialCondition set_all_gate_particles_at_rest(cantilever_body);
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		cantilever_corrected_configuration_in_strong_form(cantilever_body, {});
	cantilever_corrected_configuration_in_strong_form.parallel_exec();
	//time step size caclutation
	solid_dynamics::GetAcousticTimeStepSize computing_time_step_size(cantilever_body);
	//stress relaxation for the cantilever
	solid_dynamics::StressRelaxationFirstStep
		stress_relaxation_first_step(cantilever_body);
	solid_dynamics::StressRelaxationSecondStep
		stress_relaxation_second_step(cantilever_body);
	//constrain the cantilever holder
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_cantilever_holder(cantilever_body, new CantileverHolder(cantilever_body, "cantileverHolder"));
	//constrain the linear time varying time imposing force
	ImposingUnifromForceOnSpot
		confine_cantilever_force_imposing(cantilever_body, new ForceSpot(cantilever_body, "ForceSpot", rho0_s));
	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	WriteBodyStatesToVtu write_cantilever_states(in_output, system.real_bodies_);
	WriteObservedElasticDisplacement
		write_cantilever_displacement(in_output, cantilever_observer, { cantilever_body});
	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;
	/** apply initial condition */
	set_all_gate_particles_at_rest.exec();
	write_cantilever_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_cantilever_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real T0 = 0.1;
	Real End_Time = 0.5;
	//time step size for oupt file
	Real D_Time = T0 / 100.0;
	Real Dt = 0.01*D_Time;			//default advection time step sizes
	Real dt = 0.0; 					//default accoustic time step sizes

									//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		//integrate time (loop) until the next output time
		while (integeral_time < D_Time) {

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				stress_relaxation_first_step.parallel_exec(dt);
				constrain_cantilever_holder.parallel_exec(dt);
				confine_cantilever_force_imposing.parallel_exec(dt);
				stress_relaxation_second_step.parallel_exec(dt);

				ite++;
				dt = computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		write_cantilever_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

		tick_count t2 = tick_count::now();
		write_cantilever_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
