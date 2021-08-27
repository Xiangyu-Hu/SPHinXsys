/**
 * @file test_3d_bernoulli_beam.cpp
 * @brief Bernoulli beam for validation
 * @author Bence Rochlitz, Virtonomy
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real beam_length = 0.2;
Real cross_section_side = 0.01;
Real fixed_length = 0.01;
Real resolution_ref = cross_section_side / 6.0;		/**< Initial particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-fixed_length, 0, 0), Vecd(beam_length, cross_section_side, cross_section_side));
/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 6.45e3; // Nitinol
Real poisson = 0.3;
Real Youngs_modulus = 5e8;
Real physical_viscosity = Youngs_modulus / 100;

/** Define the geometry. */
TriangleMeshShape* CreateBeam()
{
	Vecd halfsize_beam(0.5 * (beam_length + fixed_length), 0.5 * cross_section_side, 0.5 * cross_section_side);
	Vecd translation_beam(0.5 * (beam_length - fixed_length), 0.5 * cross_section_side, 0.5 * cross_section_side);
	TriangleMeshShape* geometry_beam = new TriangleMeshShape(halfsize_beam, resolution, translation_beam);
	return geometry_beam;
}
/** Define the holder geometry. */
TriangleMeshShape* CreateHolder()
{
	Vecd halfsize_holder(0.5*fixed_length, 0.5*cross_section_side , 0.5*cross_section_side);
	Vecd translation_holder(-0.5*fixed_length, 0.5 * cross_section_side, 0.5 * cross_section_side);
	TriangleMeshShape* geometry_holder = new TriangleMeshShape(halfsize_holder, resolution, translation_holder);
	return geometry_holder;
}
/** Define the beam body. */
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addTriangleMeshShape(CreateBeam(), ShapeBooleanOps::add);
		body_shape_->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);
	}
};
/**
* @brief define the Holder base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class Holder : public BodyPartByParticle
{
public:
	Holder(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);
		tagBodyPart();
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Creat a beam body, corresponding material, particles and reaction model. */
	Beam beam_body(system, "BeamBody");
	LinearElasticSolid material(rho_0, Youngs_modulus, poisson);
	ElasticSolidParticles beam_particles(&beam_body, &material);
	/** topology */
	BodyRelationInner beam_body_inner(&beam_body);
	//-------- common particle dynamics ----------------------------------------
	TimeStepInitialization 	initialize_gravity(&beam_body);
	Real end_time = 0.15;
	Real pressure = 1e3;
	StdVec<array<Real, 2>> pressure_over_time = {
		{0.0, 0.0},
		{end_time * 0.1, pressure},
		{end_time, pressure }
	};
	beam_particles.initializeNormalDirectionFromGeometry();
	solid_dynamics::UpdateElasticNormalDirection update_normals(&beam_body);
	solid_dynamics::SurfacePressureFromSource surface_pressure(&beam_body, Vec3d(0.1, 0.5 * cross_section_side, 0.1), pressure_over_time);
	
	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form(&beam_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(&beam_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(&beam_body_inner);
	/** Setup the damping stress, if you know what you are doing. */
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(&beam_body_inner);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(&beam_body, new Holder(&beam_body, "Holder"));
	DampingWithRandomChoice<DampingBySplittingInner<indexVector, Vec3d>> muscle_damping(&beam_body_inner, 0.1, "Velocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_states(in_output, system.real_bodies_);

	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();

	corrected_configuration_in_strong_form.parallel_exec();
	write_states.writeToFile(0);
	/** Setup physical parameters. */
	int ite = 0;
	Real output_period = end_time / 100.0;		
	Real dt = 0.0; 
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period) 
		{
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			update_normals.parallel_exec();
			initialize_gravity.parallel_exec(); // gravity force
			surface_pressure.parallel_exec();

			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
