/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case by testing the relaxation with evolution method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
////----------------------------------------------------------------------
////	Complex geometry with tan=0.5 angle setup.
////----------------------------------------------------------------------
//Real L = 1.0;
//Vec2d circle_center(0.0, L);
//Real circle_radius = 0.5;
//Real resolution_ref = L / 20.0;
//Real BW = resolution_ref * 4.0;
//BoundingBox system_domain_bounds(Vec2d(-BW - 4 * L / 2, -BW - 4 * L / 2), Vec2d(BW + 4 * L / 2, BW + 4 * L / 2));
////----------------------------------------------------------------------
////	Define geometries
////----------------------------------------------------------------------
//std::vector<Vecd> createBlockTriangleShape()
//{
//	std::vector<Vecd> shape;
//	shape.push_back(Vecd(-3 * L / 2, L / 2));
//	shape.push_back(Vecd(L / 2, 4 * L / 2));
//	shape.push_back(Vecd(L / 2, L / 2));
//	shape.push_back(Vecd(3 * L / 2, L / 2));
//	shape.push_back(Vecd(3 * L / 2, -L / 2));
//	shape.push_back(Vecd(-L / 2, -L * 4 / 2));
//	shape.push_back(Vecd(-L / 2, -L / 2));
//	shape.push_back(Vecd(-3 * L / 2, -L / 2));
//	shape.push_back(Vecd(- 3 * L / 2, L / 2));
//	return shape;
//}
//
//class ShapeBody : public MultiPolygonShape
//{
//public:
//	explicit ShapeBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
//	{
//		multi_polygon_.addAPolygon(createBlockTriangleShape(), ShapeBooleanOps::add);
//		multi_polygon_.addACircle(Vecd(-L / 2, -L / 2), circle_radius, 100, ShapeBooleanOps::add);
//		multi_polygon_.addACircle(Vecd( L / 2, L / 2), circle_radius, 100, ShapeBooleanOps::add);
//	}
//};

////----------------------------------------------------------------------
////	Complex geometry with tan=1 angle setup.
////----------------------------------------------------------------------
//Real L = 1.0;
//Vec2d circle_center(0.0, L);
//Real circle_radius = 0.5;
//Real resolution_ref = L / 20.0;
//Real BW = resolution_ref * 4.0;
//BoundingBox system_domain_bounds(Vec2d(-BW - 4 * L / 2, -BW - 4 * L / 2), Vec2d(BW + 4 * L / 2, BW + 4 * L / 2));
////----------------------------------------------------------------------
////	Define geometries
////----------------------------------------------------------------------
//std::vector<Vecd> createBlockTriangleShape()
//{
//	std::vector<Vecd> shape;
//	shape.push_back(Vecd(-L, 0.0));
//	shape.push_back(Vecd(L / 2, 3 * L / 2));
//	shape.push_back(Vecd(L / 2, L / 2));
//	shape.push_back(Vecd(3 * L / 2, L / 2));
//	shape.push_back(Vecd(3 * L / 2, 0.0));
//	shape.push_back(Vecd(-L / 2, -L * 3 / 2));
//	shape.push_back(Vecd(-L / 2, -L / 2));
//	shape.push_back(Vecd(-3 * L / 2, -L / 2));
//	shape.push_back(Vecd(-3 * L / 2, 0.0));
//	shape.push_back(Vecd(-L, 0.0));
//	return shape;
//}
//
//class ShapeBody : public MultiPolygonShape
//{
//public:
//	explicit ShapeBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
//	{
//		multi_polygon_.addAPolygon(createBlockTriangleShape(), ShapeBooleanOps::add);
//		multi_polygon_.addACircle(Vecd(-L / 2, -L / 2), circle_radius, 100, ShapeBooleanOps::add);
//		multi_polygon_.addACircle(Vecd( L / 2, L / 2), circle_radius, 100, ShapeBooleanOps::add);
//	}
//};

////----------------------------------------------------------------------
////	Circle geometry setup.
////----------------------------------------------------------------------
//Vec2d center(0.0, 0.0);	/**< Location of the cylinder center. */
//Real radius = 1.0;		/**< Radius of the cylinder. */
//Real resolution_ref = radius / 20.0;
//Real BW = resolution_ref * 4.0;
//BoundingBox system_domain_bounds(Vec2d(-BW - radius, -BW - radius), Vec2d(BW + radius, BW + radius));
////----------------------------------------------------------------------
////	Define geometries
////----------------------------------------------------------------------
//class ShapeBody : public MultiPolygonShape
//{
//public:
//	explicit ShapeBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
//	{
//		multi_polygon_.addACircle(center, radius, 100, ShapeBooleanOps::add);
//	}
//}; 

//----------------------------------------------------------------------
//	Box geometry setup.
//----------------------------------------------------------------------
Real LL = 2.0;
Real LH = 2.0;
Real resolution_ref = LH / 40.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW-LL/2, -BW-LL/2), Vec2d(BW + LL/2, BW + LH/2));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.0, 0.0);
Vec2d water_block_translation = Vec2d(0.5 * LL, 0.5 * LH);
class ShapeBody : public ComplexShape
{
public:
	explicit ShapeBody(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform2d(water_block_halfsize), water_block_translation);
	}
};

////----------------------------------------------------------------------
////	Basic geometry parameters and numerical setup.
////----------------------------------------------------------------------
//Real L = 1.0;
//Real resolution_ref = L / 20.0;
//Real BW = resolution_ref * 2.0;
//BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(2 * L + BW, 2 * L + BW));
////----------------------------------------------------------------------
////	Define geometries
////----------------------------------------------------------------------
//std::vector<Vecd> createBlockTriangleShape()
//{
//	std::vector<Vecd> shape;
//	shape.push_back(Vecd(0, L / 2));
//	shape.push_back(Vecd(0, 3 * L / 2));
//	shape.push_back(Vecd(L / 2, 2 * L));
//	shape.push_back(Vecd(3 * L / 2, 2 * L));
//	shape.push_back(Vecd(2 * L, 3 * L / 2));
//	shape.push_back(Vecd(2 * L, L / 2));
//	shape.push_back(Vecd(3 * L / 2, 0));
//	shape.push_back(Vecd(L / 2, 0));
//	shape.push_back(Vecd(0, L / 2 + 0));
//	return shape;
//}
//
//class ShapeBody : public MultiPolygonShape
//{
//public:
//	explicit ShapeBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
//	{
//		multi_polygon_.addAPolygon(createBlockTriangleShape(), ShapeBooleanOps::add);
//	}
//};

int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody body(sph_system, makeShared<ShapeBody>("ShapeBody"));
	//body.sph_adaptation_->resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial();
	body.addBodyStateForRecording<Vecd>("Position");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
		: body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation shape_body_inner(body);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle posisiton. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });

		/** An relaxation process include 0th and 1st order consistency. */
		InteractionDynamics<relax_dynamics::CalculateParticleStress> calculate_particle_stress(shape_body_inner, true);
		relax_dynamics::RelaxationStepInner relaxation_0th_inner(shape_body_inner, true);
		relax_dynamics::RelaxationStepImplicitInner relaxation_0th_implicit_inner(shape_body_inner, true);
		relax_dynamics::RelaxationStepByStressInner relaxation_1st_inner(shape_body_inner, true);
		relax_dynamics::RelaxationStepByStressImplicitInner relaxation_1st_implicit_inner(shape_body_inner, true);

		/** Update relaxation residue. */
		InteractionDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_kinetic_energy(shape_body_inner);
		ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>> write_particle_averaged_kinetic_energy(io_environment, body, "particle_kinetic_energy");
		ReducedQuantityRecording<ReduceDynamics<QuantityMaximum<Real>>> write_particle_maximum_kinetic_energy(io_environment, body, "particle_kinetic_energy");

		InteractionDynamics<relax_dynamics::CheckCorrectedZeroOrderConsistency> check_corrected_zero_order_consistency(shape_body_inner, true);
		InteractionDynamics<relax_dynamics::CheckCorrectedFirstOrderConsistency> check_corrected_first_order_consistency(shape_body_inner, true);

		ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>> write_particle_averaged_first_error(io_environment, body, "corrected_first_order_error");
		ReducedQuantityRecording<ReduceDynamics<QuantityMaximum<Real>>> write_particle_maximum_first_error(io_environment, body, "corrected_first_order_error");

		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		write_insert_body_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite_p = 0;
		GlobalStaticVariables::physical_time_ = ite_p;
		while (ite_p <= 10000)
		{
			//calculate_particle_stress.exec();
			//relaxation_1st_implicit_inner.exec();
			relaxation_0th_implicit_inner.exec();

			if (ite_p % 50 == 0) 
			{
				update_kinetic_energy.exec();
				write_particle_averaged_kinetic_energy.writeToFile(ite_p);
				write_particle_maximum_kinetic_energy.writeToFile(ite_p);

				check_corrected_zero_order_consistency.exec();
				check_corrected_first_order_consistency.exec();
				write_particle_averaged_first_error.writeToFile(ite_p);
				write_particle_maximum_first_error.writeToFile(ite_p);
			}

			ite_p += 1;
			GlobalStaticVariables::physical_time_ = ite_p;

			if (ite_p % 50 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_insert_body_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physical relaxation process of body finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
}
