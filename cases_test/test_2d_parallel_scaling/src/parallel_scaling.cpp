/**
 * @file 	parallel_scaling.cpp
 * @brief 	testing paraellel scaling.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding parallel computing using TBB.
 * @author 	Xiangyu Hu
 * @version 0.1
 */
 /**
  * @brief 	SPHinXsys Library.
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
Real particle_spacing_ref = 0.0125; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0;					/**< Gravity force of fluid. */
Real U_max = 2.0*sqrt(gravity_g*LH);		/**< Characteristic velocity. */
Real c_f = 10.0* U_max;					/**< Reference sound speed. */
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &sph_system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: FluidBody(sph_system, body_name, refinement_level, op)
	{
		/** Geomtry definition. */
		std::vector<Point> water_block_shape;
		water_block_shape.push_back(Point(0.0, 0.0));
		water_block_shape.push_back(Point(0.0, LH));
		water_block_shape.push_back(Point(LL, LH));
		water_block_shape.push_back(Point(LL, 0.0));
		water_block_shape.push_back(Point(0.0, 0.0));
		body_shape_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		/** Basic material parameters*/
		rho_0_ = rho0_f;
		c_0_ = c_f;

		/** Compute the derived material parameters*/
		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(sph_system, body_name, refinement_level, op)
	{
		/** Geomtry definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));
		body_shape_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(0.0, 0.0));
		inner_wall_shape.push_back(Point(0.0, DH));
		inner_wall_shape.push_back(Point(DL, DH));
		inner_wall_shape.push_back(Point(DL, 0.0));
		inner_wall_shape.push_back(Point(0.0, 0.0));
		body_shape_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/**
 * @brief 	Fluid observer body definition.
 */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &sph_system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(sph_system, body_name, refinement_level, 1.3, op)
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
	SPHSystem sph_system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref, 4);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock* water_block
		= new WaterBlock(sph_system, "WaterBody", 0, ParticlesGeneratorOps::lattice);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary* wall_boundary
		= new WallBoundary(sph_system, "Wall", 0, ParticlesGeneratorOps::lattice);
	SolidParticles 					solid_particles(wall_boundary);
	/**
	 * @brief 	Particle and body creation of fluid observer.
	 */
	FluidObserver* fluid_observer
		= new FluidObserver(sph_system, "Fluidobserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles 	observer_particles(fluid_observer);

	/** topology */
	SPHBodyComplexRelation* water_block_complex_relation = new SPHBodyComplexRelation(water_block, { wall_boundary });
	SPHBodyComplexRelation* wall_complex_relation = new SPHBodyComplexRelation(wall_boundary, {});
	SPHBodyContactRelation* fluid_observer_contact_relation = new SPHBodyContactRelation(fluid_observer, { water_block });

	/**
 * @brief 	Define all numerical methods which are used in this case.
 */
 /** Define external force. */
	Gravity 							gravity(Vecd(0.0, -gravity_g));
	/**
	 * @brief 	Methods used only once.
	 */
	 /** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_complex_relation);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummationFreeSurface 		update_fluid_density(water_block_complex_relation);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize 			get_fluid_advection_time_step_size(water_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationFirstHalfRiemann
		pressure_relaxation_first_half(water_block_complex_relation);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block_complex_relation);

	/** Pre-simulation*/
	sph_system.InitializeSystemCellLinkedLists();
	sph_system.InitializeSystemConfigurations();
	get_wall_normal.exec();

	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();

	tick_count time_instance = tick_count::now();

	int screen_output_interval = 1000;
	for (int i = 0; i != 10000; ++i) {

//		pressure_relaxation_first_half.parallel_exec();
//		pressure_relaxation_second_half.parallel_exec();
//		get_fluid_time_step_size.parallel_exec();
//		update_fluid_density.parallel_exec();

		/** Update cell linked list and configuration. */
//		water_block->UpdateCellLinkedList();
		water_block_complex_relation->updateConfiguration();
//		fluid_observer_contact_relation->updateConfiguration();
		if (i % screen_output_interval == 0)
			cout << fixed << setprecision(9) << "N=" << i <<  "\n";
	}

	tick_count::interval_t computing_interval = tick_count::now() - time_instance;
	cout << fixed << setprecision(9) << "interval_updating_configuration = "
		<< computing_interval.seconds() << "\n";

	return 0;
}
