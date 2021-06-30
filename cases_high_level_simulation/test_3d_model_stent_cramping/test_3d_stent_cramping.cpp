/**
 * @file test_3d_stent_cramping.cpp.cpp
 * @brief This is the example of a stent cramping
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#include "sphinxsys.h"
#include "solid_structural_simulation_class.h"
using namespace SPH;

int main()
{	
	/** INPUT PARAMETERS */
	/*SCALE AND RESOLUTION */
	Real scale_stl = 0.001;
	Real resolution_stent = 0.6;
	Real resolution_cylinder = 1;
	Real resolution_vessel_cylinder = 1;

	/* SIMULATION TIME */
	Real end_time_cramping = 0.05;
	Real end_time_cramping_position = end_time_cramping * 0.8;
	Real end_time_simulation = 0.1;
	Real end_time_position = end_time_simulation * 0.8;

	Real stent_r_outer = 25 + resolution_stent / 2;
	Real stent_r_outer_end = 20;
	Real cylinder_r_inner = 30 - resolution_cylinder / 2;
	Real driving_scale = stent_r_outer_end / cylinder_r_inner;
	
	/* MATERIAL PARAMETERS */
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	// 1e7 is already stable for the stent, the higher the stiffness, the more stable the stent is
	// for development we can use 1e6
	Real Youngs_modulus = 1e6;
	Real Youngs_modulus_vessel = 1e5;
	Real physical_viscosity = 200;

	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	//string relative_input_path = "C:/SPHinXsys_Virtonomy/SPHinXsys/cases_high_level_simulation/test_3d_model_stent_cramping/input/"; //for Windows, define full path
	string stent_stl = "stent.stl";
	string cylinder_stl = "cylinder.stl";
	string vessel_cylinder_stl = "vessel_cylinder.stl";

	vector<string> imported_stl_list = { stent_stl, cylinder_stl, vessel_cylinder_stl };
	vector<Vec3d> translation_list = { Vec3d(0), Vec3d(0), Vec3d(0) };
	vector<Real> resolution_list = { resolution_stent, resolution_cylinder, resolution_vessel_cylinder };

	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	NeoHookeanSolid material_vessel = NeoHookeanSolid(rho_0, Youngs_modulus_vessel, poisson); //material model used for vessels
	vector<LinearElasticSolid> material_model_list = { material, material, material_vessel };

	vector<array<int, 2>> contacting_bodies_list = { {0, 1} };

	/** INPUT DECLERATION */
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		contacting_bodies_list,
	};
	input.position_scale_solid_body_tuple_ = { PositionScaleSolidBodyTuple(1, 0.0, end_time_cramping_position, driving_scale),
												PositionScaleSolidBodyTuple(1, end_time_cramping_position, end_time_simulation, driving_scale) };

	Vecd final_position_cylinder = Vec3d(0, 0, 150)*scale_stl;
	input.position_solid_body_tuple_ = { PositionSolidBodyTuple(1, end_time_cramping, end_time_simulation * 0.8, final_position_cylinder ),
											PositionSolidBodyTuple(1, end_time_simulation * 0.8, end_time_simulation, final_position_cylinder ) };

	// spring damper constraint
	Real spring_stiffness = 500;
	Real damping_ratio = 0.1;
	input.spring_damper_tuple_ = { SpringDamperTuple(2, Vec3d(spring_stiffness), damping_ratio)};

	/** SIMULATION MODEL */
	StructuralSimulation sim(input);
	/** START SIMULATION */
	sim.RunSimulation(end_time_simulation);

	return 0;
}
