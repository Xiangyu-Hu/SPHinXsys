#include "sphinxsys.h"
#include "solid_structural_simulation_class.h"
using namespace SPH;

int main()
{	
	/** INPUT PARAMETERS */
	Real scale_stl = 0.001 / 4; // diameter of 0.025 m
	Real resolution_mass = 8;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e5;
	Real physical_viscosity = 200;

	// SI setup - designed
	Real rho_0 = 122231;
	Real end_time_simulation = 0.2;
	Real end_time_position = end_time_simulation * 0.8;
	//Vecd final_position_center = Vecd(0, 0, 0.1);

	/** STL IMPORT PARAMETERS */
	string relative_input_path = "./input/"; //path definition for linux
	//string relative_input_path = "C:/SPHinXsys_Virtonomy/SPHinXsys/cases_high_level_simulation/test_3d_ball_position_solid_body/input/";
	vector<std::string> imported_stl_list = { "ball_mass.stl" };
	vector<Vec3d> translation_list = { Vec3d(0) };
	vector<Real> resolution_list = { resolution_mass};

	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	vector<LinearElasticSolid> material_model_list = { material };

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
		{}
	};
	Vecd translation_vector = Vec3d(0,0,200)*scale_stl;
	input.scale_system_boundaries_ = 10;
	input.translation_solid_body_tuple_ = { TranslateSolidBodyTuple(0, 0.0, end_time_position, translation_vector ) };
											
	/** SIMULATION MODEL */
	StructuralSimulation sim (input);
	/** START SIMULATION */
	sim.RunSimulation(end_time_simulation);

	return 0;
}
