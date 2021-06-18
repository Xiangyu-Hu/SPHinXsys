/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#include "sphinxsys.h"
#include "solid_structural_simulation_class.h"
using namespace SPH;

class SimTotalArtificialHeart
{
public:
	SimTotalArtificialHeart()
	{
		//set up struct

		/** INPUT PARAMETERS */
		Real scale_stl = 0.001;
		Real resolution_tah = 8;
		Real resolution_aorta = 8;
		Real resolution_diaphragm = 8;
		Real resolution_latrium = 8;
		Real resolution_partery = 8;
		Real resolution_ratrium = 8;
		Real rho_0 = 1000.0;
		Real poisson = 0.35;
		Real Youngs_modulus = 1e5;
		Real Youngs_modulus_tah = 1e6;
		Real physical_viscosity = 200;
		Real translation_tah_in[] = {0,-200,0};
   		std::string stls[] = { "TAH_basic2_pos.stl", "Aorta.stl", "Diaphragm.stl", "LA.stl", "PA.stl", "RA.stl" };
		/** STL IMPORT PARAMETERS */
		string relative_input_path = "./input/"; //path definition for linux

		Vec3d translation_tah = Vec3d(translation_tah_in[0], translation_tah_in[1], translation_tah_in[2]);
		vector<string> imported_stl_list(std::begin(stls), std::end(stls));
		vector<Vec3d> translation_list = {translation_tah, Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0)};
		vector<Real> resolution_list = {resolution_tah, resolution_aorta, resolution_diaphragm, resolution_latrium, 
			resolution_partery, resolution_ratrium};

		LinearElasticSolid material_tah = LinearElasticSolid(rho_0, Youngs_modulus_tah, poisson);
		NeoHookeanSolid material_vessel = NeoHookeanSolid(rho_0, Youngs_modulus, poisson);
		vector<LinearElasticSolid> material_model_list = {material_tah, 
			material_vessel, material_vessel, material_vessel, material_vessel, material_vessel};

		vector<IndexPair> contacting_bodies_list = {IndexPair(0, 1), IndexPair(0, 2), IndexPair(0, 3), IndexPair(0, 4), 
													IndexPair(0, 5), IndexPair(1, 4), IndexPair(3, 4), IndexPair(4, 5)};
		/** CONTACT ORGANS WITH ORGANS*/
		//IndexPair(1, 4); //Aorta with PA
		//IndexPair(2, 5); //Diaphragm with RA
		//IndexPair(3, 4); //LA with PA
		//IndexPair(4, 5); //PA with RA

		/** INPUT DECLERATION */
		StructuralSimulationInput input{
			relative_input_path,
			imported_stl_list,
			scale_stl,
			translation_list,
			resolution_list,
			material_model_list,
			physical_viscosity,
			contacting_bodies_list};
		input.non_zero_gravity_ = vector<GravityPair>{GravityPair(0, Vec3d(0.0, 45.0, 0.0))}; // gravity for TAH
		input.spring_damper_tuple_ = {SpringDamperTuple(1, Vec3d(0.1, 0.1, 0.1), 0.01),
									  SpringDamperTuple(2, Vec3d(0.1, 0.1, 0.1), 0.01),
									  SpringDamperTuple(3, Vec3d(0.1, 0.1, 0.1), 0.01),
									  SpringDamperTuple(4, Vec3d(0.1, 0.1, 0.1), 0.01),
									  SpringDamperTuple(5, Vec3d(0.1, 0.1, 0.1), 0.01)};

		/** SIMULATION MODEL */
		sim = new StructuralSimulation(&input);
	};
	~SimTotalArtificialHeart()
	{
		delete sim;
	};	

public: //C++ Backend functions
	void initAndRunCompleteSimulation(float endTime) { sim->RunSimulation(SPH::Real(endTime)); };
public: //WASM functions
	void initSimulationJS() { sim->InitSimulationJS(); };
	void runSimulationFixedDurationJS(float duration) { sim->RunSimulationFixedDurationJS(SPH::Real(duration)); };

private:
	StructuralSimulation *sim;
};

int main()
{
	SimTotalArtificialHeart simTotalArtificialHeart;
	simTotalArtificialHeart.initSimulationJS();
	for (int step=0;step<5;step++){
		simTotalArtificialHeart.runSimulationFixedDurationJS(0.1*step);
	}
	return 0;
}