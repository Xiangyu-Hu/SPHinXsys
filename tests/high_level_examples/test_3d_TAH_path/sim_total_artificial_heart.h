/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#ifndef SIM_TOTAL_ARTIFICIAL_HEART_H
#define SIM_TOTAL_ARTIFICIAL_HEART_H

#include "sphinxsys.h"
#include "solid_structural_simulation_class.h"
using namespace SPH;

struct SimTotalArtificialHeartInput
{
	std::vector<std::string> material_model_name;
	double scale_stl;
	std::vector<double> resolution;
	double rho_0;
	double poisson;
	double Youngs_modulus;
	double Youngs_modulus_tah;
	double physical_viscosity;	
	std::array<double, 3> translation_tah;
	std::vector<std::string> stls;
	std::string relative_input_path;
	std::vector<std::array<int, 2>> contacting_bodies_list;
};

class SimTotalArtificialHeart
{
public:
	SimTotalArtificialHeart(SimTotalArtificialHeartInput &input)
	{
		std::vector<string> imported_stl_list(std::begin(input.stls), std::end(input.stls));
		std::vector<Vec3d> translation_list = {Vec3d(input.translation_tah[0], input.translation_tah[1],
													 input.translation_tah[2]),
											   Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0)};
		std::vector<Real> resolution_list = input.resolution;

		LinearElasticSolid material_tah = LinearElasticSolid(input.rho_0, input.Youngs_modulus_tah, input.poisson);
		NeoHookeanSolid material_vessel = NeoHookeanSolid(input.rho_0, input.Youngs_modulus, input.poisson);
		std::vector<LinearElasticSolid> material_model_list = {material_tah,
															   material_vessel, material_vessel, material_vessel,
															   material_vessel, material_vessel};

		/** INPUT DECLERATION */
		StructuralSimulationInput inputStructuralSim{
			input.relative_input_path,
			imported_stl_list,
			input.scale_stl,
			translation_list,
			resolution_list,
			material_model_list,
			input.physical_viscosity,
			input.contacting_bodies_list};
		inputStructuralSim.non_zero_gravity_ = std::vector<GravityPair>{GravityPair(0, Vec3d(0.0, 45.0, 0.0))}; // gravity for TAH
		inputStructuralSim.spring_damper_tuple_ = {SpringDamperTuple(1, Vec3d(0.1, 0.1, 0.1), 0.01),
												   SpringDamperTuple(2, Vec3d(0.1, 0.1, 0.1), 0.01),
												   SpringDamperTuple(3, Vec3d(0.1, 0.1, 0.1), 0.01),
												   SpringDamperTuple(4, Vec3d(0.1, 0.1, 0.1), 0.01),
												   SpringDamperTuple(5, Vec3d(0.1, 0.1, 0.1), 0.01)};

		/** SIMULATION MODEL */
		sim.reset(new StructuralSimulation(inputStructuralSim));
	};

	~SimTotalArtificialHeart(){};

public: //C++ Backend functions
	void runCompleteSimulation(double endTime) { sim->runSimulation(SPH::Real(endTime)); };

public: //WASM functions
	void runSimulationFixedDurationJS(int number_of_steps) { sim->runSimulationFixedDurationJS(number_of_steps); };

private:
	std::unique_ptr<StructuralSimulation> sim;
};

#endif //SIM_TOTAL_ARTIFICIAL_HEART_H