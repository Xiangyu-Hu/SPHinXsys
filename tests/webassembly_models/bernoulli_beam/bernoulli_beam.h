/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#ifndef SIM_TOTAL_ARTIFICIAL_HEART_H
#define SIM_TOTAL_ARTIFICIAL_HEART_H

#include "sphinxsys.h"
#include "structural_simulation_class.h"
using namespace SPH;

#ifdef __EMSCRIPTEN__
#include <emscripten/val.h>
#endif

struct BernoulliBeamInput
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
	StlList stls;
	std::string relative_input_path;
	StdVec<IndexVector> contacting_bodies_list;
};

StructuralSimulationInput createSimulationInput(const BernoulliBeamInput& input, std::shared_ptr<LinearElasticSolid> material_tah, std::shared_ptr<NeoHookeanSolid> material_vessel)
{
	StlList imported_stl_list = input.stls;
	std::vector<Vec3d> translation_list = {Vec3d(input.translation_tah[0], input.translation_tah[1],
													input.translation_tah[2]),
											Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0), Vec3d(0)};
	std::vector<Real> resolution_list = input.resolution;

	std::vector<shared_ptr<LinearElasticSolid>> material_model_list = {
		material_tah,
		material_vessel, material_vessel, material_vessel,
		material_vessel, material_vessel
	};

	/** INPUT DECLERATION */
	StructuralSimulationInput inputStructuralSim = {
		input.relative_input_path,
		imported_stl_list,
		input.scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		StdVec<Real>(6, input.physical_viscosity),
		input.contacting_bodies_list};
	inputStructuralSim.non_zero_gravity_ = std::vector<GravityPair>{GravityPair(0, Vec3d(0.0, 45.0, 0.0))}; // gravity for TAH
	
	return inputStructuralSim;
}

class BernoulliBeam
{
public:	
	SimTotalArtificialHeart(const BernoulliBeamInput& input):
	material_tah_(make_shared<LinearElasticSolid>(input.rho_0, input.Youngs_modulus_tah, input.poisson)),
	material_vessel_(make_shared<NeoHookeanSolid>(input.rho_0, input.Youngs_modulus, input.poisson))
	{
		sim.reset(new StructuralSimulation(createSimulationInput(input, material_tah_, material_vessel_)));
	}

public: //C++ Backend functions
	void runCompleteSimulation(double endTime) { sim->runSimulation(SPH::Real(endTime)); };
	
private:
	std::shared_ptr<LinearElasticSolid> material_tah_;
	std::shared_ptr<NeoHookeanSolid> material_vessel_;
	std::unique_ptr<StructuralSimulation> sim;
};

class BernoulliBeamJS
{
public:
	BernoulliBeamJS(const BernoulliBeamInput& input):
	material_tah_(make_shared<LinearElasticSolid>(input.rho_0, input.Youngs_modulus_tah, input.poisson)),
	material_vessel_(make_shared<NeoHookeanSolid>(input.rho_0, input.Youngs_modulus, input.poisson))
	{
		sim_js_.reset(new StructuralSimulationJS(createSimulationInput(input, material_tah_, material_vessel_)));
	}

	void runSimulation(int number_of_steps) 
	{ 
		try
		{
			sim_js_->runSimulationFixedDuration(number_of_steps);
		}
		catch (const std::exception& e)
		{
			if (on_error_)
				on_error_(e.what());
		}
	}

	BodyStatesRecordingToVtuString::VtuStringData getVtuData() const { return sim_js_->getVtuData(); }

	void onError(emscripten::val on_error) { on_error_ = [on_error](const std::string& error_message) { on_error(error_message); }; }

private:
	std::shared_ptr<LinearElasticSolid> material_tah_;
	std::shared_ptr<NeoHookeanSolid> material_vessel_;
	std::unique_ptr<StructuralSimulationJS> sim_js_;
	std::function<void(const std::string&)> on_error_;
};

#endif //SIM_TOTAL_ARTIFICIAL_HEART_H