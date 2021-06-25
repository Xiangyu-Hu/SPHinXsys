/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#include "sim_total_artificial_heart.h"

int main(int argc, char *argv[])
{
	//use simulation collection id to get json file for simulation definittion
	string ticket_id = argv[0];
	//download JSON file and fill in the following struct	

	SimTotalArtificialHeartInput input;
	input.scale_stl = 0.001;
	input.resolution = {8.0,8.0,8.0,8.0,8.0,8.0};
	// in order
	// resolution_tah = 8.0;
	// resolution_aorta = 8.0;
	// resolution_diaphragm = 8.0;
	// resolution_latrium = 8.0;
	// resolution_partery = 8.0;
	// resolution_ratrium = 8.0;
	input.rho_0 = 1000.0;
	input.poisson = 0.35;
	input.Youngs_modulus = 1e5;
	input.Youngs_modulus_tah = 1e6;
	input.physical_viscosity = 200.0;
	input.translation_tah = {0, -200.0, 0};

	string tah_stl = "TAH_basic2_pos.stl";
	string aorta_stl = "Aorta.stl";
	string diaphragm_stl = "Diaphragm.stl";
	string latrium_stl = "LA.stl";
	string partery_stl = "PA.stl";
	string ratrium_stl = "RA.stl";

	input.stls = { tah_stl, aorta_stl, diaphragm_stl, latrium_stl, partery_stl, ratrium_stl };
	input.relative_input_path = "./input/";
	input.contacting_bodies_list = {{0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {1,4}, {3,4}, {4,5}};

	/* DOWNLOAD STLs files at this point */

	// set up the simulation
	SimTotalArtificialHeart simTotalArtificialHeart(input);
	int number_of_steps = 700;
	simTotalArtificialHeart.runSimulationFixedDurationJS(number_of_steps);
	
	return 0;
}