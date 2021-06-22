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
	input.resolution_tah = 8.0;
	input.resolution_aorta = 8.0;
	input.resolution_diaphragm = 8.0;
	input.resolution_latrium = 8.0;
	input.resolution_partery = 8.0;
	input.resolution_ratrium = 8.0;
	input.rho_0 = 1000.0;
	input.poisson = 0.35;
	input.Youngs_modulus = 1e5;
	input.Youngs_modulus_tah = 1e6;
	input.physical_viscosity = 200.0;
	input.translation_tah = {0, -200.0, 0};
	input.stls = {"TAH_basic2_pos.stl", "Aorta.stl", "Diaphragm.stl", "LA.stl", "PA.stl", "RA.stl"};
	input.relative_input_path = "./input/";
	input.contacting_bodies_list = {std::pair<int, int>(0, 1), std::pair<int, int>(0, 2), std::pair<int, int>(0, 3),
									std::pair<int, int>(0, 4), std::pair<int, int>(0, 5), std::pair<int, int>(1, 4),
									std::pair<int, int>(3, 4), std::pair<int, int>(4, 5)};
	/** CONTACT ORGANS WITH ORGANS*/
	//IndexPair(1, 4); //Aorta with PA
	//IndexPair(2, 5); //Diaphragm with RA
	//IndexPair(3, 4); //LA with PA
	//IndexPair(4, 5); //PA with RA

	//download STLs file

	SimTotalArtificialHeart simTotalArtificialHeart(input);
	for (int step = 0; step < 5; step++)
	{
		simTotalArtificialHeart.runSimulationFixedDurationJS(0.02 * step);
	}
	return 0;
}