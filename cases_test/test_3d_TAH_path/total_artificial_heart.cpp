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

	//download STL file here for c++ application

	SimTotalArtificialHeart simTotalArtificialHeart;
	simTotalArtificialHeart.initSimulationJS();
	for (int step=0;step<5;step++){
		simTotalArtificialHeart.runSimulationFixedDurationJS(0.1*step);
	}
	return 0;
}