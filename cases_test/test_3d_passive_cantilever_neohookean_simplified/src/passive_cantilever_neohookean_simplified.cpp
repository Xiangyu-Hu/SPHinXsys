/**
 * @file passive_cantilever_neohookean.cpp
 * @brief This is the example of myocardium with simple neohookean tissue model 
 * @author Bence Rochlitz, Chi Zhang  and Xiangyu Hu
 * @version 0.1.0
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
#include "cantilever_class.h"

/** Name space. */
using namespace SPH;

/**
 *  The main program
 */
int main()
{	
	//Simulation setup
	PassiveCantilever *passive_cantilever = new PassiveCantilever()
	passive_cantilever->initialize_simulation();

	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_displacement("Displacement", in_output, myocardium_observer_contact);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 1.0;
	Real output_period = end_time / 100.0;		

	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period) 
		{
			if (ite % 100 == 0) {
				cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< passive_cantilever->dt << "\n";
			}

			passive_cantilever->perform_simulation_step();
			
			GlobalStaticVariables::physical_time_ += passive_cantilever->dt;
			integration_time += passive_cantilever->dt;
			ite++;
		}
		write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t2 = tick_count::now();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
