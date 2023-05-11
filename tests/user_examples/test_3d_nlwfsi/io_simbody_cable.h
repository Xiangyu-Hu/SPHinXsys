/**
 * @file 	io_simbody_cable.h
 * @brief 	Classes for simbody cable output files.
 * @author	Nicolò Salis
 */

#pragma once

#include "io_base.h"

namespace SPH
{
	/**
	 * @class SimBodyCableStatesIO
	 * @brief base class for write and read SimBody states.
	 */
	template <class CableSpring>
	class SimBodyCableStatesIO
	{
	protected:
		IOEnvironment &io_environment_;
		SimTK::RungeKuttaMersonIntegrator &integ_;
		CableSpring &cable_;

	public:
		SimBodyCableStatesIO(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, CableSpring &cable1)
			: io_environment_(io_environment), integ_(integ), cable_(cable1){};
		virtual ~SimBodyCableStatesIO(){};
	};

	/**
	 * @class WriteSimBodyCableStates
	 * @brief base class for write SimBody states.
	 */
	template <class CableSpring>
	class WriteSimBodyCableStates : public SimBodyCableStatesIO<CableSpring>
	{
	public:
		WriteSimBodyCableStates(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, CableSpring &cable1)
			: SimBodyCableStatesIO<CableSpring>(io_environment, integ, cable1){};
		virtual ~WriteSimBodyCableStates(){};

		virtual void writeToFile(size_t iteration_step) = 0;
	};

	/**
	 * @class ReadSimBodyCableStates
	 * @brief base class for read SimBody states.
	 */
	template <class CableSpring>
	class ReadSimBodyCableStates : public SimBodyCableStatesIO<CableSpring>
	{
	public:
		ReadSimBodyCableStates(IOEnvironment &io_environment, CableSpring *cable1)
			: SimBodyCableStatesIO<CableSpring>(io_environment, cable1){};
		ReadSimBodyCableStates(IOEnvironment &io_environment, StdVec<CableSpring *> cables)
			: SimBodyCableStatesIO<CableSpring>(io_environment, cables){};
		virtual ~ReadSimBodyCableStates(){};

		virtual void readFromFile(size_t iteration_step) = 0;
	};

	/**
	 * @class WriteSimBodyCableData
	 * @brief Write total force acting a single cable element.
	 */
	class WriteSimBodyCableData : public WriteSimBodyCableStates<SimTK::CableSpring>
	{
	protected:
		std::string filefullpath_;

	public:
		WriteSimBodyCableData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::CableSpring &cable1, std::string cableinf);
		virtual ~WriteSimBodyCableData(){};
		virtual void writeToFile(size_t iteration_step = 0) override;
	};
	//=============================================================================================//
	WriteSimBodyCableData::
		WriteSimBodyCableData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::CableSpring &cable1, std::string cableinf)
		: WriteSimBodyCableStates<SimTK::CableSpring>(io_environment, integ, cable1),
		  filefullpath_(io_environment_.output_folder_ + "/cable" + cableinf + ".dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "length"
				 << " ";
		out_file << "  "
				 << "rate"
				 << " ";
        out_file << "  "
				 << "integ-rate"
				 << " ";
        out_file << "  "
        		 << "unitpow"
				 << " ";
        out_file << "  "
				 << "tension"
				 << " ";
        out_file << "  "
				 << "disswork"
				 << " "; 
		out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyCableData::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		
        const SimTK::State &state = integ_.getState();
        const SimTK::CablePath& path1 = cable_.getCablePath();
		
		out_file << "  " << path1.getCableLength(state) << "  ";
        out_file << "  " << path1.getCableLengthDot(state) << "  ";
        out_file << "  " << path1.getIntegratedCableLengthDot(state)<< "  ";
        out_file << "  " << path1.calcCablePower(state, 1)<< "  "; // unit power
        out_file << "  " << cable_.getTension(state)<< "  ";
        out_file << "  " << cable_.getDissipatedEnergy(state)<< "  ";
		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
}
