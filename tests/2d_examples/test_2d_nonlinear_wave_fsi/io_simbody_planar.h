/**
 * @file 	io_simbody_planar.h
 * @brief 	Classes for simbody planar output (displacement and rotation).
 * @author	Nicolò Salis
 */

#pragma once

#include "io_base.h"

namespace SPH
{
	/**
	 * @class WriteSimBodyPlanarData
	 * @brief Write displacement and rotation of planar solid body.
	 */
	class WriteSimBodyPlanarData : public WriteSimBodyStates<SimTK::MobilizedBody::Planar>
	{
	protected:
		std::string filefullpath_;

	public:
		WriteSimBodyPlanarData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Planar &planarbody);
		virtual ~WriteSimBodyPlanarData(){};
		virtual void writeToFile(size_t iteration_step = 0) override;
	};
	//=============================================================================================//
	WriteSimBodyPlanarData::
		WriteSimBodyPlanarData(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Planar &planarbody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Planar>(io_environment, integ, planarbody),
		  filefullpath_(io_environment_.output_folder_ + "/mb_planar_data.dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "translation x"
				 << " ";
		out_file << "  "
				 << "translation y"
				 << " ";
		out_file << "  "
				 << "angle"
				 << " ";
		out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyPlanarData::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		const SimTK::State &state = integ_.getState();

		out_file << "  " << mobody_.getTranslation(state)[0] << "  "; 
		out_file << "  " << mobody_.getTranslation(state)[1] << "  ";
		out_file << "  " << mobody_.getAngle(state) << "  ";
		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
}
