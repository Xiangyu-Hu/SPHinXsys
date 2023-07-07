/**
 * @file 	io_simbody_free.h
 * @brief 	Classes for simbody free output.
 * @author	Nicolò Salis
 */

#pragma once

#include "io_base.h"

namespace SPH
{
	/**
	 * @class WriteSimBodyFreeRotationMatrix
	 * @brief Write displacement and rotation of planar solid body.
	 */
	class WriteSimBodyFreeRotationMatrix : public WriteSimBodyStates<SimTK::MobilizedBody::Free>
	{
	protected:
		std::string filefullpath_;

	public:
		WriteSimBodyFreeRotationMatrix(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Free &freebody);
		virtual ~WriteSimBodyFreeRotationMatrix(){};
		virtual void writeToFile(size_t iteration_step = 0) override;
	};
	//=============================================================================================//
	WriteSimBodyFreeRotationMatrix::
		WriteSimBodyFreeRotationMatrix(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Free &freebody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Free>(io_environment, integ, freebody),
		  filefullpath_(io_environment_.output_folder_ + "/RotationMatrix.dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "rotation [1,1]"
				 << " ";
		out_file << "  "
				 << "rotation [1,2]"
				 << " ";
		out_file << "  "
				 << "rotation [1,3]"
				 << " ";
		out_file << "  "
				 << "rotation [2,1]"
				 << " ";
		out_file << "  "
				 << "rotation [2,2]"
				 << " ";
		out_file << "  "
				 << "rotation [2,3]"
				 << " ";
		out_file << "  "
				 << "rotation [3,1]"
				 << " ";
		out_file << "  "
				 << "rotation [3,2]"
				 << " ";
		out_file << "  "
				 << "rotation [3,3]"
				 << " ";
				 out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyFreeRotationMatrix::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		const SimTK::State &state = integ_.getState();

		out_file << "  " << mobody_.getBodyRotation(state)[0][0] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[0][1] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[0][2] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[1][0] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[1][1] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[1][2] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[2][0] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[2][1] << "  ";
		out_file << "  " << mobody_.getBodyRotation(state)[2][2] << "  ";

		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
/**
	 * @class WriteSimBodyVelocity
	 * @brief Write displacement and rotation of planar solid body.
	 */
	class WriteSimBodyVelocity : public WriteSimBodyStates<SimTK::MobilizedBody::Free>
	{
	protected:
		std::string filefullpath_;

	public:
		WriteSimBodyVelocity(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Free &freebody);
		virtual ~WriteSimBodyVelocity(){};
		virtual void writeToFile(size_t iteration_step = 0) override;
	};
	//=============================================================================================//
	WriteSimBodyVelocity::
		WriteSimBodyVelocity(IOEnvironment &io_environment, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Free &freebody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Free>(io_environment, integ, freebody),
		  filefullpath_(io_environment_.output_folder_ + "/BodyVelocity.dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "velocity [0]"
				 << " ";
		out_file << "  "
				 << "velocity [1]"
				 << " ";
		out_file << "  "
				 << "velocity [2]"
				 << " ";
		out_file << "  ";
				 out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyVelocity::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		const SimTK::State &state = integ_.getState();

		out_file << "  " << mobody_.getBodyOriginVelocity(state)[0] << "  ";
		out_file << "  " << mobody_.getBodyOriginVelocity(state)[1] << "  ";
		out_file << "  " << mobody_.getBodyOriginVelocity(state)[2] << "  ";

		out_file << "\n";
		out_file.close();
	};
	//=================================================================================================//
}
