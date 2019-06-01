/**
 * @file 	output.h
 * @brief 	Base calss for all output functions.
 * @author	Chi Zhang and Xiangyu Hu
 * @version	0.1
 */

#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_physical_dynamics.h"
 
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

#include <fstream>
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH {

	class SPHBody;
	class SPHSystem;
	class WeaklyCompressibleFluidBody;
	class ObserverBody;
	class VelocityBoundCheck;

	/**
	 * @class Output
	 * @brief The base class which defines folders for output and restart.
	 */
	class Output
	{
	protected:
		SPHSystem &sph_system_;

	public:
		Output(SPHSystem &sph_system);
		virtual ~Output() {};

		std::string output_folder_;
		std::string restart_folder_;

		void WriteCaseSetup(Real start_time, Real D_Time, Real end_Time);
	};

	/**
	 * @class WriteBodyStates
	 * @brief base class for write body states.
	 */
	class WriteBodyStates
	{
	protected:
		std::string output_folder_;
		std::string restart_folder_;
		Output &output_;
		SPHBodyVector bodies_;
		SPHBody *body_;

	public:
		/** write for one body */
		WriteBodyStates(Output &output, SPHBody *body);
		/** write for several bodies */
		WriteBodyStates(Output &output, SPHBodyVector bodies);
		virtual ~WriteBodyStates() {};

		std::string filefullpath_;
		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class WriteSimBodyStates
	 * @brief base class for write SimBody states.
	 */
	template<class MobilizedBodyType>
	class WriteSimBodyStates
	{
	protected:
		std::string output_folder_;
		std::string restart_folder_;
		Output& output_;
		StdVec<MobilizedBodyType*> mobodies_;

	public:
		/** write for several bodies */
		WriteSimBodyStates(Output& output, StdVec<MobilizedBodyType*> mobodies) 
			: output_(output), mobodies_(mobodies)
		{
			output_folder_ = output_.output_folder_;
			restart_folder_ = output_.restart_folder_;
		};
		virtual ~WriteSimBodyStates() {};

		std::string filefullpath_;
		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class WriteBodyStatesToVtu
	 * @brief  Write files for bodies
	 * the output file is VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class WriteBodyStatesToVtu : public WriteBodyStates
	{
	public:
		WriteBodyStatesToVtu(Output &output, SPHBodyVector bodies);
		virtual ~WriteBodyStatesToVtu() {};

		virtual void WriteToFile(Real time) override;
	};
	
	/**
	 * @class WriteBodyStatesToPlt
	 * @brief  Write files for bodies
	 * the output file is dat format can visualized by TecPlot
	 */
	class WriteBodyStatesToPlt : public WriteBodyStates
	{
	public:
		WriteBodyStatesToPlt(Output &output, SPHBodyVector bodies);
		virtual ~WriteBodyStatesToPlt() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteToPltIfVelocityOutOfBound
	 * @brief  output body sates if paritcle velocity
	 * out of a bound
	 */
	class WriteToPltIfVelocityOutOfBound
		: public WriteBodyStatesToPlt
	{
	protected:
		StdVec<VelocityBoundCheck *> check_bodies_;
	public:
		WriteToPltIfVelocityOutOfBound(Output &output,
			SPHBodyVector bodies, Real velocity_bound);
		virtual ~WriteToPltIfVelocityOutOfBound() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteRelaxBodyMeshToPlt
	 * @brief  write the background mesh data for relax body
	 */
	class WriteRelaxBodyMeshToPlt : public WriteBodyStates
	{
	protected:
		RelaxBody * relax_body_;

	public:
		WriteRelaxBodyMeshToPlt(Output &output, RelaxBody *relax_body);
		virtual ~WriteRelaxBodyMeshToPlt() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteObservedFluidPressure
	 * @brief write files for observed fluid pressure 
	 */
	class WriteObservedFluidPressure
		: public WriteBodyStates,
		public observer_dynamics::ObserveFluidPressure
	{
	protected:
		ObserverBody *observer_;
	public:
		WriteObservedFluidPressure(Output &output,
			ObserverBody* observer, WeaklyCompressibleFluidBody *interacting_body);
		virtual ~WriteObservedFluidPressure() {};
		virtual void WriteToFile(Real time) override;

	};

	/**
	* @class WriteObservedFluidPressure
	* @brief write files for observed fluid pressure
	*/
	class WriteObservedFluidVelocity
		: public WriteBodyStates,
		public observer_dynamics::ObserveFluidVelocity
	{
	protected:
		size_t dimension_;
		ObserverBody * observer_;
	public:
		WriteObservedFluidVelocity(Output &output,
			ObserverBody* observer, WeaklyCompressibleFluidBody *interacting_body);
		virtual ~WriteObservedFluidVelocity() {};
		virtual void WriteToFile(Real time) override;

	};

	/**
	 * @class WriteObservedElasticDisplacement
	 * @brief write files for observed elastic displacement
	 */
	class WriteObservedElasticDisplacement
		: public WriteBodyStates, 
		public observer_dynamics::ObserveElasticDisplacement
	{
	protected:
		size_t dimension_;
		ObserverBody *observer_;
	public:
		WriteObservedElasticDisplacement(Output &output,
			ObserverBody* observer, ElasticBody *interacting_body);
		virtual ~WriteObservedElasticDisplacement() {};
		virtual void WriteToFile(Real time) override;

	};

	/**
	 * @class WriteWaterMechanicalEnergy
	 * @brief write files for the total mechanical energy of a weakly compressible fluid body
	 */
	class WriteWaterMechanicalEnergy 
		: public WriteBodyStates, public fluid_dynamics::TotalMechanicalEnergy
	{
	public:
		WriteWaterMechanicalEnergy(Output &output, 
			WeaklyCompressibleFluidBody* water_block, ExternalForce* external_force);
		virtual ~WriteWaterMechanicalEnergy() {};
		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteTotalViscousForceOnSolid
	 * @brief write total viscous force acting a solid body
	 */
	class WriteTotalViscousForceOnSolid
		: public WriteBodyStates, public solid_dynamics::TotalViscousForceOnSolid
	{
	protected:
		size_t dimension_;

	public:
		WriteTotalViscousForceOnSolid(Output &output, SolidBody *solid_body);
		virtual ~WriteTotalViscousForceOnSolid() {};
		virtual void WriteToFile(Real time) override;
	};

	/** 
	 * @class WriteTotalForceOnSolid
	 * @brief Write total force acting a solid body.
	 */
	class WriteTotalForceOnSolid
		: public WriteBodyStates, public solid_dynamics::TotalForceOnSolid
	{
	protected:
		size_t dimension_;

	public:
		WriteTotalForceOnSolid(Output &output, SolidBody *solid_body);
		virtual ~WriteTotalForceOnSolid() {};
		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteWaterFront
	 * @brief write files for water front in free surface flow
	 */
	class WriteWaterFront
		: public WriteBodyStates, public fluid_dynamics::DambreakWaterFront
	{
	public:
		WriteWaterFront(Output &output, WeaklyCompressibleFluidBody* water_block);
		virtual ~WriteWaterFront() {};
		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteBodyStatesToXml
	 * @brief write files for bodies 
	 * the output file is XML format can be read by XML parser
	 */
	class WriteBodyStatesToXml : public WriteBodyStates
	{
	public:
		WriteBodyStatesToXml(Output &output, SPHBodyVector bodies);
		virtual ~WriteBodyStatesToXml() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteRestartFileToXml
	 * @brief Write the restart file in XML format.
	 */
	class WriteRestartFileToXml : public WriteBodyStates
	{
	public:
		WriteRestartFileToXml(Output &output, SPHBodyVector bodies);
		virtual ~WriteRestartFileToXml() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteTotalForceOnSolid
	* @brief Write total force acting a solid body.
	*/
	class WriteSimBodyPinAngle : public WriteSimBodyStates<SimTK::MobilizedBody::Pin>
	{
	protected:
		SimTK::RungeKuttaMersonIntegrator &integ_;
	public:
		WriteSimBodyPinAngle(Output& output, StdVec<SimTK::MobilizedBody::Pin*> mobodies, SimTK::RungeKuttaMersonIntegrator &integ);
		virtual ~WriteSimBodyPinAngle() {};
		virtual void WriteToFile(Real time) override;
	};
}
