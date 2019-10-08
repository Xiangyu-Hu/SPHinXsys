/**
 * @file 	in_output.h
 * @brief 	Classes for input and output functions.
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
	/**
	 * @class In_Output
	 * @brief The base class which defines folders for output, 
	 * restart and particle reload folders.
	 */
	class In_Output
	{
	public:
		In_Output(SPHSystem &sph_system);
		virtual ~In_Output() {};

		SPHSystem &sph_system_;
		std::string output_folder_;
		std::string restart_folder_;
		std::string particle_reload_folder_;
		std::string restart_step_;
	};

	/**
	 * @class BodyStatesIO
	 * @brief base class for write and read body states.
	 */
	class BodyStatesIO
	{
	protected:
		In_Output &in_output_;
		SPHBody *body_;
		SPHBodyVector bodies_;
	public:
		BodyStatesIO(In_Output &in_output, SPHBody *body)
			: in_output_(in_output), body_(body) {};
		BodyStatesIO(In_Output &in_output, SPHBodyVector bodies)
			: in_output_(in_output), bodies_(bodies) {};
		virtual ~BodyStatesIO() {};
	};

	/**
	 * @class WriteBodyStates
	 * @brief base class for write body states.
	 */
	class WriteBodyStates : public BodyStatesIO
	{
	public:
		WriteBodyStates(In_Output &in_output, SPHBody *body)
			: BodyStatesIO(in_output, body) {};
		WriteBodyStates(In_Output &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies) {};
		virtual ~WriteBodyStates() {};

		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class ReadBodyStates
	 * @brief base class for read body states.
	 */
	class ReadBodyStates : public BodyStatesIO
	{
	public:
		ReadBodyStates(In_Output &in_output, SPHBody *body)
			: BodyStatesIO(in_output, body) {};
		ReadBodyStates(In_Output &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies) {};
		virtual ~ReadBodyStates() {};

		virtual void ReadFromFile(size_t iteration_step) = 0;
	};
	/**
	 * @class SimBodyStatesIO
	 * @brief base class for write and read SimBody states.
	 */
	template<class MobilizedBodyType>
	class SimBodyStatesIO
	{
	protected:
		In_Output &in_output_;
		MobilizedBodyType* mobody_;
		StdVec<MobilizedBodyType*> mobodies_;
	public:
		SimBodyStatesIO(In_Output& in_output, MobilizedBodyType* mobody)
			: in_output_(in_output), mobody_(mobody) {};
		SimBodyStatesIO(In_Output& in_output, StdVec<MobilizedBodyType*> mobodies)
			: in_output_(in_output), mobodies_(mobodies) {};
		virtual ~SimBodyStatesIO() {};
	};

	/**
	 * @class WriteSimBodyStates
	 * @brief base class for write SimBody states.
	 */
	template<class MobilizedBodyType>
	class WriteSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		WriteSimBodyStates(In_Output& in_output, MobilizedBodyType* mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobody) {};
		WriteSimBodyStates(In_Output& in_output, StdVec<MobilizedBodyType*> mobodies)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobodies) {};
		virtual ~WriteSimBodyStates() {};

		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class ReadSimBodyStates
	 * @brief base class for read SimBody states.
	 */
	template<class MobilizedBodyType>
	class ReadSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		ReadSimBodyStates(In_Output& in_output, MobilizedBodyType* mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobody) {};
		ReadSimBodyStates(In_Output& in_output, StdVec<MobilizedBodyType*> mobodies)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobodies) {};
		virtual ~ReadSimBodyStates() {};

		virtual void ReadFromFile(size_t iteration_step) = 0;
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
		WriteBodyStatesToVtu(In_Output& in_output, SPHBodyVector bodies)
			: WriteBodyStates(in_output, bodies) {};
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
		WriteBodyStatesToPlt(In_Output& in_output, SPHBodyVector bodies)
			: WriteBodyStates(in_output, bodies) {};;
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
		WriteToPltIfVelocityOutOfBound(In_Output& in_output,
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
		RelaxBody *relax_body_;
		std::string filefullpath_;
	public:
		WriteRelaxBodyMeshToPlt(In_Output& in_output, RelaxBody *relax_body);
		virtual ~WriteRelaxBodyMeshToPlt() {};

		virtual void WriteToFile(Real time = 0.0) override;
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
		std::string filefullpath_;
	public:
		WriteObservedFluidPressure(In_Output& in_output,
			ObserverBody* observer, FluidBody *interacting_body);
		virtual ~WriteObservedFluidPressure() {};
		virtual void WriteToFile(Real time = 0.0) override;
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
		std::string filefullpath_;
	public:
		WriteObservedFluidVelocity(In_Output& in_output,
			ObserverBody* observer, FluidBody *interacting_body);
		virtual ~WriteObservedFluidVelocity() {};
		virtual void WriteToFile(Real time = 0.0) override;
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
		std::string filefullpath_;
	public:
		WriteObservedElasticDisplacement(In_Output& in_output,
			ObserverBody* observer, SolidBody *interacting_body);
		virtual ~WriteObservedElasticDisplacement() {};
		virtual void WriteToFile(Real time = 0.0) override;

	};

	/**
	 * @class WriteObservedVoltage
	 * @brief write the observed voltage of electrophysiology to files.
	 */
	class WriteObservedVoltage
		: public WriteBodyStates,
		public observer_dynamics::ObserveMuscleVoltage
	{
	protected:
		ObserverBody *observer_;
		std::string filefullpath_;
	public:
		/** Constructor and Destructor. */
		WriteObservedVoltage(In_Output& in_output,
			ObserverBody* observer, SolidBody *interacting_body);
		virtual ~WriteObservedVoltage() {};
		/**
		 * @brief Output data to files.
		 * @param[in] time Physical time.
		 */
		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteWaterMechanicalEnergy
	 * @brief write files for the total mechanical energy of a weakly compressible fluid body
	 */
	class WriteWaterMechanicalEnergy 
		: public WriteBodyStates, public fluid_dynamics::TotalMechanicalEnergy
	{
	protected:
		std::string filefullpath_;
	public:
		WriteWaterMechanicalEnergy(In_Output& in_output, FluidBody* water_block, ExternalForce* external_force);
		virtual ~WriteWaterMechanicalEnergy() {};
		virtual void WriteToFile(Real time = 0.0) override;
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
		std::string filefullpath_;
	public:
		WriteTotalViscousForceOnSolid(In_Output& in_output, SolidBody *solid_body);
		virtual ~WriteTotalViscousForceOnSolid() {};
		virtual void WriteToFile(Real time = 0.0) override;
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
		std::string filefullpath_;
	public:
		WriteTotalForceOnSolid(In_Output& in_output, SolidBody *solid_body);
		virtual ~WriteTotalForceOnSolid() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteWaterFront
	 * @brief write files for water front in free surface flow
	 */
	class WriteUpperBoundInXDirection
		: public WriteBodyStates, public UpperBoundInXDirection
	{
	protected:
		std::string filefullpath_;
	public:
		WriteUpperBoundInXDirection(In_Output& in_output, SPHBody* body);
		virtual ~WriteUpperBoundInXDirection() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class ReloadParticleIO
	 * @brief For write  and read particle reload.
	 */
	class ReloadParticleIO
	{
	protected:
		StdVec<std::string> file_paths_;

	public:
		ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~ReloadParticleIO() {};
	};

	/**
	  * @class WriteReloadParticle
	  * @brief Write the reload particles file in XML format.
	  */
	class WriteReloadParticle : public ReloadParticleIO, public WriteBodyStates
	{
	public:
		WriteReloadParticle(In_Output& in_output, SPHBodyVector bodies)
			: ReloadParticleIO(in_output, bodies), WriteBodyStates(in_output, bodies) {};
		virtual ~WriteReloadParticle() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	  * @class ReadReloadParticle
	  * @brief Write the reload particles file in XML format.
	  */
	class ReadReloadParticle : public ReloadParticleIO, public ReadBodyStates
	{
	public:
		ReadReloadParticle(In_Output& in_output, SPHBodyVector bodies, StdVec<std::string> reload_body_names);
		virtual ~ReadReloadParticle() {};

		virtual void ReadFromFile(size_t iteration_step = 0) override;
	};

	/**
	 * @class RestartIO
	 * @brief Write the restart file in XML format.
	 */
	class RestartIO
	{
	protected:
		std::string overall_file_path_;
		StdVec<std::string> file_paths_;

	public:
		RestartIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~RestartIO() {};
	};

	/**
	 * @class WriteRestart
	 * @brief Write the restart file in XML format.
	 */
	class WriteRestart : public RestartIO, public WriteBodyStates
	{
	public:
		WriteRestart(In_Output& in_output, SPHBodyVector bodies)
			: RestartIO(in_output, bodies), WriteBodyStates(in_output, bodies) {};
		virtual ~WriteRestart() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteRestart
	 * @brief Write the restart file in XML format.
	 */
	class ReadRestart : public RestartIO, public ReadBodyStates
	{
	protected:
		Real ReadRestartTime(size_t restart_step);
	public:
		ReadRestart(In_Output& in_output, SPHBodyVector bodies)
			: RestartIO(in_output, bodies), ReadBodyStates(in_output, bodies) {};
		virtual ~ReadRestart() {};
		virtual Real ReadRestartFiles(size_t restart_step) {
			ReadFromFile(restart_step);
			return ReadRestartTime(restart_step);
		};
		virtual void ReadFromFile(size_t iteration_step = 0) override;
	};

	/**
	 * @class WriteTotalForceOnSolid
	* @brief Write total force acting a solid body.
	*/
	class WriteSimBodyPinAngleAndAngleRate : public WriteSimBodyStates<SimTK::MobilizedBody::Pin>
	{
	protected:
		SimTK::RungeKuttaMersonIntegrator &integ_;
		std::string filefullpath_;
		SimTK::Visualizer *visulizer;
	public:
		WriteSimBodyPinAngleAndAngleRate(In_Output& in_output, StdVec<SimTK::MobilizedBody::Pin*> mobodies, SimTK::RungeKuttaMersonIntegrator &integ);
		virtual ~WriteSimBodyPinAngleAndAngleRate() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};
}
