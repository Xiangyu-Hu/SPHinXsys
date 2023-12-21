#include "io_simbody.h"

namespace SPH
{
//=============================================================================================//
WriteSimBodyPinData::
    WriteSimBodyPinData(SPHSystem &sph_system,
                        SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody)
    : WriteSimBodyStates<SimTK::MobilizedBody::Pin>(sph_system, integ, pinbody),
      filefullpath_(io_environment_.output_folder_ + "/mb_pinbody_data.dat")
{
    std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

    out_file << "\"time\""
             << "   ";
    out_file << "  "
             << "angles"
             << " ";
    out_file << "  "
             << "angle_rates"
             << " ";
    out_file << "\n";

    out_file.close();
};
//=============================================================================================//
void WriteSimBodyPinData::writeToFile(size_t iteration_step)
{
    std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
    out_file << GlobalStaticVariables::physical_time_ << "   ";
    const SimTK::State &state = integ_.getState();

    out_file << "  " << mobody_.getAngle(state) << "  " << mobody_.getRate(state) << "  ";

    out_file << "\n";
    out_file.close();
};
//=============================================================================================//
WriteSimBodyCableData::
    WriteSimBodyCableData(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                          SimTK::CableSpring &cable1, std::string cable_inf)
    : WriteSimBodyStates<SimTK::CableSpring>(sph_system, integ, cable1),
      filefullpath_(io_environment_.output_folder_ + "/cable" + cable_inf + ".dat")
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
    const SimTK::CablePath &path1 = mobody_.getCablePath();

    out_file << "  " << path1.getCableLength(state) << "  ";
    out_file << "  " << path1.getCableLengthDot(state) << "  ";
    out_file << "  " << path1.getIntegratedCableLengthDot(state) << "  ";
    out_file << "  " << path1.calcCablePower(state, 1) << "  "; // unit power
    out_file << "  " << mobody_.getTension(state) << "  ";
    out_file << "  " << mobody_.getDissipatedEnergy(state) << "  ";
    out_file << "\n";
    out_file.close();
};
//=============================================================================================//
WriteSimBodyPlanarData::
    WriteSimBodyPlanarData(SPHSystem &sph_system, SimTK::RungeKuttaMersonIntegrator &integ,
                           SimTK::MobilizedBody::Planar &planar_body)
    : WriteSimBodyStates<SimTK::MobilizedBody::Planar>(sph_system, integ, planar_body),
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
} // namespace SPH
