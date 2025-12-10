#include "io_base.h"

#include "sph_system.hpp"

namespace SPH
{
//=============================================================================================//
BaseIO::BaseIO(SPHSystem &sph_system)
    : sph_system_(sph_system), io_environment_(sph_system.getIOEnvironment()),
      sv_physical_time_(sph_system_.getSystemVariableByName<Real>("PhysicalTime")) {}
//=============================================================================================//
std::string BaseIO::convertPhysicalTimeToString(Real convertPhysicalTimeToStream)
{
    size_t i_time = size_t(sv_physical_time_->getValue() * 1.0e6);
    return padValueWithZeros(i_time);
}
//=============================================================================================/
bool BaseIO::isBodyIncluded(const SPHBodyVector &bodies, SPHBody *sph_body)
{
    auto result = std::find_if(bodies.begin(), bodies.end(),
                               [&](auto &body) -> bool
                               { return body == sph_body; });
    return result != bodies.end() ? true : false;
}
//=============================================================================================//
BodyStatesRecording::BodyStatesRecording(SPHSystem &sph_system)
    : BaseIO(sph_system), bodies_(sph_system.getRealBodies()),
      state_recording_(sph_system_.StateRecording()) {}
//=============================================================================================//
BodyStatesRecording::BodyStatesRecording(SPHBody &body)
    : BaseIO(body.getSPHSystem()), bodies_({&body}),
      state_recording_(sph_system_.StateRecording()) {}
//=============================================================================================//
void BodyStatesRecording::writeToFile()
{
    for (auto &derived_variable : derived_variables_)
    {
        derived_variable->exec();
    }
    writeWithFileName(convertPhysicalTimeToString(sv_physical_time_->getValue()));
}
//=============================================================================================//
void BodyStatesRecording::writeToFile(size_t iteration_step)
{
    for (auto &derived_variable : derived_variables_)
    {
        derived_variable->exec();
    }
    writeWithFileName(padValueWithZeros(iteration_step));
};
//=============================================================================================//
RestartIO::RestartIO(SPHSystem &sph_system)
    : BaseIO(sph_system), bodies_(sph_system.getSPHBodies()),
      overall_file_path_(io_environment_.RestartFolder() + "/Restart_time_")
{
    if (sph_system_.RestartStep() == 0)
    {
        io_environment_.resetForRestart();
    }

    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        file_names_.push_back(io_environment_.RestartFolder() + "/" + bodies_[i]->getName() + "_rst_");
    }
}
//=============================================================================================//
void RestartIO::writeToFile(size_t iteration_step)
{
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(iteration_step) + ".dat";
    if (fs::exists(overall_filefullpath))
    {
        fs::remove(overall_filefullpath);
    }
    std::ofstream out_file(overall_filefullpath.c_str(), std::ios::app);
    out_file << std::fixed << std::setprecision(9) << sv_physical_time_->getValue() << "   \n";
    out_file.close();

    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i] + padValueWithZeros(iteration_step) + ".xml";

        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }
        BaseParticles &base_particles = bodies_[i]->getBaseParticles();
        base_particles.writeParticlesToXmlForRestart(filefullpath);
    }
}
//=============================================================================================//
Real RestartIO::readRestartTime(size_t restart_step)
{
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".dat";
    if (!fs::exists(overall_filefullpath))
    {
        std::cout << "\n Error: the input file:" << overall_filefullpath << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    Real restart_time;
    std::ifstream in_file(overall_filefullpath.c_str());
    in_file >> restart_time;
    in_file.close();

    return restart_time;
}
//=============================================================================================//
void RestartIO::readFromFile(size_t restart_step)
{
    std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i] + padValueWithZeros(restart_step) + ".xml";

        if (!fs::exists(filefullpath))
        {
            std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        BaseParticles &base_particles = bodies_[i]->getBaseParticles();
        base_particles.readParticlesFromXmlForRestart(filefullpath);
    }
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBodyVector bodies)
    : BaseIO(bodies[0]->getSPHSystem()), bodies_(bodies)
{
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        file_names_.push_back(io_environment_.ReloadFolder() + "/" + bodies_[i]->getName() + "_rld.xml");
    }
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHSystem &sph_system)
    : ReloadParticleIO(sph_system.getRealBodies()) {}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBody &sph_body, const std::string &given_body_name)
    : BaseIO(sph_body.getSPHSystem()), bodies_({&sph_body})
{
    file_names_.push_back(io_environment_.ReloadFolder() + "/" + given_body_name + "_rld.xml");
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBody &sph_body)
    : ReloadParticleIO(sph_body, sph_body.getName()) {}
//=============================================================================================//
void ReloadParticleIO::writeToFile(size_t iteration_step)
{
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i];

        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }
        BaseParticles &base_particles = bodies_[i]->getBaseParticles();
        base_particles.writeParticlesToXmlForReload(filefullpath);
    }
}
//=============================================================================================//
ParticleGenerationRecording::ParticleGenerationRecording(SPHBody &sph_body)
    : BaseIO(sph_body.getSPHSystem()), sph_body_(sph_body),
      state_recording_(sph_system_.StateRecording()) {}
//=============================================================================================//
void ParticleGenerationRecording::writeToFile(size_t iteration_step)
{
    writeWithFileName(padValueWithZeros(iteration_step));
}
//=================================================================================================//
} // namespace SPH
