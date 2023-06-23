
/**
 * @file 	io_base.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "io_base.h"

#include "sph_system.h"

namespace SPH
{
//=============================================================================================//
IOEnvironment::IOEnvironment(SPHSystem &sph_system, bool delete_output)
    : sph_system_(sph_system),
      input_folder_("./input"), output_folder_("./output"),
      restart_folder_("./restart"), reload_folder_("./reload")
{
    if (!fs::exists(input_folder_))
    {
        fs::create_directory(input_folder_);
    }

    if (!fs::exists(output_folder_))
    {
        fs::create_directory(output_folder_);
    }

    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }

    if (!fs::exists(reload_folder_))
    {
        fs::create_directory(reload_folder_);
    }

    if (sph_system.RestartStep() == 0)
    {
        fs::remove_all(restart_folder_);
        fs::create_directory(restart_folder_);
        if (delete_output == true)
        {
            fs::remove_all(output_folder_);
            fs::create_directory(output_folder_);
        }
    }

    sph_system.io_environment_ = this;
}
//=============================================================================================//
ParameterizationIO &IOEnvironment::defineParameterizationIO()
{
    return parameterization_io_ptr_keeper_.createRef<ParameterizationIO>(input_folder_);
}
//=============================================================================================//
std::string BaseIO::convertPhysicalTimeToString(Real convertPhysicalTimeToStream)
{
    int i_time = int(GlobalStaticVariables::physical_time_ * 1.0e6);
    return padValueWithZeros(i_time);
}
//=============================================================================================//
void BodyStatesRecording::writeToFile()
{
    writeWithFileName(convertPhysicalTimeToString(GlobalStaticVariables::physical_time_));
}
//=============================================================================================//
void BodyStatesRecording::writeToFile(size_t iteration_step)
{
    writeWithFileName(padValueWithZeros(iteration_step));
};
//=============================================================================================//
RestartIO::RestartIO(IOEnvironment &io_environment, SPHBodyVector bodies)
    : BaseIO(io_environment), bodies_(bodies),
      overall_file_path_(io_environment.restart_folder_ + "/Restart_time_")
{
    std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_names_),
                   [&](SPHBody *body) -> std::string
                   { return io_environment.restart_folder_ + "/SPHBody_" + body->getName() + "_rst_"; });
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
    out_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "   \n";
    out_file.close();

    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i] + padValueWithZeros(iteration_step) + ".xml";

        if (fs::exists(filefullpath))
        {
            fs::remove(filefullpath);
        }
        bodies_[i]->writeParticlesToXmlForRestart(filefullpath);
    }
}
//=============================================================================================//
Real RestartIO::readRestartTime(size_t restart_step)
{
    std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
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
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i] + padValueWithZeros(restart_step) + ".xml";

        if (!fs::exists(filefullpath))
        {
            std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        bodies_[i]->readParticlesFromXmlForRestart(filefullpath);
    }
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(IOEnvironment &io_environment, SPHBodyVector bodies)
    : BaseIO(io_environment), bodies_(bodies)
{
    std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_names_),
                   [&](SPHBody *body) -> std::string
                   { return io_environment.reload_folder_ + "/" + body->getName() + "_rld.xml"; });
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(IOEnvironment &io_environment, SPHBody &sph_body,
                                   const std::string &given_body_name)
    : BaseIO(io_environment), bodies_({&sph_body})
{
    file_names_.push_back(io_environment.reload_folder_ + "/" + given_body_name + "_rld.xml");
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(IOEnvironment &io_environment, SPHBody &sph_body)
    : ReloadParticleIO(io_environment, sph_body, sph_body.getName()) {}
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
        bodies_[i]->writeToXmlForReloadParticle(filefullpath);
    }
}
//=============================================================================================//
void ReloadParticleIO::readFromFile(size_t restart_step)
{
    std::cout << "\n Reloading particles from files." << std::endl;
    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        std::string filefullpath = file_names_[i];

        if (!fs::exists(filefullpath))
        {
            std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        bodies_[i]->readFromXmlForReloadParticle(filefullpath);
    }
}
//=================================================================================================//
} // namespace SPH
