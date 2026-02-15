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
    return result != bodies.end();
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
    writeWithFileName("ite_" + padValueWithZeros(iteration_step));
};
//=============================================================================================//
RestartIO::RestartIO(SPHSystem &sph_system)
    : BaseIO(sph_system), real_bodies_(sph_system.getRealBodies()),
      overall_file_path_(io_environment_.RestartFolder() + "/Restart_")
{
    if (sph_system_.RestartStep() == 0)
    {
        io_environment_.resetForRestart();
    }

    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        file_names_.push_back(io_environment_.RestartFolder() + "/" + real_bodies_[i]->getName() + "_rst_");
    }
}
//=============================================================================================//
void RestartIO::writeToFile(size_t iteration_step)
{
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(iteration_step) + ".xml";
    if (fs::exists(overall_filefullpath))
    {
        fs::remove(overall_filefullpath);
    }

    // Create a new XML document for restart
    XmlParser restart_xml("xml_restart", "restart_data");
    
    // Add restart time as an attribute to the root element
    restart_xml.setAttributeToElement(restart_xml.first_element_, "restart_time", sv_physical_time_->getValue());

    // Write all bodies to the single XML file
    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
        std::string body_name = real_bodies_[i]->getName();
        
        std::cout << "\n Total real particles of body " << body_name
                  << " write to restart is " << base_particles.TotalRealParticles() << "\n";
        
        // Add a body element
        restart_xml.addNewElement(restart_xml.first_element_, "body");
        
        // Get the last added body element
        tinyxml2::XMLElement *body_element = restart_xml.first_element_->LastChildElement("body");
        
        // Set body name attribute
        restart_xml.setAttributeToElement(body_element, "name", body_name);
        
        // Write particles to this body element
        base_particles.writeParticlesToXmlForRestart(restart_xml, body_element);
    }
    
    // Write the consolidated XML file
    restart_xml.writeToXmlFile(overall_filefullpath);
}
//=============================================================================================//
Real RestartIO::readRestartTime(size_t restart_step)
{
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".xml";
    
    // Check for new format first
    if (fs::exists(overall_filefullpath))
    {
        XmlParser restart_xml("xml_restart");
        restart_xml.loadXmlFile(overall_filefullpath);
        
        Real restart_time;
        restart_xml.queryAttributeValue(restart_xml.first_element_, "restart_time", restart_time);
        return restart_time;
    }
    
    // Fallback to old format for backward compatibility
    std::string old_filefullpath = io_environment_.RestartFolder() + "/Restart_time_" + padValueWithZeros(restart_step) + ".dat";
    if (!fs::exists(old_filefullpath))
    {
        std::cout << "\n Error: the input file:" << overall_filefullpath << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    Real restart_time;
    std::ifstream in_file(old_filefullpath.c_str());
    in_file >> restart_time;
    in_file.close();

    return restart_time;
}
//=============================================================================================//
void RestartIO::readFromFile(size_t restart_step)
{
    std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
    
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".xml";
    
    // Check for new consolidated format first
    if (fs::exists(overall_filefullpath))
    {
        XmlParser restart_xml("xml_restart");
        restart_xml.loadXmlFile(overall_filefullpath);
        
        // Iterate through all body elements in the XML
        for (size_t i = 0; i < real_bodies_.size(); ++i)
        {
            std::string body_name = real_bodies_[i]->getName();
            BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
            
            // Find the body element by iterating through child elements
            tinyxml2::XMLElement *body_element = restart_xml.first_element_->FirstChildElement("body");
            bool found = false;
            
            while (body_element != nullptr)
            {
                const char *name_attr = nullptr;
                body_element->QueryAttribute("name", &name_attr);
                
                if (name_attr != nullptr && std::string(name_attr) == body_name)
                {
                    found = true;
                    base_particles.readParticlesFromXmlForRestart(restart_xml, body_element);
                    std::cout << "\n Total real particles of body " << body_name
                              << " from restart is " << base_particles.TotalRealParticles() << "\n";
                    break;
                }
                
                body_element = body_element->NextSiblingElement("body");
            }
            
            if (!found)
            {
                std::cout << "\n Error: body " << body_name << " not found in restart file: " 
                          << overall_filefullpath << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
        }
    }
    else
    {
        // Fallback to old format for backward compatibility
        for (size_t i = 0; i < real_bodies_.size(); ++i)
        {
            std::string filefullpath = file_names_[i] + padValueWithZeros(restart_step) + ".xml";

            if (!fs::exists(filefullpath))
            {
                std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
            BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
            base_particles.readParticlesFromXmlForRestart(filefullpath);
        }
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
