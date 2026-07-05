#include "io_base.hpp"

#include "base_dynamics.h"
#include "io_environment.h"
#include "sph_system.h"

namespace SPH
{
//=============================================================================================//
BaseIO::BaseIO(SPHSystem &sph_system)
    : sph_system_(sph_system), io_environment_(IO::getEnvironment()),
      sv_physical_time_(&sph_system.svPhysicalTime()) {}
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
BodyStatesRecording::~BodyStatesRecording() = default;
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
RestartIO::RestartIO(SPHSystem &sph_system, bool summary_enabled)
    : BaseIO(sph_system), summary_enabled_(summary_enabled),
      real_bodies_(sph_system.getRealBodies()),
      overall_file_path_(io_environment_.RestartFolder() + "/Restart_")
{
    if (sph_system_.RestartStep() == 0)
    {
        io_environment_.resetForRestart();
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
        std::string body_name = real_bodies_[i]->Name();

        // Add a body element
        restart_xml.addNewElement(restart_xml.first_element_, "body");

        // Get the last added body element
        tinyxml2::XMLElement *body_element = restart_xml.first_element_->LastChildElement("body");

        // Set body name attribute
        restart_xml.setAttributeToElement(body_element, "name", body_name);

        // Write particles to this body element
        base_particles.writeParticlesToXml(restart_xml, body_element);
    }

    // Write the consolidated XML file
    restart_xml.writeToXmlFile(overall_filefullpath);

    if (summary_enabled_)
    {
        reportRestartSummary(iteration_step);
    }
}
//=============================================================================================//
void RestartIO::reportRestartSummary(size_t restart_step)
{
    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();
        std::string body_name = real_bodies_[i]->Name();

        std::cout << "Restart Information Summary:\n";
        std::cout << "---------------------------------------------\n";
        std::cout << "Total real particles of body " << body_name
                  << " written to restart: " << base_particles.TotalRealParticles() << "\n";
        std::cout << "---------------------------------------------\n";
    }
}
//=============================================================================================//
Real RestartIO::readRestartTime(size_t restart_step)
{
    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".xml";
    XmlParser restart_xml("xml_restart");
    restart_xml.loadXmlFile(overall_filefullpath);

    Real restart_time;
    restart_xml.queryAttributeValue(restart_xml.first_element_, "restart_time", restart_time);
    sv_physical_time_->setValue(restart_time);
    return restart_time;
}
//=============================================================================================//
void RestartIO::readFromFile(size_t restart_step)
{
    std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;

    std::string overall_filefullpath = overall_file_path_ + padValueWithZeros(restart_step) + ".xml";
    XmlParser restart_xml("xml_restart");
    restart_xml.loadXmlFile(overall_filefullpath);

    // Iterate through all body elements in the XML
    for (size_t i = 0; i < real_bodies_.size(); ++i)
    {
        std::string body_name = real_bodies_[i]->Name();
        BaseParticles &base_particles = real_bodies_[i]->getBaseParticles();

        // Find the body element by iterating through child elements
        tinyxml2::XMLElement *body_element = restart_xml.first_element_->FirstChildElement("body");
        bool found = false;

        while (body_element != nullptr)
        {
            const char *name_attr = body_element->Attribute("name");

            if (name_attr != nullptr && std::string(name_attr) == body_name)
            {
                found = true;
                base_particles.readParticlesFromXml(restart_xml, body_element);
                std::cout << "\n Total real particles of body " << body_name
                          << " read from restart: " << base_particles.TotalRealParticles() << "\n";
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
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBodyVector bodies)
    : BaseIO(bodies[0]->getSPHSystem()), bodies_(bodies),
      overall_file_path_(io_environment_.ReloadFolder() + "/Reload.xml")
{
    for (size_t i = 0; i < bodies_.size(); ++i)
        body_names_.push_back(bodies_[i]->Name());
}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHSystem &sph_system)
    : ReloadParticleIO(sph_system.getRealBodies()) {}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBody &sph_body, const std::string &given_body_name)
    : BaseIO(sph_body.getSPHSystem()), bodies_({&sph_body}),
      body_names_({given_body_name}),
      overall_file_path_(io_environment_.ReloadFolder() + "/Reload.xml") {}
//=============================================================================================//
ReloadParticleIO::ReloadParticleIO(SPHBody &sph_body)
    : ReloadParticleIO(sph_body, sph_body.Name()) {}
//=============================================================================================//
void ReloadParticleIO::writeToFile(size_t iteration_step)
{
    if (fs::exists(overall_file_path_))
    {
        fs::remove(overall_file_path_);
    }

    XmlParser reload_xml("xml_particle_reload", "reload_data");

    for (size_t i = 0; i < bodies_.size(); ++i)
    {
        BaseParticles &base_particles = bodies_[i]->getBaseParticles();
        std::string body_name = body_names_[i];

        reload_xml.addNewElement(reload_xml.first_element_, "body");
        tinyxml2::XMLElement *body_element = reload_xml.first_element_->LastChildElement("body");
        reload_xml.setAttributeToElement(body_element, "name", body_name);

        base_particles.writeParticlesToXml(reload_xml, body_element);
    }

    reload_xml.writeToXmlFile(overall_file_path_);
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
