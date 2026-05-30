#include "base_particles.hpp"

#include "base_body.h"
#include "base_material.h"

namespace SPH
{
//=================================================================================================//
BaseParticles::BaseParticles(SPHBody &sph_body)
    : sv_total_real_particles_(nullptr),
      particles_bound_(0), original_id_(nullptr), sorted_id_(nullptr),
      dv_pos_(nullptr), Vol_(nullptr), rho_(nullptr), mass_(nullptr),
      sph_body_(sph_body), body_name_(sph_body.Name()),
      reload_xml_parser_(*xml_parser_ptrs_.createPtr<XmlParser>("xml_particle_reload", "particles")),
      total_body_parts_(0)
{
    sph_body.assignBaseParticles(this);
    sv_total_real_particles_ = registerSingleVariable<UnsignedInt>("TotalRealParticles");
}
//=================================================================================================//
BaseParticles::~BaseParticles() = default;
//=================================================================================================//
SPHAdaptation &BaseParticles::getSPHAdaptation()
{
    return sph_body_.getSPHAdaptation();
}
//=================================================================================================//
std::string BaseParticles::getBodyName()
{
    return sph_body_.Name();
};
//=================================================================================================//
void BaseParticles::initializeBasicDiscreteVariables()
{
    addEvolvingVariable<Vecd>("Position");
    addEvolvingVariable<Real>("VolumetricMeasure");
    //----------------------------------------------------------------------
    //		register non-geometric variables
    //----------------------------------------------------------------------
    rho_ = registerStateVariableData<Real>("Density", sph_body_.getMatterMaterial().ReferenceDensity());
    mass_ = registerStateVariableData<Real>("Mass",
                                            [&](UnsignedInt i) -> Real
                                            { return rho_[i] * ParticleVolume(i); });
    //----------------------------------------------------------------------
    //		unregistered variables and data
    //----------------------------------------------------------------------
    original_id_ = registerDiscreteVariableData<UnsignedInt>("OriginalID", particles_bound_, AssignIndex());
    addEvolvingVariable<UnsignedInt>("OriginalID");
    addVariableToWrite<UnsignedInt>("OriginalID");
    sorted_id_ = registerDiscreteVariableData<UnsignedInt>("SortedID", particles_bound_, AssignIndex());
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasure(StdVec<Vecd> &pos, StdVec<Real> &Vol)
{
    dv_pos_ = registerStateVariableFrom<Vecd>("Position", pos);
    Vol_ = registerStateVariableDataFrom<Real>("VolumetricMeasure", Vol);
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasureFromReload()
{
    dv_pos_ = registerStateVariableFromReload<Vecd>("Position");
    Vol_ = registerStateVariableDataFromReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBounds(UnsignedInt number_of_particles)
{
    sv_total_real_particles_->setValue(number_of_particles);
    particles_bound_ = number_of_particles;
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBoundsFromReloadXml()
{
    initializeAllParticlesBounds(reload_xml_parser_.Size(reload_xml_parser_.first_element_));
}
//=================================================================================================//
void BaseParticles::increaseParticlesBounds(UnsignedInt extra_size)
{
    particles_bound_ += extra_size;
}
//=================================================================================================//
void BaseParticles::checkEnoughReserve()
{
    if (TotalRealParticles() >= particles_bound_)
    {
        std::cout << "\n Error: Not enough particle reserve! \n"
                  << " Please ensure the particle reserve size when generating particles. \n";
        exit(EXIT_FAILURE);
    }
}
//=================================================================================================//
void BaseParticles::copyFromAnotherParticle(UnsignedInt index, UnsignedInt another_index)
{
    copy_particle_state_(all_state_data_, index, another_index);
}
//=================================================================================================//
UnsignedInt BaseParticles::allocateGhostParticles(UnsignedInt ghost_size)
{
    UnsignedInt ghost_lower_bound = particles_bound_;
    particles_bound_ += ghost_size;
    return ghost_lower_bound;
}
//=================================================================================================//
void BaseParticles::updateGhostParticle(UnsignedInt ghost_index, UnsignedInt index)
{
    copyFromAnotherParticle(ghost_index, index);
    /** For a ghost particle, its sorted id is that of corresponding real particle. */
    sorted_id_[ghost_index] = index;
}
//=================================================================================================//
void BaseParticles::switchToBufferParticle(UnsignedInt index)
{
    UnsignedInt last_real_particle_index = TotalRealParticles() - 1;
    if (index < last_real_particle_index)
    {
        copyFromAnotherParticle(index, last_real_particle_index);
        // update original and sorted_id as well
        std::swap(original_id_[index], original_id_[last_real_particle_index]);
        sorted_id_[original_id_[index]] = index;
    }
    decrementTotalRealParticles();
}
//=================================================================================================//
UnsignedInt BaseParticles::createRealParticleFrom(UnsignedInt index)
{
    UnsignedInt new_original_id = TotalRealParticles();
    original_id_[new_original_id] = new_original_id;
    /** Buffer Particle state copied from real particle. */
    copyFromAnotherParticle(new_original_id, index);
    /** Realize the buffer particle by increasing the number of real particle by one.  */
    incrementTotalRealParticles();
    return new_original_id;
}
//=================================================================================================//
int BaseParticles::getNewBodyPartID()
{
    total_body_parts_++;
    return total_body_parts_;
};
//=================================================================================================//
void BaseParticles::resizeXmlDocForParticles(XmlParser &xml_parser)
{
    UnsignedInt total_elements = xml_parser.Size(xml_parser.first_element_);

    UnsignedInt total_real_particles = TotalRealParticles();
    if (total_elements != total_real_particles)
    {
        xml_parser.resize(xml_parser.first_element_, total_real_particles, "particle");
    }
}
//=================================================================================================//
void BaseParticles::resetTotalRealParticlesFromXmlDoc(XmlParser &xml_parser)
{
    sv_total_real_particles_->setValue(xml_parser.Size(xml_parser.first_element_));
}
//=================================================================================================//
void BaseParticles::readReloadXmlFile(const std::string &filefullpath, const std::string &body_name)
{
    reload_xml_parser_.loadXmlFile(filefullpath);
    // Navigate first_element_ to the matching body element so that downstream
    // functions (initializeAllParticlesBoundsFromReloadXml, registerStateVariableFromReload)
    // continue to work unchanged by iterating its <particle> children directly.
    tinyxml2::XMLElement *body_element = reload_xml_parser_.first_element_->FirstChildElement("body");
    while (body_element != nullptr)
    {
        const char *name_attr = body_element->Attribute("name");
        if (name_attr != nullptr && std::string(name_attr) == body_name)
        {
            reload_xml_parser_.first_element_ = body_element;
            return;
        }
        body_element = body_element->NextSiblingElement("body");
    }
    std::cout << "\n Error: body " << body_name << " not found in reload file: "
              << filefullpath << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
}
//=================================================================================================//
void BaseParticles::writeParticlesToXml(XmlParser &xml_parser, TinyXMLElement *body_element)
{
    // Resize the body element to have the correct number of particle children
    UnsignedInt total_real_particles = TotalRealParticles();
    UnsignedInt total_elements = xml_parser.Size(body_element);

    if (total_elements != total_real_particles)
    {
        xml_parser.resize(body_element, total_real_particles, "particle");
    }

    // Write all evolving variables to the body element's particle children
    OperationOnDataAssemble<DiscreteVariables, WriteParticleVariableToXmlElement>
        write_variable_to_element(body_element);
    write_variable_to_element(evolving_variables_, xml_parser);
}
//=================================================================================================//
void BaseParticles::readParticlesFromXml(XmlParser &xml_parser, TinyXMLElement *body_element)
{
    // Reset total real particles from the body element's particle count
    sv_total_real_particles_->setValue(xml_parser.Size(body_element));

    // Read all evolving variables from the body element's particle children
    OperationOnDataAssemble<DiscreteVariables, ReadParticleVariableFromXmlElement>
        read_variable_from_element(body_element);
    read_variable_from_element(evolving_variables_, this, xml_parser);
}
//=================================================================================================//
} // namespace SPH
