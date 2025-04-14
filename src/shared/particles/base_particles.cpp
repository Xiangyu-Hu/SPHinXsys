#include "base_particles.hpp"

#include "base_body.h"
#include "base_body_part.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "xml_parser.h"

namespace SPH
{
//=================================================================================================//
BaseParticles::BaseParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : sv_total_real_particles_(nullptr),
      particles_bound_(0), original_id_(nullptr), sorted_id_(nullptr),
      dv_pos_(nullptr), Vol_(nullptr), rho_(nullptr), mass_(nullptr),
      sph_body_(sph_body), body_name_(sph_body.getName()),
      base_material_(*base_material),
      restart_xml_parser_("xml_restart", "particles"),
      reload_xml_parser_("xml_particle_reload", "particles")
{
    sph_body.assignBaseParticles(this);
    sv_total_real_particles_ = registerSingularVariable<UnsignedInt>("TotalRealParticles");
}
//=================================================================================================//
SPHAdaptation &BaseParticles::getSPHAdaptation()
{
    return sph_body_.getSPHAdaptation();
}
//=================================================================================================//
void BaseParticles::initializeBasicParticleVariables()
{
    addEvolvingVariable<Vecd>("Position");
    addEvolvingVariable<Real>("VolumetricMeasure");
     //----------------------------------------------------------------------
    //		register non-geometric variables
    //----------------------------------------------------------------------
    rho_ = registerStateVariable<Real>("Density", base_material_.ReferenceDensity());
    mass_ = registerStateVariable<Real>("Mass",
                                        [&](size_t i) -> Real
                                        { return rho_[i] * ParticleVolume(i); });
    //----------------------------------------------------------------------
    //		unregistered variables and data
    //----------------------------------------------------------------------
    original_id_ = registerDiscreteVariable<UnsignedInt>("OriginalID", particles_bound_, AssignIndex());
    addEvolvingVariable<UnsignedInt>("OriginalID");
    addVariableToWrite<UnsignedInt>("OriginalID");
    sorted_id_ = registerDiscreteVariable<UnsignedInt>("SortedID", particles_bound_, AssignIndex());
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasure(StdLargeVec<Vecd> &pos, StdLargeVec<Real> &Vol)
{
    dv_pos_ = registerStateVariableOnlyFrom<Vecd>("Position", pos);
    Vol_ = registerStateVariableFrom<Real>("VolumetricMeasure", Vol);
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasureFromReload()
{
    dv_pos_ = registerStateVariableOnlyFromReload<Vecd>("Position");
    Vol_ = registerStateVariableFromReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBounds(size_t number_of_particles)
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
void BaseParticles::increaseParticlesBounds(size_t extra_size)
{
    particles_bound_ += extra_size;
}
//=================================================================================================//
void BaseParticles::copyFromAnotherParticle(size_t index, size_t another_index)
{
    copy_particle_state_(all_state_data_, index, another_index);
}
//=================================================================================================//
size_t BaseParticles::allocateGhostParticles(size_t ghost_size)
{
    size_t ghost_lower_bound = particles_bound_;
    particles_bound_ += ghost_size;
    return ghost_lower_bound;
}
//=================================================================================================//
void BaseParticles::updateGhostParticle(size_t ghost_index, size_t index)
{
    copyFromAnotherParticle(ghost_index, index);
    /** For a ghost particle, its sorted id is that of corresponding real particle. */
    sorted_id_[ghost_index] = index;
}
//=================================================================================================//
void BaseParticles::switchToBufferParticle(size_t index)
{
    size_t last_real_particle_index = TotalRealParticles() - 1;
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
    size_t new_original_id = TotalRealParticles();
    original_id_[new_original_id] = new_original_id;
    /** Buffer Particle state copied from real particle. */
    copyFromAnotherParticle(new_original_id, index);
    /** Realize the buffer particle by increasing the number of real particle by one.  */
    incrementTotalRealParticles();
    return new_original_id;
}
//=================================================================================================//
void BaseParticles::resizeXmlDocForParticles(XmlParser &xml_parser)
{
    size_t total_elements = xml_parser.Size(xml_parser.first_element_);

    UnsignedInt total_real_particles = TotalRealParticles();
    if (total_elements <= total_real_particles)
    {
        xml_parser.resize(xml_parser.first_element_, total_real_particles, "particle");
    }
}
//=================================================================================================//
void BaseParticles::writeParticlesToXmlForRestart(const std::string &filefullpath)
{
    resizeXmlDocForParticles(restart_xml_parser_);
    write_restart_variable_to_xml_(evolving_variables_, restart_xml_parser_);
    restart_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readParticlesFromXmlForRestart(const std::string &filefullpath)
{
    restart_xml_parser_.loadXmlFile(filefullpath);
    read_restart_variable_from_xml_(evolving_variables_, this, restart_xml_parser_);
}
//=================================================================================================//
void BaseParticles::writeParticlesToXmlForReload(const std::string &filefullpath)
{
    resizeXmlDocForParticles(reload_xml_parser_);
    write_reload_variable_to_xml_(evolving_variables_, reload_xml_parser_);
    reload_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readReloadXmlFile(const std::string &filefullpath)
{
    reload_xml_parser_.loadXmlFile(filefullpath);
}
//=================================================================================================//
} // namespace SPH
