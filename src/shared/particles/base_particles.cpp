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
    : sv_total_real_particles_(nullptr), real_particles_bound_(0), particles_bound_(0),
      original_id_(nullptr), sorted_id_(nullptr),
      pos_(nullptr), Vol_(nullptr), rho_(nullptr), mass_(nullptr),
      sph_body_(sph_body), body_name_(sph_body.getName()),
      base_material_(*base_material),
      restart_xml_parser_("xml_restart", "particles"),
      reload_xml_parser_("xml_particle_reload", "particles"),
      copy_particle_state_(all_state_data_),
      write_restart_variable_to_xml_(variables_to_restart_, restart_xml_parser_),
      write_reload_variable_to_xml_(variables_to_reload_, reload_xml_parser_),
      read_restart_variable_from_xml_(variables_to_restart_, restart_xml_parser_)
{
    sph_body.assignBaseParticles(this);
    sv_total_real_particles_ = registerSingularVariable<UnsignedInt>("TotalRealParticles");
}
//=================================================================================================//
void BaseParticles::initializeBasicParticleVariables()
{
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
    original_id_ = registerDiscreteVariable<UnsignedInt>("OriginalID", particles_bound_, getAssignIndex());
    sorted_id_ = registerDiscreteVariable<UnsignedInt>("SortedID", particles_bound_, getAssignIndex());
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasure(StdLargeVec<Vecd> &pos, StdLargeVec<Real> &Vol)
{
    pos_ = registerStateVariableFrom<Vecd>("Position", pos);
    Vol_ = registerStateVariableFrom<Real>("VolumetricMeasure", Vol);
    addVariableToReload<Vecd>("Position");
    addVariableToReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasureFromReload()
{
    pos_ = registerStateVariableFromReload<Vecd>("Position");
    Vol_ = registerStateVariableFromReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBounds(size_t number_of_particles)
{
    UnsignedInt *total_real_particles = sv_total_real_particles_->Data();
    *total_real_particles = number_of_particles;
    real_particles_bound_ = number_of_particles;
    particles_bound_ = real_particles_bound_;
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBoundsFromReloadXml()
{
    initializeAllParticlesBounds(reload_xml_parser_.Size(reload_xml_parser_.first_element_));
}
//=================================================================================================//
void BaseParticles::increaseAllParticlesBounds(size_t buffer_size)
{
    real_particles_bound_ += buffer_size;
    particles_bound_ += buffer_size;
}
//=================================================================================================//
void BaseParticles::copyFromAnotherParticle(size_t index, size_t another_index)
{
    copy_particle_state_(index, another_index);
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
void BaseParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
{
    resizeXmlDocForParticles(restart_xml_parser_);
    write_restart_variable_to_xml_();
    restart_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readParticleFromXmlForRestart(std::string &filefullpath)
{
    restart_xml_parser_.loadXmlFile(filefullpath);
    read_restart_variable_from_xml_(this);
}
//=================================================================================================//
void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
{
    resizeXmlDocForParticles(reload_xml_parser_);
    write_reload_variable_to_xml_();
    reload_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
XmlParser &BaseParticles::readReloadXmlFile(const std::string &filefullpath)
{
    is_reload_file_read_ = true;
    reload_xml_parser_.loadXmlFile(filefullpath);
    return reload_xml_parser_;
}
//=================================================================================================//
BaseParticles::CreateRealParticleFrom::CreateRealParticleFrom(BaseParticles &base_particles)
    : copyable_states_()
{
    OperationBetweenDataAssembles<ParticleVariables, DiscreteVariableArrays, DiscreteVariableArraysInitialization>
        initialize_discrete_variable_array(base_particles.VariablesToSort(), copyable_states_);
    initialize_discrete_variable_array();
}
//=================================================================================================//
} // namespace SPH
