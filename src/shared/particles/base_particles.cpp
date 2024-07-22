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
    : total_real_particles_(0), real_particles_bound_(0), particles_bound_(0),
      original_id_(nullptr), sorted_id_(nullptr), sequence_(nullptr),
      particle_sorting_(nullptr),
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
}
//=================================================================================================//
void BaseParticles::initializeBasicParticleVariables()
{
    //----------------------------------------------------------------------
    //		register non-geometric variables
    //----------------------------------------------------------------------
    rho_ = registerSharedVariable<Real>("Density", base_material_.ReferenceDensity());
    mass_ = registerSharedVariable<Real>("Mass",
                                         [&](size_t i) -> Real
                                         { return (*rho_)[i] * ParticleVolume(i); });
    //----------------------------------------------------------------------
    //		unregistered variables and data
    //----------------------------------------------------------------------
    original_id_ = registerSharedVariable<size_t>("OriginalID",
                                                  [&](size_t i) -> size_t
                                                  { return i; });
    sorted_id_ = registerSharedVariableFrom<size_t>("SortedID", "OriginalID");
    sequence_ = registerSharedVariable<size_t>("Sequence");
    particle_sorting_ = particle_sort_ptr_keeper_.createPtr<ParticleSorting>(*this);
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasure(StdLargeVec<Vecd> &pos, StdLargeVec<Real> &Vol)
{
    pos_ = registerSharedVariableFrom<Vecd>("Position", pos);
    Vol_ = registerSharedVariableFrom<Real>("VolumetricMeasure", Vol);
    addVariableToReload<Vecd>("Position");
    addVariableToReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::registerPositionAndVolumetricMeasureFromReload()
{
    pos_ = registerSharedVariableFromReload<Vecd>("Position");
    Vol_ = registerSharedVariableFromReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBounds(size_t total_real_particles)
{
    total_real_particles_ = total_real_particles;
    real_particles_bound_ = total_real_particles_;
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
    (*sorted_id_)[ghost_index] = index;
}
//=================================================================================================//
void BaseParticles::switchToBufferParticle(size_t index)
{
    size_t last_real_particle_index = total_real_particles_ - 1;
    if (index < last_real_particle_index)
    {
        copyFromAnotherParticle(index, last_real_particle_index);
        // update original and sorted_id as well
        std::swap((*original_id_)[index], (*original_id_)[last_real_particle_index]);
        (*sorted_id_)[(*original_id_)[index]] = index;
    }
    total_real_particles_ -= 1;
}
//=================================================================================================//
void BaseParticles::createRealParticleFrom(size_t index)
{
    size_t new_original_id = total_real_particles_;
    (*original_id_)[new_original_id] = new_original_id;
    /** Buffer Particle state copied from real particle. */
    copyFromAnotherParticle(new_original_id, index);
    /** Realize the buffer particle by increasing the number of real particle in the body.  */
    total_real_particles_ += 1;
}
//=================================================================================================//
void BaseParticles::writePltFileHeader(std::ofstream &output_file)
{
    output_file << " VARIABLES = \"x\",\"y\",\"z\",\"ID\"";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\""
                    << ",\"" << variable_name << "_z\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };
}
//=================================================================================================//
void BaseParticles::writePltFileParticleData(std::ofstream &output_file, size_t index)
{
    // write particle positions and index first
    Vec3d particle_position = upgradeToVec3d((*pos_)[index]);
    output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
                << index << " ";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *variable->DataField();
        output_file << variable_data[index] << " ";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *variable->DataField();
        Vec3d vector_value = upgradeToVec3d(variable_data[index]);
        output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *variable->DataField();
        output_file << variable_data[index] << " ";
    };
}
//=================================================================================================//
void BaseParticles::writeParticlesToPltFile(std::ofstream &output_file)
{
    writePltFileHeader(output_file);
    output_file << "\n";

    size_t total_real_particles = total_real_particles_;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        writePltFileParticleData(output_file, i);
        output_file << "\n";
    };
}
//=================================================================================================//
void BaseParticles::resizeXmlDocForParticles(XmlParser &xml_parser)
{
    size_t total_elements = xml_parser.Size(xml_parser.first_element_);

    if (total_elements <= total_real_particles_)
    {
        xml_parser.resize(xml_parser.first_element_, total_real_particles_, "particle");
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
} // namespace SPH
