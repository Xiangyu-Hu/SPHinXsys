#include "base_particles.hpp"

#include "base_body.h"
#include "base_body_part.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "xml_parser.h"

//=====================================================================================================//
namespace SPH
{
//=================================================================================================//
BaseParticles::BaseParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : total_real_particles_(0), real_particles_bound_(0), particles_bound_(0),
      particle_sorting_(*this),
      sph_body_(sph_body), body_name_(sph_body.getName()),
      base_material_(*base_material),
      restart_xml_parser_("xml_restart", "particles"),
      reload_xml_parser_("xml_particle_reload", "particles"),
      resize_particles_(all_particle_data_),
      copy_particle_data_(all_particle_data_),
      write_restart_variable_to_xml_(variables_to_restart_, restart_xml_parser_),
      write_reload_variable_to_xml_(variables_to_reload_, reload_xml_parser_),
      read_restart_variable_from_xml_(variables_to_restart_, restart_xml_parser_),
      read_reload_variable_from_xml_(variables_to_reload_, reload_xml_parser_)
{
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    registerVariable(pos_, "Position");
    registerVariable(Vol_, "VolumetricMeasure");
    //----------------------------------------------------------------------
    //		add particle reload data on geometries
    //----------------------------------------------------------------------
    addVariableToReload<Vecd>("Position");
    addVariableToReload<Real>("VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeOtherVariables()
{
    //----------------------------------------------------------------------
    //		register non-geometric data
    //----------------------------------------------------------------------
    registerVariable(rho_, "Density", base_material_.ReferenceDensity());
    registerVariable(mass_, "Mass",
                     [&](size_t i) -> Real
                     { return rho_[i] * ParticleVolume(i); });
    //----------------------------------------------------------------------
    //		initialize unregistered data
    //----------------------------------------------------------------------
    for (size_t i = 0; i != particles_bound_; ++i)
    {
        sorted_id_.push_back(sequence_.size());
        sequence_.push_back(0);
    }
}
//=================================================================================================//
void BaseParticles::initializeAllParticlesBounds()
{
    real_particles_bound_ = total_real_particles_;
    particles_bound_ = real_particles_bound_;
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
    copy_particle_data_(index, another_index);
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
    size_t last_real_particle_index = total_real_particles_ - 1;
    if (index < last_real_particle_index)
    {
        copyFromAnotherParticle(index, last_real_particle_index);
        // update unsorted and sorted_id as well
        std::swap(unsorted_id_[index], unsorted_id_[last_real_particle_index]);
        sorted_id_[unsorted_id_[index]] = index;
    }
    total_real_particles_ -= 1;
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
    Vec3d particle_position = upgradeToVec3d(pos_[index]);
    output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
                << index << " ";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *(std::get<type_index_int>(all_particle_data_)[variable->IndexInContainer()]);
        output_file << variable_data[index] << " ";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *(std::get<type_index_Vecd>(all_particle_data_)[variable->IndexInContainer()]);
        Vec3d vector_value = upgradeToVec3d(variable_data[index]);
        output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *(std::get<type_index_Real>(all_particle_data_)[variable->IndexInContainer()]);
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
void BaseParticles::writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface &surface_particles)
{
    size_t total_surface_particles = surface_particles.body_part_particles_.size();

    // write current/final particle positions first
    output_file << "   <Points>\n";
    output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
    output_file << "    ";
    for (size_t i = 0; i != total_surface_particles; ++i)
    {
        size_t particle_i = surface_particles.body_part_particles_[i];
        Vec3d particle_position = upgradeToVec3d(pos_[particle_i]);
        output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
    }
    output_file << std::endl;
    output_file << "    </DataArray>\n";
    output_file << "   </Points>\n";

    // write header of particles data
    output_file << "   <PointData  Vectors=\"vector\">\n";

    // write sorted particles ID
    output_file << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_file << "    ";
    for (size_t i = 0; i != total_surface_particles; ++i)
    {
        size_t particle_i = surface_particles.body_part_particles_[i];
        output_file << particle_i << " ";
    }
    output_file << std::endl;
    output_file << "    </DataArray>\n";

    // write unsorted particles ID
    output_file << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_file << "    ";
    for (size_t i = 0; i != total_surface_particles; ++i)
    {
        size_t particle_i = surface_particles.body_part_particles_[i];
        output_file << unsorted_id_[particle_i] << " ";
    }
    output_file << std::endl;
    output_file << "    </DataArray>\n";

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write_))
    {
        StdLargeVec<Matd> &variable_data = *(std::get<type_index_Matd>(all_particle_data_)[variable->IndexInContainer()]);
        output_file << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_file << "    ";
        for (size_t i = 0; i != total_surface_particles; ++i)
        {
            size_t particle_i = surface_particles.body_part_particles_[i];
            Mat3d matrix_value = upgradeToMat3d(variable_data[particle_i]);
            for (int k = 0; k != 3; ++k)
            {
                Vec3d col_vector = matrix_value.col(k);
                output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
            }
        }
        output_file << std::endl;
        output_file << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *(std::get<type_index_Vecd>(all_particle_data_)[variable->IndexInContainer()]);
        output_file << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_file << "    ";
        for (size_t i = 0; i != total_surface_particles; ++i)
        {
            size_t particle_i = surface_particles.body_part_particles_[i];
            Vec3d vector_value = upgradeToVec3d(variable_data[particle_i]);
            output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_file << std::endl;
        output_file << "    </DataArray>\n";
    }

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *(std::get<type_index_Real>(all_particle_data_)[variable->IndexInContainer()]);
        output_file << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_file << "    ";
        for (size_t i = 0; i != total_surface_particles; ++i)
        {
            size_t particle_i = surface_particles.body_part_particles_[i];
            output_file << std::fixed << std::setprecision(9) << variable_data[particle_i] << " ";
        }
        output_file << std::endl;
        output_file << "    </DataArray>\n";
    }

    // write integers
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *(std::get<type_index_int>(all_particle_data_)[variable->IndexInContainer()]);
        output_file << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
        output_file << "    ";
        for (size_t i = 0; i != total_surface_particles; ++i)
        {
            size_t particle_i = surface_particles.body_part_particles_[i];
            output_file << std::fixed << std::setprecision(9) << variable_data[particle_i] << " ";
        }
        output_file << std::endl;
        output_file << "    </DataArray>\n";
    }
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
    write_restart_variable_to_xml_(all_particle_data_);
    restart_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readParticleFromXmlForRestart(std::string &filefullpath)
{
    restart_xml_parser_.loadXmlFile(filefullpath);
    read_restart_variable_from_xml_(all_particle_data_);
}
//=================================================================================================//
void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
{
    resizeXmlDocForParticles(reload_xml_parser_);
    write_reload_variable_to_xml_(all_particle_data_);
    reload_xml_parser_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readFromXmlForReloadParticle(std::string &filefullpath)
{
    reload_xml_parser_.loadXmlFile(filefullpath);
    total_real_particles_ = reload_xml_parser_.Size(reload_xml_parser_.first_element_);
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        unsorted_id_.push_back(i);
    };
    resize_particles_(total_real_particles_);
    read_reload_variable_from_xml_(all_particle_data_);
}
//=================================================================================================//
} // namespace SPH
  //=====================================================================================================//
