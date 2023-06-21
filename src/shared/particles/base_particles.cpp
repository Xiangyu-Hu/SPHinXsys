#include "base_particles.hpp"

#include "base_body.h"
#include "base_body_part.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "xml_engine.h"

//=====================================================================================================//
namespace SPH
{
//=================================================================================================//
BaseParticles::BaseParticles(SPHBody &sph_body, BaseMaterial *base_material)
    : total_real_particles_(0), real_particles_bound_(0), total_ghost_particles_(0),
      particle_sorting_(*this),
      sph_body_(sph_body), body_name_(sph_body.getName()),
      base_material_(*base_material),
      restart_xml_engine_("xml_restart", "particles"),
      reload_xml_engine_("xml_particle_reload", "particles")
{
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    registerVariable(pos_, "Position");
    registerVariable(Vol_, "VolumetricMeasure");
    //----------------------------------------------------------------------
    //		add particle reload data on geometries
    //----------------------------------------------------------------------
    addVariableToList<Vecd>(variables_to_reload_, "Position");
    addVariableToList<Real>(variables_to_reload_, "VolumetricMeasure");
}
//=================================================================================================//
void BaseParticles::initializeOtherVariables()
{
    //----------------------------------------------------------------------
    //		register non-geometric data
    //----------------------------------------------------------------------
    registerVariable(vel_, "Velocity");
    registerVariable(acc_, "Acceleration");
    registerVariable(acc_prior_, "PriorAcceleration");
    registerVariable(rho_, "Density", base_material_.ReferenceDensity());
    registerVariable(mass_, "MassiveMeasure",
                     [&](size_t i) -> Real
                     { return rho_[i] * Vol_[i]; });
    registerVariable(surface_indicator_, "SurfaceIndicator");
    /**
     *	add basic output particle data
     */
    addVariableToWrite<Vecd>("Velocity");
    /**
     *	add restart output particle data
     */
    addVariableToList<Vecd>(variables_to_restart_, "Position");
    addVariableToList<Vecd>(variables_to_restart_, "Velocity");
    addVariableToList<Vecd>(variables_to_restart_, "Acceleration");
    addVariableToList<Real>(variables_to_restart_, "VolumetricMeasure");
    //----------------------------------------------------------------------
    //		initialize unregistered data
    //----------------------------------------------------------------------
    for (size_t i = 0; i != real_particles_bound_; ++i)
    {
        sorted_id_.push_back(sequence_.size());
        sequence_.push_back(0);
    }
}
//=================================================================================================//
void BaseParticles::addAParticleEntry()
{
    unsorted_id_.push_back(sequence_.size());
    sorted_id_.push_back(sequence_.size());
    sequence_.push_back(0);

    add_particle_data_with_default_value_(all_particle_data_);
}
//=================================================================================================//
void BaseParticles::addBufferParticles(size_t buffer_size)
{
    for (size_t i = 0; i != buffer_size; ++i)
    {
        addAParticleEntry();
    }
    real_particles_bound_ += buffer_size;
}
//=================================================================================================//
void BaseParticles::copyFromAnotherParticle(size_t index, size_t another_index)
{
    updateFromAnotherParticle(index, another_index);
}
//=================================================================================================//
void BaseParticles::updateFromAnotherParticle(size_t index, size_t another_index)
{
    copy_particle_data_(all_particle_data_, index, another_index);
}
//=================================================================================================//
size_t BaseParticles::insertAGhostParticle(size_t index)
{
    total_ghost_particles_ += 1;
    size_t expected_size = real_particles_bound_ + total_ghost_particles_;
    size_t expected_particle_index = expected_size - 1;
    if (expected_size <= pos_.size())
    {
        copyFromAnotherParticle(expected_particle_index, index);
        /** For a ghost particle, its sorted id is that of corresponding real particle. */
        sorted_id_[expected_particle_index] = index;
    }
    else
    {
        addAParticleEntry();
        copyFromAnotherParticle(expected_particle_index, index);
        /** For a ghost particle, its sorted id is that of corresponding real particle. */
        sorted_id_[expected_particle_index] = index;
    }
    return expected_particle_index;
}
//=================================================================================================//
void BaseParticles::switchToBufferParticle(size_t index)
{
    size_t last_real_particle_index = total_real_particles_ - 1;
    if (index < last_real_particle_index)
    {
        updateFromAnotherParticle(index, last_real_particle_index);
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

    // compute derived particle variables
    for (auto &derived_variable : derived_variables_)
    {
        derived_variable->exec();
    }

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
void BaseParticles::resizeXmlDocForParticles(XmlEngine &xml_engine)
{
    size_t total_elements = xml_engine.SizeOfXmlDoc();

    if (total_elements <= total_real_particles_)
    {
        for (size_t i = total_elements; i != total_real_particles_; ++i)
            xml_engine.addElementToXmlDoc("particle");
    }
}
//=================================================================================================//
void BaseParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
{
    resizeXmlDocForParticles(restart_xml_engine_);
    WriteAParticleVariableToXml write_variable_to_xml(restart_xml_engine_, total_real_particles_);
    DataAssembleOperation<loopParticleVariables> loop_variable_namelist;
    loop_variable_namelist(all_particle_data_, variables_to_restart_, write_variable_to_xml);
    restart_xml_engine_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readParticleFromXmlForRestart(std::string &filefullpath)
{
    restart_xml_engine_.loadXmlFile(filefullpath);
    ReadAParticleVariableFromXml read_variable_from_xml(restart_xml_engine_, total_real_particles_);
    DataAssembleOperation<loopParticleVariables> loop_variable_namelist;
    loop_variable_namelist(all_particle_data_, variables_to_restart_, read_variable_from_xml);
}
//=================================================================================================//
void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
{
    resizeXmlDocForParticles(reload_xml_engine_);
    WriteAParticleVariableToXml write_variable_to_xml(reload_xml_engine_, total_real_particles_);
    DataAssembleOperation<loopParticleVariables> loop_variable_namelist;
    loop_variable_namelist(all_particle_data_, variables_to_reload_, write_variable_to_xml);
    reload_xml_engine_.writeToXmlFile(filefullpath);
}
//=================================================================================================//
void BaseParticles::readFromXmlForReloadParticle(std::string &filefullpath)
{
    reload_xml_engine_.loadXmlFile(filefullpath);
    total_real_particles_ = reload_xml_engine_.SizeOfXmlDoc();
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        unsorted_id_.push_back(i);
    };
    resize_particle_data_(all_particle_data_, total_real_particles_);
    ReadAParticleVariableFromXml read_variable_from_xml(reload_xml_engine_, total_real_particles_);
    DataAssembleOperation<loopParticleVariables> loop_variable_namelist;
    loop_variable_namelist(all_particle_data_, variables_to_reload_, read_variable_from_xml);
}
//=================================================================================================//
} // namespace SPH
  //=====================================================================================================//