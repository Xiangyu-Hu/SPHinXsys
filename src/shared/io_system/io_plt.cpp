
/**
 * @file 	io_plt.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "io_plt.h"

namespace SPH
{
//=============================================================================================//
void PltEngine::
    writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name)
{
    out_file << "\"" << quantity_name << "\""
             << "   ";
}
//=============================================================================================//
void PltEngine::
    writeAQuantityHeader(std::ofstream &out_file, const Vecd &quantity, const std::string &quantity_name)
{
    for (int i = 0; i != Dimensions; ++i)
        out_file << "\"" << quantity_name << "[" << i << "]\""
                 << "   ";
}
//=============================================================================================//
void PltEngine::writeAQuantity(std::ofstream &out_file, const Real &quantity)
{
    out_file << std::fixed << std::setprecision(9) << quantity << "   ";
}
//=============================================================================================//
void PltEngine::writeAQuantity(std::ofstream &out_file, const Vecd &quantity)
{
    for (int i = 0; i < Dimensions; ++i)
        out_file << std::fixed << std::setprecision(9) << quantity[i] << "   ";
}
//=================================================================================================//
void BodyStatesRecordingToPlt::writePltFileHeader(
    std::ofstream &output_file, ParticleVariables &variables_to_write)
{
    output_file << " VARIABLES = \"x\",\"y\",\"z\",\"ID\"";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write))
    {
        std::string variable_name = variable->Name();
        output_file << ",\"" << variable_name << "_x\""
                    << ",\"" << variable_name << "_y\""
                    << ",\"" << variable_name << "_z\"";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write))
    {
        output_file << ",\"" << variable->Name() << "\"";
    };
}
//=================================================================================================//
void BodyStatesRecordingToPlt::writePltFileParticleData(
    std::ofstream &output_file, ParticleVariables &variables_to_write, Vecd *position, size_t index)
{
    // write particle positions and index first
    Vec3d particle_position = upgradeToVec3d(position[index]);
    output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
                << index << " ";

    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write))
    {
        int *data_field = variable->Data();
        output_file << data_field[index] << " ";
    };

    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write))
    {
        Vecd *data_field = variable->Data();
        Vec3d vector_value = upgradeToVec3d(data_field[index]);
        output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
    };

    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write))
    {
        Real *data_field = variable->Data();
        output_file << data_field[index] << " ";
    };
}
//=============================================================================================//
void BodyStatesRecordingToPlt::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        BaseParticles &particles = body->getBaseParticles();
        ParticleVariables &variables_to_write = particles.VariablesToWrite();
        if (body->checkNewlyUpdated())
        {
            if (state_recording_)
            {
                std::string filefullpath = io_environment_.output_folder_ +
                                           "/SPHBody_" + body->getName() + "_" + sequence + ".plt";
                if (fs::exists(filefullpath))
                {
                    fs::remove(filefullpath);
                }
                std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
                writePltFileHeader(out_file, variables_to_write);
                out_file << "\n";

                Vecd *position = particles.ParticlePositions();
                for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
                {
                    writePltFileParticleData(out_file, variables_to_write, position, i);
                    out_file << "\n";
                };
                out_file.close();
            }
        }
        body->setNotNewlyUpdated();
    }
}
//=============================================================================================//
MeshRecordingToPlt ::MeshRecordingToPlt(SPHSystem &sph_system, BaseMeshField &mesh_field)
    : BaseIO(sph_system), mesh_field_(mesh_field),
      filefullpath_(io_environment_.output_folder_ + "/" + mesh_field.Name() + ".dat") {}
//=============================================================================================//
void MeshRecordingToPlt::writeToFile(size_t iteration_step)
{
    std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
    mesh_field_.writeMeshFieldToPlt(out_file);
    out_file.close();
}
//=================================================================================================//
} // namespace SPH
