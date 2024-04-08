/**
 * @file 	test_3d_unstructured_mesh_reader.cpp
 * @brief 	This is a test to show the mesh parser for ICEM and Fluent .msh files.
 * @details We consider a 3D channel mesh for ICEM and Fluent.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "test_3d_unstructured_mesh_reader.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh_3d read_mesh_data(mesh_fullpath);
    //----------------------------------------------------------------------
   //	Build up the environment of a SPHSystem.
   //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody wave_block(sph_system, makeShared<AirBody>("AirBody"));
    wave_block.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho, heat_capacity_ratio);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    wave_block.generateParticlesWithReserve<UnstructuredMesh_3d>(ghost_boundary, read_mesh_data);
    GhostCreationFromMesh_3d ghost_creation(wave_block, read_mesh_data, ghost_boundary);

    // Visualization in FVM with date in cell.
    BodyStatesRecordingInMeshToVtu write_real_body_states(wave_block, read_mesh_data);
    write_real_body_states.writeToFile(0);
    return 0;
}