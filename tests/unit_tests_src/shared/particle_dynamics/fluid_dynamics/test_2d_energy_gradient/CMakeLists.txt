cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake)
set(CMAKE_VERBOSE_MAKEFILE on)
PROJECT("test_2d_energy_gradient")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/include
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/bodies
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/adaptations
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/common
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/common
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/geometries
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/meshes
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/particles
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/particle_generator
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/for_2D_build/particle_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/kernels
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/meshes
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/execution
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/geometries
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/physical_closure/materials
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_generator
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particles
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/tools
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_neighborhood
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/physical_closure
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/sphinxsys_system
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/io_system
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/simbody_sphinxsys
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/bodies/complex_bodies
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/body_relations
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/general_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/continuum_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/fluid_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/solid_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/diffusion_reaction_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/fluid_dynamics/eulerian_fluid_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/fluid_dynamics/boundary_condition
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/fluid_dynamics/non_newtonian_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/general_dynamics/domian_bouding
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/general_dynamics/surface_indication
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/general_dynamics/boundary_condition
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/relax_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/shared_ck/particle_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/shared_ck/particle_dynamics/configuration_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/shared_ck/neighborhood
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/shared_ck/smoothing_kernel
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/shared_ck/body_relation
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/external_force
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/active_muscle_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/configuration_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/solid_dynamics/contact_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/diffusion_optimization_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/dissipation_dynamics
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/electro_physiology
  /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/regression_test
)

include_directories(/opt/homebrew/include)
link_directories(/opt/homebrew/lib)
link_directories(/Users/kursat/Desktop/repo/sphinxsys-vscode/build/src)
include_directories(/opt/homebrew/include)
include_directories(/opt/homebrew/include/eigen3)
include_directories(/opt/simbody/include/simbody)
include_directories(/Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/sphinxsys_system)
link_directories(/opt/simbody/lib)
include_directories(/opt/simbody/include)

set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")
file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})
aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(test_2d_energy_gradient ${DIR_SRCS})
add_test(NAME test_2d_energy_gradient COMMAND test_2d_energy_gradient --state_recording=${TEST_STATE_RECORDING}
WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
set_target_properties(test_2d_energy_gradient PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
target_sources(test_2d_energy_gradient PRIVATE
    /Users/kursat/Desktop/repo/sphinxsys-vscode/src/shared/particle_dynamics/fluid_dynamics/energy_gradient.cpp
)

find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
find_package(TBB REQUIRED)
include_directories(${TBB_INCLUDE_DIRS})
target_link_libraries(test_2d_energy_gradient 
    sphinxsys_2d 
    SimTKsimbody 
    SimTKmath 
    SimTKcommon 
    gtest 
    ${Boost_LIBRARIES}
    TBB::tbb TBB::tbbmalloc
)
