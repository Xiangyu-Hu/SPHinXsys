/**
 * @file 	2d_interpolation.cpp
 * @brief 	test the linear reproducing approximation of interpolation.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real mu_f = 1.0e-1;      /**< Viscosity. */
Real U_max = 1.0;        // make sure the maximum anticipated speed
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real width = 1.0;
Real height = 0.5;
Real particle_spacing = 0.01;
Real boundary_width = particle_spacing * 4; // boundary width
//----------------------------------------------------------------------
//	Google test item.
//----------------------------------------------------------------------
Vecd approximated_coordinate(0.0, 0.0);
Vecd reference_coordinate(rand_uniform(0.0, width), rand_uniform(0.0, height));
TEST(Interpolation, Error)
{
    EXPECT_LT((approximated_coordinate - reference_coordinate).norm(), 1.0e-6);
    std::cout << "Reference Coordinate: " << reference_coordinate << " and "
              << "Predicted Coordinate: " << approximated_coordinate << std::endl;
};

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(scaled_container);
        add<GeometricShapeBox>(Transform(translate_to_origin), scaled_container);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-boundary_width * 2, -boundary_width * 2),
                                     Vecd(width + boundary_width * 2, height + boundary_width * 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_airfoil_particles(water_block);
    random_airfoil_particles.exec(0.5);

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {reference_coordinate};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Contact<> fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    ObservedQuantityRecording<MainExecutionPolicy, Vecd, RestoringCorrection>
        fluid_observer_position("Position", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_cell_linked_list.exec();
    fluid_observer_contact_relation.exec();
    fluid_observer_position.writeToFile(0);
    approximated_coordinate = *fluid_observer_position.getObservedQuantity();
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
