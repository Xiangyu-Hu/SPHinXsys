#include "structural_simulation_class.h"
#include <gtest/gtest.h>

TEST(BernoulliBeam20x, Pressure)
{
    Real scale_stl = 0.001;
    Real end_time = 0.15;

    Real rho_0 = 6.45e3; // Nitinol
    Real poisson = 0.3;
    Real Youngs_modulus = 5e8;
    Real physical_viscosity = Youngs_modulus / 100;
    Real pressure = 1e3;

    /** STL IMPORT PARAMETERS */
    std::string relative_input_path = "./input/"; // path definition for linux
    std::vector<std::string> imported_stl_list = {"bernoulli_beam_20x.stl"};
    std::vector<Vec3d> translation_list = {Vec3d::Zero()};
    std::vector<Real> resolution_list = {10.0 / 6.0};
    SharedPtr<SaintVenantKirchhoffSolid> material = makeShared<SaintVenantKirchhoffSolid>(rho_0, Youngs_modulus, poisson);
    std::vector<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list = {material};

    SharedPtr<TriangleMeshShapeSTL> specimen = makeShared<TriangleMeshShapeSTL>("./input/bernoulli_beam_20x.stl", Vec3d::Zero(), scale_stl, "bernoulli_beam_20x");
    BoundingBoxd fixation = specimen->getBounds();
    fixation.upper_[0] = fixation.lower_[0] + 0.01;

    StructuralSimulationInput input{
        relative_input_path,
        imported_stl_list,
        scale_stl,
        translation_list,
        resolution_list,
        material_model_list,
        {physical_viscosity},
        {}};
    input.body_indices_fixed_constraint_region_ = StdVec<ConstrainedRegionPair>{ConstrainedRegionPair(0, fixation)};
    StdVec<std::array<Real, 2>> pressure_over_time = {
        {Real(0), Real(0)},
        {Real(end_time * 0.1), Real(pressure)},
        {Real(end_time), Real(pressure)}};
    input.surface_pressure_tuple_ = StdVec<PressureTuple>{PressureTuple(0, specimen, Vec3d(0.1, 0.0, 0.1), pressure_over_time)};
    input.particle_relaxation_list_ = StdVec<bool>{true};

    //=================================================================================================//
    StructuralSimulation sim(input);
    sim.runSimulation(end_time);

    Real displ_max = sim.getMaxDisplacement(0);
    Real displ_max_analytical = 4.8e-3; // in mm, absolute max displacement
    EXPECT_NEAR(displ_max, displ_max_analytical, displ_max_analytical * 0.1);
    std::cout << "displ_max: " << displ_max << std::endl;
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
