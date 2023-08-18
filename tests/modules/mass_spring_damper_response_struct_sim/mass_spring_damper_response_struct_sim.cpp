#include "structural_simulation_class.h"
#include <gtest/gtest.h>

TEST(StructuralSimulation, MassSpringDamperResponse)
{
    /** INPUT PARAMETERS */
    Real scale_stl = 0.001 / 4; // diameter of 0.025 m
    Real resolution_mass = 8;
    Real poisson = 0.35;
    Real Youngs_modulus = 1e7;
    Real physical_viscosity = 200;
    // SI setup - designed
    Real rho_0 = 122231;
    Real gravity_force = 10;
    Real spring_coeff = 200;
    Real spring_damper_ratio = 0.05;
    Real end_time = 2;
    /** STL IMPORT PARAMETERS */
    std::string relative_input_path = "./input/"; // path definition for linux
    std::vector<std::string> imported_stl_list = {"ball_mass.stl"};
    std::vector<Vec3d> translation_list = {Vec3d::Zero()};
    std::vector<Real> resolution_list = {resolution_mass};
    SharedPtr<SaintVenantKirchhoffSolid> material = makeShared<SaintVenantKirchhoffSolid>(rho_0, Youngs_modulus, poisson);
    std::vector<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list = {material};
    /** INPUT DECLARATION */
    StructuralSimulationInput input{
        relative_input_path,
        imported_stl_list,
        scale_stl,
        translation_list,
        resolution_list,
        material_model_list,
        {physical_viscosity},
        {}};
    input.non_zero_gravity_ = std::vector<GravityPair>{GravityPair(0, Vec3d(0.0, 0.0, gravity_force))};
    input.spring_damper_tuple_ = {SpringDamperTuple(0, Vec3d(0, 0, spring_coeff), spring_damper_ratio)};

    //=================================================================================================//
    StructuralSimulation sim(input);
    sim.runSimulation(end_time);
    //=================================================================================================//

    StdLargeVec<Vecd> &pos_0 = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos0_;
    StdLargeVec<Vecd> &pos_n = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_;
    Real end_displ = 0.05;

    for (size_t index = 0; index < pos_0.size(); index++)
    {
        Real end_pos = pos_0[index][2] + end_displ;
        EXPECT_NEAR(pos_n[index][2], end_pos, end_displ * 1e-2);
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
