/**
 * @file 	2d_gradient.cpp
 * @brief 	test the linear and parabolic reproducing gradient and hessian
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
//	Geometric shapes used in the test
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(container);
        add<GeometricShapeBox>(Transform(translate_to_origin), container);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd container_outer(0.5 * width + boundary_width, 0.5 * height + boundary_width);
        Vecd container(0.5 * width + 2.0 * boundary_width, 0.5 * height);
        Transform translate_to_origin_outer(Vec2d(-boundary_width, -boundary_width) + container_outer);
        Transform translate_to_origin_inner(Vec2d(-boundary_width, 0.0) + container);

        add<GeometricShapeBox>(Transform(translate_to_origin_outer), container_outer);
        subtract<GeometricShapeBox>(Transform(translate_to_origin_inner), container);
    }
};
//----------------------------------------------------------------------
//	Google test items
//----------------------------------------------------------------------
Vec2d random_observation(rand_uniform(0.0, width), rand_uniform(0.0, height));

Matd approximated_gradient(Matd::Zero());
Matd reference_gradient(Matd::Identity());
TEST(LinearGradient, Error)
{
    EXPECT_LT((reference_gradient - approximated_gradient).norm(), 1.0e-5);
    std::cout << "Reference Gradient: " << reference_gradient << " and "
              << "Predicted Gradient: " << approximated_gradient << std::endl;
};

Vec2d first_coefficient(rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0));
VecMat2d second_coefficient{rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0), rand_uniform(-1.0, 1.0)};
class ParabolicProfile : public ReturnFunction<Real>
{
    Vec2d first_coefficient_;
    VecMat2d second_coefficient_;

  public:
    ParabolicProfile() : first_coefficient_(first_coefficient),
                         second_coefficient_(second_coefficient) {};

    Real operator()(const Vec2d &position)
    {
        return first_coefficient_.dot(position) +
               second_coefficient_.dot(vectorizeTensorSquare(position));
    }
};

VecMat2d approximated_hessian = VecMat2d::Zero();
VecMat2d reference_hessian = 2.0 * second_coefficient;
Real hessian_rmse = 0.0;
TEST(Hessian, Error)
{
    EXPECT_LT(hessian_rmse, 6.0e-4);
    std::cout << "Reference Hessian: " << reference_hessian << " and "
              << "Predicted Hessian: " << approximated_hessian
              << ", Hessian RMSE (multi-point): " << hessian_rmse << std::endl;
};

Vec2d approximated_2nd_order_gradient = Vec2d::Zero();
Vec2d reference_2nd_order_gradient =
    first_coefficient +
    Vec2d(2.0 * random_observation[0] * second_coefficient[0] + random_observation[1] * second_coefficient[2],
          2.0 * random_observation[1] * second_coefficient[1] + random_observation[0] * second_coefficient[2]);
TEST(SecondOrderGradient, Error)
{
    EXPECT_LT((reference_2nd_order_gradient - approximated_2nd_order_gradient).norm(), 1.0e-5);
    std::cout << "Reference Second Order Gradient: " << reference_2nd_order_gradient << " and "
              << "Predicted Second Order Gradient: " << approximated_2nd_order_gradient << std::endl;
};

VecMat2d vector_x_second_coefficient{0.34, -0.27, 0.49};
VecMat2d vector_y_second_coefficient{-0.42, 0.31, -0.23};
class VectorXParabolicProfile : public ReturnFunction<Real>
{
    VecMat2d vector_x_second_coefficient_;

  public:
    VectorXParabolicProfile() : vector_x_second_coefficient_(vector_x_second_coefficient) {};

    Real operator()(const Vec2d &position)
    {
        VecMat2d vectorized_position = vectorizeTensorSquare(position);
        return vector_x_second_coefficient_.dot(vectorized_position);
    }
};

class VectorYParabolicProfile : public ReturnFunction<Real>
{
    VecMat2d vector_y_second_coefficient_;

  public:
    VectorYParabolicProfile() : vector_y_second_coefficient_(vector_y_second_coefficient) {};

    Real operator()(const Vec2d &position)
    {
        VecMat2d vectorized_position = vectorizeTensorSquare(position);
        return vector_y_second_coefficient_.dot(vectorized_position);
    }
};

template <typename TransportType, typename... OtherTransportType>
class UnitTestInterfaceModel
{
  public:
    struct Coefficient
    {
        Real operator()(size_t, size_t) const { return 1.0; }
    };

    template <class ExecutionPolicy>
    Coefficient getInterfaceCoeff(const ExecutionPolicy &) const
    {
        return Coefficient();
    }
};

template <typename TransportType, typename... OtherTransportType>
class UnitTestVariableInterfaceModel
{
  public:
    struct Coefficient
    {
        Real operator()(size_t index_i, size_t index_j) const
        {
            Real mode = Real((index_i + 2 * index_j) % 11) / 10.0;
            return 0.8 + 0.4 * mode;
        }
    };

    template <class ExecutionPolicy>
    Coefficient getInterfaceCoeff(const ExecutionPolicy &) const
    {
        return Coefficient();
    }
};

Vec2d approximated_double_curl = Vec2d::Zero();
Vec2d reference_double_curl =
    Vec2d(vector_y_second_coefficient[2] - 2.0 * vector_x_second_coefficient[1],
          -2.0 * vector_y_second_coefficient[0] + vector_x_second_coefficient[2]);
TEST(DoubleCurl, Error)
{
    EXPECT_LT((reference_double_curl - approximated_double_curl).norm(), 1.0e-4);
    std::cout << "Reference Double Curl: " << reference_double_curl << " and "
              << "Predicted Double Curl: " << approximated_double_curl << std::endl;
};

Vec2d reference_double_curl_identity_rhs =
    Vec2d(2.0 * vector_x_second_coefficient[0] + vector_y_second_coefficient[2] -
              (2.0 * vector_x_second_coefficient[0] + 2.0 * vector_x_second_coefficient[1]),
          vector_x_second_coefficient[2] + 2.0 * vector_y_second_coefficient[1] -
              (2.0 * vector_y_second_coefficient[0] + 2.0 * vector_y_second_coefficient[1]));
TEST(DoubleCurlIdentity, Error)
{
    EXPECT_LT((reference_double_curl_identity_rhs - approximated_double_curl).norm(), 1.0e-4);
    std::cout << "Reference grad(div)-laplacian: " << reference_double_curl_identity_rhs << " and "
              << "Predicted Double Curl: " << approximated_double_curl << std::endl;
};

Real interior_double_curl_rmse = 0.0;
Real boundary_double_curl_rmse = 0.0;
Vec2d approximated_double_curl_interface = Vec2d::Zero();
Vec2d approximated_double_curl_interface_variable = Vec2d::Zero();
Real interior_double_curl_rmse_variable_coeff = 0.0;
Real boundary_double_curl_rmse_variable_coeff = 0.0;
TEST(DoubleCurlBoundaryBand, Error)
{
    EXPECT_LT(interior_double_curl_rmse, 5.0e-4);
    EXPECT_LT(boundary_double_curl_rmse, 2.0e-3);
    EXPECT_LT(boundary_double_curl_rmse, 10.0 * interior_double_curl_rmse + 1.0e-6);
    std::cout << "Interior RMSE: " << interior_double_curl_rmse << ", Boundary RMSE: " << boundary_double_curl_rmse
              << std::endl;
};

TEST(DoubleCurlInterfacePath, Error)
{
    EXPECT_LT((reference_double_curl - approximated_double_curl_interface).norm(), 1.0e-4);
    std::cout << "Reference Double Curl: " << reference_double_curl << ", Interface-path Prediction: "
              << approximated_double_curl_interface << std::endl;
};

TEST(DoubleCurlVariableInterfaceCoeff, Behavior)
{
    Real variable_coeff_delta = (approximated_double_curl_interface_variable - approximated_double_curl_interface).norm();
    EXPECT_LT((approximated_double_curl_interface_variable - reference_double_curl).norm(), 3.0e-1);
    EXPECT_GT(boundary_double_curl_rmse_variable_coeff, 10.0 * boundary_double_curl_rmse);
    EXPECT_LT(interior_double_curl_rmse_variable_coeff, 2.5e-1);
    EXPECT_LT(boundary_double_curl_rmse_variable_coeff, 2.5e-1);
    std::cout << "Variable-coeff Double Curl: " << approximated_double_curl_interface_variable
              << ", delta vs constant-coeff: " << variable_coeff_delta
              << ", interior RMSE: " << interior_double_curl_rmse_variable_coeff
              << ", boundary RMSE: " << boundary_double_curl_rmse_variable_coeff << std::endl;
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
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_fluid_particles(water_block);
    random_fluid_particles.exec(0.25); // randomize particle to avoid the symmetry

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vec2d> observation_location = {random_observation};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    ObserverBody profile_observer(sph_system, "ProfileObserver");
    ObserverBody hessian_observer(sph_system, "HessianObserver");
    Real offset = 1.5 * particle_spacing;
    StdVec<Vec2d> profile_points = {
        Vec2d(0.25 * width, 0.25 * height), Vec2d(0.75 * width, 0.25 * height),
        Vec2d(0.25 * width, 0.75 * height), Vec2d(0.75 * width, 0.75 * height),
        Vec2d(offset, 0.5 * height),         Vec2d(width - offset, 0.5 * height),
        Vec2d(0.5 * width, offset),          Vec2d(0.5 * width, height - offset)};
    profile_observer.generateParticles<ObserverParticles>(profile_points);
    StdVec<Vec2d> hessian_points = {
        Vec2d(0.2 * width, 0.2 * height), Vec2d(0.8 * width, 0.2 * height),
        Vec2d(0.2 * width, 0.8 * height), Vec2d(0.8 * width, 0.8 * height),
        Vec2d(0.5 * width, 0.2 * height), Vec2d(0.5 * width, 0.8 * height),
        Vec2d(0.2 * width, 0.5 * height), Vec2d(0.8 * width, 0.5 * height)};
    hessian_observer.generateParticles<ObserverParticles>(hessian_points);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall});
    Contact<> fluid_observer_contact(fluid_observer, {&water_block});
    Contact<> profile_observer_contact(profile_observer, {&water_block});
    Contact<> hessian_observer_contact(hessian_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_water_block_inner(water_block_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_water_wall_contact(water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_fluid_observer_contact(fluid_observer_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_profile_observer_contact(profile_observer_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_hessian_observer_contact(hessian_observer_contact);
    UnitTestInterfaceModel<Real> interface_model;
    StdVec<UnitTestInterfaceModel<Real> *> interface_models = {&interface_model};
    UnitTestVariableInterfaceModel<Real> variable_interface_model;
    StdVec<UnitTestVariableInterfaceModel<Real> *> variable_interface_models = {&variable_interface_model};

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        linear_correction_matrix(DynamicsArgs(water_block_inner, 0.0), water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>, Contact<Vecd>>>
        position_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Position")),
            DynamicsArgs(water_wall_contact, std::string("Position")));
    ObservedQuantityRecording<MainExecutionPolicy, Matd, RestoringCorrection>
        observed_position_gradient("PositionGradient", fluid_observer_contact);

    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>, Contact<>>>
        displacement_matrix_gradient(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>, Contact<>>>
        hessian_correction_matrix(DynamicsArgs(water_block_inner, 0.0), water_wall_contact);

    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ParabolicProfile>>>
        water_block_initial_condition(water_block, "Phi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<ParabolicProfile>>>
        wall_initial_condition(wall, "Phi");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<VectorXParabolicProfile>>>
        water_block_vector_x_initial_condition(water_block, "PhiX");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<VectorXParabolicProfile>>>
        wall_vector_x_initial_condition(wall, "PhiX");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<VectorYParabolicProfile>>>
        water_block_vector_y_initial_condition(water_block, "PhiY");
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, SpatialDistribution<VectorYParabolicProfile>>>
        wall_vector_y_initial_condition(wall, "PhiY");
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        variable_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>, Contact<UnitTestInterfaceModel<Real>>>>
        variable_hessian(
            DynamicsArgs(water_block_inner, std::string("Phi")),
            DynamicsArgs(water_wall_contact, std::string("Phi"), interface_models));
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian("PhiHessian", fluid_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_multi("PhiHessian", hessian_observer_contact);

    InteractionDynamicsCK<MainExecutionPolicy, SecondOrderGradient<Inner<Real>>>
        variable_2nd_order_gradient(
            DynamicsArgs(water_block_inner, std::string("Phi")));
    ObservedQuantityRecording<MainExecutionPolicy, Vec2d, RestoringCorrection>
        observed_2nd_order_gradient("PhiGradient", fluid_observer_contact);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        variable_x_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("PhiX")),
            DynamicsArgs(water_wall_contact, std::string("PhiX")));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>, Contact<Real>>>
        variable_y_linear_gradient(
            DynamicsArgs(water_block_inner, std::string("PhiY")),
            DynamicsArgs(water_wall_contact, std::string("PhiY")));
    std::string phi_x_interface_name = "PhiX";
    std::string phi_y_interface_name = "PhiY";
    InteractionDynamicsCK<
        MainExecutionPolicy,
        Hessian<Inner<Real>, Contact<UnitTestInterfaceModel<Real>>>>
        variable_x_hessian_interface(
            DynamicsArgs(water_block_inner, phi_x_interface_name),
            DynamicsArgs(water_wall_contact, phi_x_interface_name, interface_models));
    InteractionDynamicsCK<
        MainExecutionPolicy,
        Hessian<Inner<Real>, Contact<UnitTestInterfaceModel<Real>>>>
        variable_y_hessian_interface(
            DynamicsArgs(water_block_inner, phi_y_interface_name),
            DynamicsArgs(water_wall_contact, phi_y_interface_name, interface_models));
    InteractionDynamicsCK<
        MainExecutionPolicy,
        Hessian<Inner<Real>, Contact<UnitTestVariableInterfaceModel<Real>>>>
        variable_x_hessian_variable_interface(
            DynamicsArgs(water_block_inner, phi_x_interface_name),
            DynamicsArgs(water_wall_contact, phi_x_interface_name, variable_interface_models));
    InteractionDynamicsCK<
        MainExecutionPolicy,
        Hessian<Inner<Real>, Contact<UnitTestVariableInterfaceModel<Real>>>>
        variable_y_hessian_variable_interface(
            DynamicsArgs(water_block_inner, phi_y_interface_name),
            DynamicsArgs(water_wall_contact, phi_y_interface_name, variable_interface_models));
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_x_profile("PhiXHessian", profile_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_y_profile("PhiYHessian", profile_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_x_interface("PhiXHessian", fluid_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_y_interface("PhiYHessian", fluid_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_x_variable_interface("PhiXHessian", fluid_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_y_variable_interface("PhiYHessian", fluid_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_x_profile_variable_interface("PhiXHessian", profile_observer_contact);
    ObservedQuantityRecording<MainExecutionPolicy, VecMat2d, RestoringCorrection>
        observed_hessian_y_profile_variable_interface("PhiYHessian", profile_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    update_water_block_inner.exec();
    update_water_wall_contact.exec();
    update_fluid_observer_contact.exec();
    update_profile_observer_contact.exec();
    update_hessian_observer_contact.exec();

    water_block_initial_condition.exec();
    wall_initial_condition.exec();
    water_block_vector_x_initial_condition.exec();
    wall_vector_x_initial_condition.exec();
    water_block_vector_y_initial_condition.exec();
    wall_vector_y_initial_condition.exec();

    linear_correction_matrix.exec();
    position_linear_gradient.exec();
    observed_position_gradient.writeToFile(0);
    approximated_gradient = *observed_position_gradient.getObservedQuantity();

    displacement_matrix_gradient.exec();
    hessian_correction_matrix.exec();
    variable_linear_gradient.exec();
    variable_hessian.exec();
    observed_hessian.writeToFile(0);
    approximated_hessian = *observed_hessian.getObservedQuantity();
    observed_hessian_multi.writeToFile(0);
    VecMat2d *hessian_multi = observed_hessian_multi.getObservedQuantity();
    Real hessian_error_acc = 0.0;
    for (size_t i = 0; i < hessian_points.size(); ++i)
    {
        hessian_error_acc += (hessian_multi[i] - reference_hessian).squaredNorm();
    }
    hessian_rmse = std::sqrt(hessian_error_acc / Real(hessian_points.size()));

    variable_2nd_order_gradient.exec();
    observed_2nd_order_gradient.writeToFile(0);
    approximated_2nd_order_gradient = *observed_2nd_order_gradient.getObservedQuantity();

    variable_x_linear_gradient.exec();
    variable_y_linear_gradient.exec();
    variable_x_hessian_interface.exec();
    variable_y_hessian_interface.exec();
    observed_hessian_x_interface.writeToFile(0);
    observed_hessian_y_interface.writeToFile(0);
    VecMat2d approximated_hessian_x_interface = *observed_hessian_x_interface.getObservedQuantity();
    VecMat2d approximated_hessian_y_interface = *observed_hessian_y_interface.getObservedQuantity();
    approximated_double_curl_interface =
        Vec2d(0.5 * approximated_hessian_y_interface[2] - approximated_hessian_x_interface[1],
              -approximated_hessian_y_interface[0] + 0.5 * approximated_hessian_x_interface[2]);
    approximated_double_curl = approximated_double_curl_interface;

    observed_hessian_x_profile.writeToFile(0);
    observed_hessian_y_profile.writeToFile(0);

    VecMat2d *hessian_x_profile = observed_hessian_x_profile.getObservedQuantity();
    VecMat2d *hessian_y_profile = observed_hessian_y_profile.getObservedQuantity();
    Real interior_error_acc = 0.0;
    Real boundary_error_acc = 0.0;
    for (size_t i = 0; i < profile_points.size(); ++i)
    {
        Vec2d profile_double_curl = Vec2d(0.5 * hessian_y_profile[i][2] - hessian_x_profile[i][1],
                                          -hessian_y_profile[i][0] + 0.5 * hessian_x_profile[i][2]);
        Real squared_error = (profile_double_curl - reference_double_curl).squaredNorm();
        if (i < 4)
        {
            interior_error_acc += squared_error;
        }
        else
        {
            boundary_error_acc += squared_error;
        }
    }
    interior_double_curl_rmse = std::sqrt(interior_error_acc / 4.0);
    boundary_double_curl_rmse = std::sqrt(boundary_error_acc / 4.0);

    variable_x_hessian_variable_interface.exec();
    variable_y_hessian_variable_interface.exec();
    observed_hessian_x_variable_interface.writeToFile(0);
    observed_hessian_y_variable_interface.writeToFile(0);
    VecMat2d approximated_hessian_x_variable_interface = *observed_hessian_x_variable_interface.getObservedQuantity();
    VecMat2d approximated_hessian_y_variable_interface = *observed_hessian_y_variable_interface.getObservedQuantity();
    approximated_double_curl_interface_variable =
        Vec2d(0.5 * approximated_hessian_y_variable_interface[2] - approximated_hessian_x_variable_interface[1],
              -approximated_hessian_y_variable_interface[0] + 0.5 * approximated_hessian_x_variable_interface[2]);

    observed_hessian_x_profile_variable_interface.writeToFile(0);
    observed_hessian_y_profile_variable_interface.writeToFile(0);

    VecMat2d *hessian_x_profile_variable_interface = observed_hessian_x_profile_variable_interface.getObservedQuantity();
    VecMat2d *hessian_y_profile_variable_interface = observed_hessian_y_profile_variable_interface.getObservedQuantity();
    Real interior_error_acc_variable_coeff = 0.0;
    Real boundary_error_acc_variable_coeff = 0.0;
    for (size_t i = 0; i < profile_points.size(); ++i)
    {
        Vec2d profile_double_curl_variable = Vec2d(
            0.5 * hessian_y_profile_variable_interface[i][2] - hessian_x_profile_variable_interface[i][1],
            -hessian_y_profile_variable_interface[i][0] + 0.5 * hessian_x_profile_variable_interface[i][2]);
        Real squared_error_variable = (profile_double_curl_variable - reference_double_curl).squaredNorm();
        if (i < 4)
        {
            interior_error_acc_variable_coeff += squared_error_variable;
        }
        else
        {
            boundary_error_acc_variable_coeff += squared_error_variable;
        }
    }
    interior_double_curl_rmse_variable_coeff = std::sqrt(interior_error_acc_variable_coeff / 4.0);
    boundary_double_curl_rmse_variable_coeff = std::sqrt(boundary_error_acc_variable_coeff / 4.0);

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
