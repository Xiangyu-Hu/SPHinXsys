#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "polar_decomposition_3x3.h"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;

void aortic_valve_pulsatile_fsi(size_t res_factor);

int main(int argc, char *argv[]) { aortic_valve_pulsatile_fsi(1); }

struct parameters
{
    const Real scale = 0.001;
    const Real stent_diameter = 25.0 * scale;
    const Real leaflet_thickness = 1.0 * scale;

    const Real fluid_upstream_length = 1.5 * stent_diameter;
    const Real fluid_downstream_length = 1.5 * stent_diameter;

    // FLUID MATERIAL PARAMETERS
    Real rho0_f = 1100.0; /**< Density. */
    const Real U_max = 1.0;
    Real c_f = 10.0 * U_max; /**< Speed of sound. */
    Real mu_f = 3.6e-3;      /**< Infinite viscosity. */

    // TIME PARAMETERS
    const Real end_time = 0.01;

    // SOLID MATERIAL PARAMETERS
    Real rho0_s = 1000.0; /**< Reference density.*/
    Real bulk_modulus = 1.6e6 * std::pow(0.15 * scale / leaflet_thickness, 3);
    Real shear_modulus = 2.4e6 * std::pow(0.15 * scale / leaflet_thickness, 3);

    // FILEPATH
    std::string leaflet_1 = "leaflet_1";
    std::string leaflet_2 = "leaflet_2";
    std::string leaflet_3 = "leaflet_3";
    const std::vector<std::string> leaflet_names = {"leaflet_1", "leaflet_2", "leaflet_3"};
    const std::vector<std::string> seperate_leaflet_files = {
        "input/simplified_leaflets/leaflet_1.stl",
        "input/simplified_leaflets/leaflet_2.stl",
        "input/simplified_leaflets/leaflet_3.stl"};
};

//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    explicit LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct RightInflowPressure
{
    const Real outlet_pressure = 0.0;

    template <class BoundaryConditionType>
    explicit RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real physical_time)
    {
        return outlet_pressure;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    const Real time_cycle = 1.0;
    const Real radius = 0.5 * parameters{}.stent_diameter;

    // Define inlet average velocity
    const std::vector<Real> popt = {
        3.89614517e+01,
        -9.10133856e+02,
        7.02814566e+03,
        -3.71915499e+00,
        2.18342609e+02,
        -8.19949931e+03,
        1.87965868e+01,
        -9.84301025e+02,
        3.34032930e+03,
        -1.58196177e+01,
        -6.11735822e+01,
        -1.39487710e+03,
        -3.11323367e+00,
        -1.29113904e+02,
        -7.40179967e+02,
        2.50150794e+01,
        1.12384813e+02,
        1.18419142e+02};

    inline Real f1(Real t)
    {
        return popt[0] * t + popt[1] * t * t + popt[2] * t * t * t;
    };
    inline Real f2(Real t)
    {
        return popt[3] * (t - 0.06) + popt[4] * (t - 0.06) * (t - 0.06) + popt[5] * (t - 0.06) * (t - 0.06) * (t - 0.06) + f1(t);
    };
    inline Real f3(Real t)
    {
        return popt[6] * (t - 0.17) + popt[7] * (t - 0.17) * (t - 0.17) + popt[8] * (t - 0.17) * (t - 0.17) * (t - 0.17) + f2(t);
    };
    inline Real f4(Real t)
    {
        return popt[9] * (t - 0.266) + popt[10] * (t - 0.266) * (t - 0.266) + popt[11] * (t - 0.266) * (t - 0.266) * (t - 0.266) + f3(t);
    };
    inline Real f5(Real t)
    {
        return popt[12] * (t - 0.42) + popt[13] * (t - 0.42) * (t - 0.42) + popt[14] * (t - 0.42) * (t - 0.42) * (t - 0.42) + f4(t);
    };
    inline Real f6(Real t)
    {
        return popt[15] * (t - 0.93) + popt[16] * (t - 0.93) * (t - 0.93) + popt[17] * (t - 0.93) * (t - 0.93) * (t - 0.93) + f5(t);
    };
    inline Real flow_rate(Real run_time)
    {
        Real q;
        // Subtract the time for solid simulation
        Real t = run_time - std::floor(run_time / time_cycle);
        if (t <= 0.06)
            q = f1(t);
        else if (t <= 0.17)
            q = f2(t);
        else if (t <= 0.266)
            q = f3(t);
        else if (t <= 0.42)
            q = f4(t);
        else if (t <= 0.93)
            q = f5(t);
        else
            q = f6(t);
        // from ml/min to m^3/s
        return 1.e-3 / 60. * q;
    };

    template <class BoundaryConditionType>
    explicit InflowVelocity(BoundaryConditionType &boundary_condition) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Real q = flow_rate(current_time);
        Real u_ave_ = q / (Pi * radius * radius);

        // Position and velocity in local frame: normal in x direction
        Real r_sqaure = position.y() * position.y() + position.z() * position.z();
        Real u = SMAX(2.0 * u_ave_ * (1.0 - r_sqaure / (radius * radius)), 0.0);

        return u * Vec3d::UnitX();
    }
};

struct solid_object_input
{
    // mandatory input for construction
    std::string name;
    SPH::SharedPtr<SPH::Shape> mesh;
    SPH::NeoHookeanSolid *material;
    double resolution;
    double length_scale;
};

//------------------------------------------------------------------------------
void relax_solid(RealBody &body, BaseInnerRelation &inner)
{
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepInner relaxation_step_inner(inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    body.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
        }
    }
    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
}

class SolidObject
{
  private:
    SolidBody body_;
    std::unique_ptr<InnerRelation> inner_relation_;

  public:
    SolidObject(SPHSystem &system, const solid_object_input &input)
        : body_(system, input.mesh, input.name)
    {
        body_.defineBodyLevelSetShape()->correctLevelSetSign();
        body_.defineAdaptationRatios(1.15, system.ReferenceResolution() / input.resolution);
        body_.defineMaterial<NeoHookeanSolid>(*input.material);
        body_.generateParticles<BaseParticles, Lattice>();
        inner_relation_ = std::make_unique<InnerRelation>(body_);
        // relax_solid(body_, *inner_relation_);

        SimpleDynamics<NormalDirectionFromBodyShape> normal_direction(body_);
        normal_direction.exec();
    };

    SolidBody &get_body() { return body_; }
    BaseParticles &get_particles() { return body_.getBaseParticles(); }
    InnerRelation &get_inner_relation() { return *inner_relation_; }
};

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for Beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

class ElasticSolidObject : public SolidObject
{
  private:
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
    ReduceDynamics<solid_dynamics::AcousticTimeStep> time_step_size_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping;
    SPH::SimpleDynamics<SPH::solid_dynamics::UpdateElasticNormalDirection> update_elastic_normal;

  public:
    ElasticSolidObject(SPHSystem &system, const solid_object_input &input)
        : SolidObject(system, input),
          corrected_configuration_(get_inner_relation()),
          stress_relaxation_first_half_(get_inner_relation()),
          stress_relaxation_second_half_(get_inner_relation()),
          time_step_size_(get_body()),
          damping(0.2, get_inner_relation(), "Velocity", get_physical_viscosity_general(input.material->getDensity(), input.material->getYoungsModulus(), input.length_scale)),
          update_elastic_normal(get_body()) {
          };
    Real get_time_step_size() { return time_step_size_.exec(); }
    void init_corrected_config() { corrected_configuration_.exec(); }
    void exec_stress_first(Real dt) { stress_relaxation_first_half_.exec(dt); }
    void exec_stress_second(Real dt) { stress_relaxation_second_half_.exec(dt); }
    void exec_damping(Real dt) { damping.exec(dt); }
    void normal_update() { update_elastic_normal.exec(); }
};

class FSIDynamics
{
  private:
    ContactRelation contact_relation_;
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration;
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal;
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid;
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<fluid_dynamics::Integration2ndHalfWithWallNoRiemann>> pressure_force_from_fluid;

  public:
    FSIDynamics(SPHBody &solid_body, RealBodyVector contact_bodies)
        : contact_relation_(solid_body, contact_bodies),
          average_velocity_and_acceleration(contact_relation_.getSPHBody()),
          insert_body_update_normal(contact_relation_.getSPHBody()),
          viscous_force_from_fluid(contact_relation_),
          pressure_force_from_fluid(contact_relation_) {}

    ContactRelation &get_contact_relation() { return contact_relation_; }
    void exec_viscous_force() { viscous_force_from_fluid.exec(); }
    void exec_pressure_force() { pressure_force_from_fluid.exec(); }
    void update_normal() { insert_body_update_normal.exec(); }
    void initialize_displacement() { average_velocity_and_acceleration.initialize_displacement_.exec(); }
    void update_averages(Real dt) { average_velocity_and_acceleration.update_averages_.exec(dt); }
};

struct solid_solid_contact
{
    SurfaceContactRelation contact_relation_;
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_density_;
    InteractionWithUpdate<solid_dynamics::ContactForce> contact_force_;

    solid_solid_contact(SPHBody &body, RealBodyVector contact_bodies)
        : contact_relation_(body, contact_bodies),
          contact_density_(contact_relation_),
          contact_force_(contact_relation_) {}
};

void aortic_valve_pulsatile_fsi(size_t res_factor)
{
    const parameters params;
    const Real resolution_leaflet = params.leaflet_thickness / 4.0;
    const Real resolution_fluid = params.stent_diameter / (25.0 * res_factor);
    const Real resolution_wall = resolution_fluid;
    const Real wall_thickness = 4 * resolution_fluid;

    // load leaflet shape
    std::cout << "leaflet loading" << std::endl;
    std::vector<SharedPtr<TriangleMeshShapeSTL>> leaflet_shape_vec;
    for (const auto &file : params.seperate_leaflet_files)
    {
        leaflet_shape_vec.emplace_back(
            makeShared<TriangleMeshShapeSTL>(file, Vec3d::Zero(), params.scale));
    }
    std::cout << "leaflet loaded" << std::endl;

    // load fluid
    BoundingBox bbox_leaflet_1 = leaflet_shape_vec[0]->getBounds();
    const Real leaflet_length = bbox_leaflet_1.upper_.y() - bbox_leaflet_1.lower_.y();
    const Real fluid_fulllength = leaflet_length + params.fluid_upstream_length + params.fluid_downstream_length;
    Vec3d fluid_translation =
        (bbox_leaflet_1.upper_.y() + params.fluid_downstream_length - 0.5 * fluid_fulllength) * Vec3d::UnitY();
    auto axis_y = SimTK::UnitVec3(0, 1, 0);

    std::cout << "fluid loading" << std::endl;
    auto fluid_shape = makeShared<ComplexShape>("blood");
    fluid_shape->add<TriangleMeshShapeCylinder>(
        axis_y, 0.5 * params.stent_diameter, 0.5 * fluid_fulllength, 10, fluid_translation);
    for (const auto &file : params.seperate_leaflet_files)
    {
        fluid_shape->subtract<TriangleMeshShapeSTL>(file, Vec3d::Zero(), params.scale);
    }
    std::cout << "fluid loaded" << std::endl;

    std::cout << "wall loading" << std::endl;
    auto wall_shape = makeShared<ComplexShape>("vessel");
    wall_shape->add<TriangleMeshShapeCylinder>(
        axis_y, 0.5 * params.stent_diameter + wall_thickness, 0.5 * fluid_fulllength, 10, fluid_translation);
    wall_shape->subtract<TriangleMeshShapeCylinder>(
        axis_y, 0.5 * params.stent_diameter, 0.5 * fluid_fulllength, 10, fluid_translation);
    std::cout << "wall loaded" << std::endl;

    std::cout << "base loading" << std::endl;
    Vec3d base_translation = 0.5 * (bbox_leaflet_1.upper_.y() + bbox_leaflet_1.lower_.y()) * Vec3d::UnitY();
    auto base_shape = makeShared<ComplexShape>("base");
    base_shape->add<TriangleMeshShapeCylinder>(
        axis_y, 0.5 * params.stent_diameter, 0.5 * leaflet_length, 10, base_translation);
    base_shape->subtract<TriangleMeshShapeCylinder>(
        axis_y, 0.5 * params.stent_diameter - params.leaflet_thickness, 0.5 * leaflet_length, 10, base_translation);
    std::cout << "base loaded" << std::endl;

    // SYSTEM
    SPHSystem system(wall_shape->getBounds(), resolution_fluid);
    IOEnvironment in_output(system);

    // SEPERATE LEAFLET OBJECTS
    std::cout << "leaflet generation starting" << std::endl;
    const Real youngs_modulus = 9.0 * params.bulk_modulus * params.shear_modulus / (3.0 * params.bulk_modulus + params.shear_modulus);
    const Real poisson_ratio = (3.0 * params.bulk_modulus - 2.0 * params.shear_modulus) / (6.0 * params.bulk_modulus + 2.0 * params.shear_modulus);
    NeoHookeanSolid material(params.rho0_s, youngs_modulus, poisson_ratio);
    std::vector<std::unique_ptr<ElasticSolidObject>> leaflet_vec;
    for (size_t i = 0; i < params.leaflet_names.size(); ++i)
    {
        solid_object_input leaflet_input{
            .name = params.leaflet_names[i],
            .mesh = leaflet_shape_vec[i],
            .material = &material,
            .resolution = resolution_leaflet,
            .length_scale = params.leaflet_thickness};
        leaflet_vec.emplace_back(std::make_unique<ElasticSolidObject>(system, leaflet_input));
        std::cout << "the following leaflet is generated: " << params.leaflet_names[i] << std::endl;
    }
    std::cout << "leaflet generation done";

    // LEAFLET FIXATION
    std::vector<std::unique_ptr<SimpleDynamics<FixBodyPartConstraint>>> fixed_constraint_vec;
    auto fixed_particle_ids = [&](ElasticSolidObject &object)
    {
        IndexVector ids;
        auto *pos = object.get_particles().getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < object.get_particles().TotalRealParticles(); ++i)
            if (base_shape->checkContain(pos[i]))
                ids.push_back(i);
        return ids;
    };
    std::vector<std::unique_ptr<BodyPartByParticle>> base_body_parts;
    for (const auto &leaflet : leaflet_vec)
    {
        base_body_parts.emplace_back(std::make_unique<BodyPartByParticle>(leaflet->get_body()));
        auto fixed_ids = fixed_particle_ids(*leaflet);
        base_body_parts.back()->body_part_particles_ = fixed_ids;
        fixed_constraint_vec.emplace_back(
            std::make_unique<SimpleDynamics<FixBodyPartConstraint>>(*base_body_parts.back()));
    }

    // FLUID
    std::cout << "fluid generation starting" << std::endl;
    FluidBody fluid_block(system, fluid_shape);
    fluid_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(params.rho0_f, params.c_f), params.mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    fluid_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);
    auto fluid_bounding_box = fluid_shape->getBounds();
    const Real fluid_y_min = fluid_bounding_box.lower_.y();
    const Real fluid_y_max = fluid_bounding_box.upper_.y();
    std::cout << "fluid generation done" << std::endl;

    // WALL - fixed part
    std::cout << "wall generation starting" << std::endl;
    solid_object_input wall_input{
        .name = "wall",
        .mesh = wall_shape,
        .material = &material,
        .resolution = resolution_wall};
    SolidObject wall_object(system, wall_input);
    std::cout << "wall generation done" << std::endl;

    //---------------------------------------------------------------------------//
    // ALGORITHMS
    //---------------------------------------------------------------------------//
    // FLUID
    InnerRelation fluid_inner(fluid_block);
    ContactRelation fluid_contact(fluid_block, {&wall_object.get_body(), &leaflet_vec[0]->get_body(),
                                                &leaflet_vec[1]->get_body(),
                                                &leaflet_vec[2]->get_body()});
    ComplexRelation fluid_block_complex(fluid_inner, fluid_contact);
    InteractionDynamics<NablaWVComplex> kernel_summation(fluid_inner, fluid_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(fluid_inner, fluid_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(fluid_inner, fluid_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(fluid_inner, fluid_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(fluid_inner, fluid_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(fluid_inner, fluid_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(fluid_block, params.U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(fluid_block);

    //-----------------------------------Inlet and outlet
    const Real buffer_length = 3.0 * resolution_fluid;
    const Vec3d bidirectional_buffer_halfsize(buffer_length * 0.5,
                                              params.stent_diameter * 0.505,
                                              params.stent_diameter * 0.505);

    auto compute_rotation = [](const Vec3d &normal) -> Rotation3d
    {
        Vec3d source = Vec3d::UnitX();
        Vec3d target = normal.normalized();
        Eigen::Quaternion<Real> quaternion =
            Eigen::Quaternion<Real>::FromTwoVectors(source, target);
        return Rotation3d(quaternion);
    };

    const Vec3d inlet_normal = Vec3d::UnitY();
    const Vec3d inlet_pos = fluid_y_min * Vec3d::UnitY();
    const Vec3d inlet_center = inlet_pos + 0.5 * buffer_length * inlet_normal;
    AlignedBoxByCell left_emitter(fluid_block, AlignedBox(xAxis, Transform(compute_rotation(inlet_normal), inlet_center), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);

    const Vec3d outlet_normal = -Vec3d::UnitY();
    const Vec3d outlet_pos = fluid_y_max * Vec3d::UnitY();
    const Vec3d outlet_center = outlet_pos + 0.5 * buffer_length * outlet_normal;
    AlignedBoxByCell right_emitter(fluid_block, AlignedBox(xAxis, Transform(compute_rotation(outlet_normal), outlet_center), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(fluid_inner, fluid_contact);

    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //---------------------------------------------------------------------------//
    // FSI
    std::vector<std::unique_ptr<FSIDynamics>> fsi_algorithm_vec;
    fsi_algorithm_vec.reserve(leaflet_vec.size());
    for (auto &leaflet : leaflet_vec)
        fsi_algorithm_vec.emplace_back(
            std::make_unique<FSIDynamics>(leaflet->get_body(), RealBodyVector{&fluid_block}));
    //-------------------Contact between leaflets-----------------------------//
    std::vector<std::unique_ptr<solid_solid_contact>> leaflet_contact_vec;
    leaflet_contact_vec.emplace_back(std::make_unique<solid_solid_contact>(
        leaflet_vec[0]->get_body(), RealBodyVector{&leaflet_vec[1]->get_body(), &leaflet_vec[2]->get_body()}));
    leaflet_contact_vec.emplace_back(std::make_unique<solid_solid_contact>(
        leaflet_vec[1]->get_body(), RealBodyVector{&leaflet_vec[0]->get_body(), &leaflet_vec[2]->get_body()}));
    leaflet_contact_vec.emplace_back(std::make_unique<solid_solid_contact>(
        leaflet_vec[2]->get_body(), RealBodyVector{&leaflet_vec[0]->get_body(), &leaflet_vec[1]->get_body()}));
    //---------------------------------------------------------------------------//
    // OUTPUT
    BodyStatesRecordingToVtp vtp_binary_output(system);
    vtp_binary_output.addToWrite<Real>(fluid_block, "Pressure");
    vtp_binary_output.addToWrite<Vecd>(fluid_block, "Velocity");
    vtp_binary_output.addToWrite<Real>(fluid_block, "Density");
    vtp_binary_output.addToWrite<int>(fluid_block, "Indicator");
    vtp_binary_output.addToWrite<Vec3d>(fluid_block, "KernelSummation");
    vtp_binary_output.addToWrite<int>(fluid_block, "BufferIndicator");
    for (const auto &leaflet : leaflet_vec)
    {
        vtp_binary_output.addToWrite<Vecd>(leaflet->get_body(), "Velocity");
        vtp_binary_output.addToWrite<Vecd>(leaflet->get_body(), "NormalDirection");
    }

    //----------------------------------------------------------------------
    // INITIALIZATION
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    for (const auto &leaflet : leaflet_vec)
        leaflet->init_corrected_config();
    //-------------initial relaxation for fluid------------------
    auto run_fluid_relaxation = [&]()
    {
        size_t relaxation_fluid_itr = 0;
        std::cout << "Fluid relaxation starting..." << std::endl;
        while (relaxation_fluid_itr < 100)
        {
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();

            Real Dt = get_fluid_advection_time_step_size.exec();
            transport_velocity_correction.exec(Dt);

            relaxation_fluid_itr++;

            fluid_block.updateCellLinkedList();
            fluid_block_complex.updateConfiguration();
        }
        std::cout << "Fluid relaxation finished." << std::endl;
    };
    run_fluid_relaxation();
    vtp_binary_output.writeToFile(0);

    //----------------------------------------------------------------------
    //	Main loop of FSI starts here.
    //----------------------------------------------------------------------
    const Real Dt_ref = get_fluid_advection_time_step_size.exec();
    const Real dt_ref = get_fluid_time_step_size.exec();
    const Real dt_s_ref = leaflet_vec[0]->get_time_step_size();
    std::cout << "Reference time steps: "
              << "Dt_ref = " << Dt_ref
              << ", dt_ref = " << dt_ref
              << ", dt_s_ref = " << dt_s_ref
              << std::endl;

    auto *F = leaflet_vec[0]->get_particles().getVariableDataByName<Matd>("DeformationGradient");
    auto &F0 = F[0];
    Real Q[9], H[9], A[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i * 3 + j] = F0(i, j);
    polar::polar_decomposition(Q, H, A);
    std::cout << "Q: " << Q << std::endl;
    exit(0);

    auto run_fsi = [&]()
    {
        std::cout << "-----------FSI-----------" << std::endl;
        // Time stepping
        size_t number_of_iterations = 0;
        int screen_output_interval = 1;
        const int max_num_outputs = 200;
        const Real D_Time = params.end_time / max_num_outputs; /**< time stamps for output. */
        Real Dt = 0.0;
        Real dt = 0.0;
        Real dt_s = 0.0;
        Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
        while (physical_time < params.end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < D_Time)
            {
                Dt = get_fluid_advection_time_step_size.exec();
                std::cout << "Dt: " << Dt << std::endl;
                if (std::isnan(Dt))
                {
                    throw std::runtime_error("Dt is NaN!");
                }
                if (Dt < Dt_ref / 20.0)
                    throw std::runtime_error("Dt is too small!");

                boundary_indicator.exec();
                update_fluid_density.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec(Dt);

                std::cout << "line 587" << std::endl;

                // FSI
                for (auto &alg : fsi_algorithm_vec)
                    alg->exec_viscous_force();

                std::cout << "line 593" << std::endl;

                for (auto &leaflet : leaflet_vec)
                    leaflet->normal_update();

                std::cout << "line 598" << std::endl;

                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    dt = get_fluid_time_step_size.exec();
                    if (std::isnan(dt))
                    {
                        throw std::runtime_error("dt is NaN!");
                    }
                    if (dt < dt_ref / 20.0)
                        throw std::runtime_error("dt is too small!");
                    dt = SMIN(dt, Dt - relaxation_time);
                    std::cout << "dt: " << dt << std::endl;

                    // fluid
                    pressure_relaxation.exec(dt);
                    kernel_summation.exec();
                    left_inflow_pressure_condition.exec(dt);
                    right_inflow_pressure_condition.exec(dt);
                    inflow_velocity_condition.exec();
                    // fsi
                    for (auto &alg : fsi_algorithm_vec)
                        alg->exec_pressure_force();
                    // fluid
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    inner_ite_dt_s = 0;
                    Real dt_s_sum = 0.0;
                    for (auto &alg : fsi_algorithm_vec)
                        alg->initialize_displacement();
                    while (dt_s_sum < dt)
                    {
                        { // time step
                            StdVec<Real> dt_vec;
                            dt_vec.reserve(leaflet_vec.size());
                            for (auto &leaflet : leaflet_vec)
                                dt_vec.push_back(leaflet->get_time_step_size());
                            dt_vec.push_back(dt - dt_s_sum);
                            dt_s = *std::min_element(dt_vec.begin(), dt_vec.end());
                        }
                        std::cout << "dt_s: " << dt_s << std::endl;
                        // normal update
                        for (auto &leaflet : leaflet_vec)
                            leaflet->normal_update();
                        // exec solid-solid contacts
                        for (auto &contact : leaflet_contact_vec)
                            contact->contact_density_.exec();
                        for (auto &contact : leaflet_contact_vec)
                            contact->contact_force_.exec();

                        // leaflet
                        for (auto &leaflet : leaflet_vec)
                            leaflet->exec_stress_first(dt_s);
                        for (auto &constraint : fixed_constraint_vec)
                            constraint->exec();
                        for (auto &leaflet : leaflet_vec)
                            leaflet->exec_damping(dt_s);
                        for (auto &constraint : fixed_constraint_vec)
                            constraint->exec();
                        for (auto &leaflet : leaflet_vec)
                        {
                            leaflet->exec_stress_second(dt_s);
                            leaflet->get_body().updateCellLinkedList();
                        }

                        for (auto &contact : leaflet_contact_vec)
                            contact->contact_relation_.updateConfiguration();

                        dt_s_sum += dt_s;
                        ++inner_ite_dt_s;
                    }
                    for (auto &alg : fsi_algorithm_vec)
                        alg->update_averages(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                    ++inner_ite_dt;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << "N=" << number_of_iterations
                              << "\tTime = " << physical_time
                              << "\tDt = " << Dt
                              << "\tdt = " << dt
                              << "\tdt_s = " << dt_s
                              << std::endl;
                }
                ++number_of_iterations;

                // Fluid update / Config update
                left_bidirection_buffer.injection.exec();
                right_bidirection_buffer.injection.exec();
                // then do deletion for all buffers
                left_bidirection_buffer.deletion.exec();
                right_bidirection_buffer.deletion.exec();

                fluid_block.updateCellLinkedList();
                fluid_block_complex.updateConfiguration();
                for (const auto &alg : fsi_algorithm_vec)
                    alg->get_contact_relation().updateConfiguration();
                boundary_indicator.exec();
                left_bidirection_buffer.tag_buffer_particles.exec();
                right_bidirection_buffer.tag_buffer_particles.exec();
            }

            /** write run-time observation into file */
            vtp_binary_output.writeToFile();
        }
        std::cout << "-----------FSI FINISHED-----------" << std::endl;
    };

    TickCount t1 = TickCount::now();
    try
    {
        run_fsi();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        for (const auto &leaflet : leaflet_vec)
            leaflet->get_body().setNewlyUpdated();
        fluid_block.setNewlyUpdated();
        vtp_binary_output.writeToFile();
    }
    TickCount t2 = TickCount::now();
    TimeInterval tt = t2 - t1;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds" << std::endl;
}