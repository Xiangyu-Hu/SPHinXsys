/**
 * @file    water_entry_exit_3d.cpp
 * @brief   3D water entry and exit example with surface wetting considered.
 * @details This is a 3D FSI test case using M1.stl for spatial temporal identification approach.
 * @author  Siyu Zou Shuoguo Zhang and Xiangyu Hu
 */
#include "SimTKsimbody.h"
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.5; /**< Water tank length. */
Real DW = 1.0; /**< Water tank width. */
Real DH = 1.0; /**< Water tank height. */
Real LL = 1.5; /**< Water column length. */
Real LW = 1; /**< Water column width. */
Real LH = 0.5; /**< Water column height. */

Real particle_spacing_ref = 0.02;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */

//Vecd object_center(0.5 * DL, LH + 0.2, 0.5 * DW); /**< Location of the object center. */
Vecd object_center(0.46, 0.528, 0.5); /**< Location of the object center. */
// 初始速度参数
Real initial_speed = 70.0;                     /**< Initial velocity magnitude (m/s). */
Real initial_angle = -20.0 * Pi / 180.0;       /**< Initial velocity angle (radians, negative = downward). */
Real initial_pitch_angle = -20.0 * Pi / 180.0; /**< Initial pitch angle about Z-axis(radians). */

//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                         /**< Fluid density (water). */
Real rho0_s = 2063.0;                         /**< Object density. */
Real gravity_g = 9.81;                        /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH) * 10; /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                      /**< Reference sound speed. */
Real mu_f = 8.9e-4;                           /**< Water dynamics viscosity (Pa·s). */

//----------------------------------------------------------------------
//  Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 100.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real object_moisture = 0.0;                                  /**< object moisture. */
Real wall_moisture = 1.0;                                    /**< wall moisture. */

//----------------------------------------------------------------------
//  计算初始速度分量
//----------------------------------------------------------------------
Vecd calculateInitialVelocity()
{
    // 在XZ平面内的初始速度（Y方向为垂直方向）
    Real horizontal_velocity = initial_speed * cos(initial_angle);
    Real vertical_velocity = initial_speed * sin(initial_angle);

    return Vecd(horizontal_velocity, vertical_velocity, 0.0);
}

//----------------------------------------------------------------------
//  全局力到局部坐标系的转换（3D版本）
//----------------------------------------------------------------------
Vec3d transformGlobalForceToLocal(const Vec3d &global_force, Real rotation_angle_z)
{

    Real cos_theta = cos(rotation_angle_z);
    Real sin_theta = sin(rotation_angle_z);

    Vec3d local_force;
    local_force[0] = global_force[0] * cos_theta + global_force[1] * sin_theta;  // X方向
    local_force[1] = -global_force[0] * sin_theta + global_force[1] * cos_theta; // Y方向
    local_force[2] = global_force[2];

    return local_force;
}

//----------------------------------------------------------------------
//   观测点定义
//----------------------------------------------------------------------
StdVec<Vecd> observer_location = {object_center}; /**< Displacement observation point. */
StdVec<Vecd> wetting_observer_location =
    {object_center - Vecd(0.0, 0.05, 0.0)}; /**< wetting observation point. */
// 前缘观测点（入水点）
StdVec<Vecd> leading_edge_location = {object_center}; /**< Leading edge observation point. */

//----------------------------------------------------------------------
//  汇总输出类（3D版本）
//----------------------------------------------------------------------
class SummaryOutput3D
{
  public:
    SummaryOutput3D(const std::string &filename)
        : output_file_(filename)
    {
        output_file_ << "Time[s]   "
                     << "LeadingEdge_X LeadingEdge_Y LeadingEdge_Z   "
                     << "ViscousForceGlobal_X ViscousForceGlobal_Y ViscousForceGlobal_Z   "
                     << "PressureForceGlobal_X PressureForceGlobal_Y PressureForceGlobal_Z   "
                     << "ViscousForceLocal_X ViscousForceLocal_Y ViscousForceLocal_Z   "
                     << "PressureForceLocal_X PressureForceLocal_Y PressureForceLocal_Z   "
                     << "TotalForceLocal_X TotalForceLocal_Y TotalForceLocal_Z\n"; // 新增
    }

    void writeData(Real time,
                   const Vecd &leading_edge_position,
                   const Vecd &viscous_force_global,
                   const Vecd &pressure_force_global,
                   const Vecd &viscous_force_local,
                   const Vecd &pressure_force_local,
                   const Vecd &total_force_local) // 新增参数
    {
        output_file_ << std::scientific << std::setprecision(9)
                     << time << " "
                     << leading_edge_position[0] << " " << leading_edge_position[1] << " " << leading_edge_position[2] << " "
                     << viscous_force_global[0] << " " << viscous_force_global[1] << " " << viscous_force_global[2] << " "
                     << pressure_force_global[0] << " " << pressure_force_global[1] << " " << pressure_force_global[2] << " "
                     << viscous_force_local[0] << " " << viscous_force_local[1] << " " << viscous_force_local[2] << " "
                     << pressure_force_local[0] << " " << pressure_force_local[1] << " " << pressure_force_local[2] << " "
                     << total_force_local[0] << " " << total_force_local[1] << " " << total_force_local[2] << "\n"; // 新增
        output_file_.flush();
    }

  private:
    std::ofstream output_file_;
};
//----------------------------------------------------------------------
//  Definition for water body (3D)
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        // 创建水体长方体
        add<GeometricShapeBox>(Transform(Vecd(LL / 2, LH / 2, LW / 2)), Vecd(LL / 2, LH / 2, LW / 2));
    }
};

//----------------------------------------------------------------------
//  Definition for wall body (3D)
//----------------------------------------------------------------------
 class WallBoundary : public ComplexShape {
public:
explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
{
    Real eps = BW; // 微小重叠量，防止数值缝隙
    add<GeometricShapeBox>(Transform(Vecd(DL / 2, -BW / 2, DW / 2)), Vecd(DL / 2 + eps, BW / 2, DW / 2 + eps)); // 底部平板 (y从 -BW 到 0)
    add<GeometricShapeBox>(Transform(Vecd(-BW / 2, DH / 2, DW / 2)), Vecd(BW / 2, DH / 2 + eps, DW / 2 + eps)); // 左侧平板 (x从 -BW 到 0)
    add<GeometricShapeBox>(Transform(Vecd(DL + BW / 2, DH / 2, DW / 2)), Vecd(BW / 2, DH / 2 + eps, DW / 2 + eps)); // 右侧平板 (x从 DL 到 DL+BW)
    add<GeometricShapeBox>(Transform(Vecd(DL / 2, DH / 2, -BW / 2)), Vecd(DL / 2 + eps, DH / 2 + eps, BW / 2)); // 前侧平板 (z从 -BW 到 0)
    add<GeometricShapeBox>(Transform(Vecd(DL / 2, DH / 2, DW + BW / 2)), Vecd(DL / 2 + eps, DH / 2 + eps, BW / 2)); // 后侧平板 (z从 DW 到 DW+BW)
    // 顶部平板 (y从 DH 到 DH+BW)，若需要封闭可取消注释
    // add<GeometricShapeBox>(Transform(Vecd(DL/2, DH + BW/2, DW/2)),Vecd(DL/2 + eps, BW/2, DW/2 + eps));
}};
//----------------------------------------------------------------------
// 运行粒子松弛（3D版本）
//----------------------------------------------------------------------
int runParticleRelaxation(SolidBody &floating_object)
{
    std::cout << "Starting 3D particle relaxation..." << std::endl;

    InnerRelation object_inner(floating_object); // 创建内部关系

    // 松弛方法
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(floating_object);
    BodyStatesRecordingToVtp write_inserted_body_to_vtp(floating_object); // 写入Vtp文件
    ReloadParticleIO write_particle_reload_files(floating_object); // 写入粒子重启文件
    RelaxationStepInner relaxation_step_inner(object_inner);// 物理松弛步骤（3D版本）

    //----------------------------------------------------------------------
    // 粒子松弛开始
    //----------------------------------------------------------------------
    random_inserted_body_particles.exec(0.25); // 随机初始化粒子位置
    relaxation_step_inner.SurfaceBounding().exec(); // 表面约束（确保粒子在表面内）
    write_inserted_body_to_vtp.writeToFile(0); // 输出初始状态

    //----------------------------------------------------------------------
    // 松弛浮体粒子
    //----------------------------------------------------------------------
    int ite_p = 0;
    int max_iterations = 5000; // 3D可能需要更多迭代

    while (ite_p < max_iterations)
    {
        relaxation_step_inner.exec();
        ite_p += 1;

        if (ite_p % 200 == 0)
        {
            std::cout << std::fixed << std::setprecision(9)
                      << "Relaxation steps for the 3D object N = " << ite_p << "\n";
            write_inserted_body_to_vtp.writeToFile(ite_p);
        }

        // 3D中可能需要额外的检查
        if (ite_p % 500 == 0)
        {
            // 更新邻居列表
            floating_object.updateCellLinkedList();
            object_inner.updateConfiguration();
        }
    }

    std::cout << "3D physics relaxation process finished!" << std::endl;

    // 输出结果
    write_particle_reload_files.writeToFile(0);

    return 0;
}

//----------------------------------------------------------------------
//  Definition for object body using STL 
//----------------------------------------------------------------------
class FloatingObject : public ComplexShape
{
  public:
    explicit FloatingObject(const std::string &shape_name) : ComplexShape(shape_name)
    {
        std::string full_path_to_file = "./input/M1TR3.stl";
        Vecd translation(0.0, 0.0, 0.0);

        add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
    }
};
//----------------------------------------------------------------------
// 创建Simbody约束区域的形状
//----------------------------------------------------------------------
std::shared_ptr<ComplexShape> createSimbodyConstrainShape3D()
{
    auto constrain_shape = makeShared<ComplexShape>("SimbodyConstrainShape");

    std::string full_path_to_file = "./input/M1TR3.stl";

    Vecd translation(0.0, 0.0, 0.0);

    constrain_shape->add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);

    return constrain_shape;
}

//----------------------------------------------------------------------
//  湿表面初始化模板类
//----------------------------------------------------------------------
template <class BodyType>
class WettingBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingBodyInitialCondition(BodyType &sph_body, Real moisture)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)),
          moisture_(moisture) {}

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = moisture_;
    }
  protected:
    Real *phi_;
    Real moisture_;
};

//----------------------------------------------------------------------
// 获取物体的旋转角度
//----------------------------------------------------------------------
Real getObjectRotationAngle(SimTK::MobilizedBody &mob, const SimTK::State &state)
{
    SimTK::Rotation rot = mob.getBodyRotation(state);
    SimTK::Vec3 angles = rot.convertRotationToBodyFixedXYZ();
    // 在XZ平面运动，我们关注绕Y轴的旋转（pitch）
    // 根据欧拉角顺序，angles[1]对应绕Y轴的旋转
    // 2 z轴旋转
    return angles[2];
}

//----------------------------------------------------------------------
// 物体初始速度设置类（3D完整版本）
//----------------------------------------------------------------------
class ObjectInitialVelocity : public LocalDynamics
{
  public:
    explicit ObjectInitialVelocity(SPHBody &sph_body,
                                   const Vec3d &initial_vel,
                                   const Vec3d &initial_angular_vel = Vec3d(0.0, 0.0, 0.0))
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vec3d>("Position")),
          vel_(particles_->getVariableDataByName<Vec3d>("Velocity")),
          initial_vel_(initial_vel),
          initial_angular_vel_(initial_angular_vel),
          mass_(particles_->getVariableDataByName<Real>("Mass"))
    {
        calculateMassCentroid();
    }

    void calculateMassCentroid()
    {
        centroid_ = Vec3d(0.0, 0.0, 0.0);
        total_mass_ = 0.0;
        size_t total_particles = particles_->TotalRealParticles();

        for (size_t i = 0; i < total_particles; ++i)
        {
            Real mass_i = mass_[i];
            centroid_ += mass_i * pos_[i];
            total_mass_ += mass_i;
        }
        if (total_mass_ > 0)
            centroid_ /= total_mass_;
    }

    void update(size_t index_i, Real dt)
    {
        Vec3d r = pos_[index_i] - centroid_;

        // 设置平动速度
        vel_[index_i] = initial_vel_;

        // 添加绕质心的旋转速度（3D版本）
        if (initial_angular_vel_.norm() > 1e-9)
        {
            // 角速度 × 位置向量 = 切向速度
            Vec3d tangential_vel = initial_angular_vel_.cross(r);
            vel_[index_i] += tangential_vel;
        }
    }

  protected:
    Vec3d *pos_;
    Vec3d *vel_;
    Real *mass_;
    Vec3d initial_vel_;
    Vec3d initial_angular_vel_;
    Vec3d centroid_;
    Real total_mass_;
};

//// 为不同类型创建别名
using WettingFluidBodyInitialCondition = WettingBodyInitialCondition<FluidBody>;
using WettingWallBodyInitialCondition = WettingBodyInitialCondition<SolidBody>;
using WettingObjectBodyInitialCondition = WettingBodyInitialCondition<SolidBody>;

//----------------------------------------------------------------------
//  Main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{

    //----------------------------------------------------------------------
    //  Build up an SPHSystem
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);

    // 运行配置
     sph_system.setRunParticleRelaxation(false);
     sph_system.setReloadParticles(true);
    //sph_system.setRunParticleRelaxation(true);
    //sph_system.setReloadParticles(false);

    sph_system.handleCommandlineOptions(ac, av);

    if (sph_system.RunParticleRelaxation())
    {
        // 只创建浮体进行松弛
        SolidBody floating_object(sph_system, makeShared<FloatingObject>("FloatingObject"));
        floating_object.defineAdaptationRatios(1.15, 1.0);
        floating_object.defineBodyLevelSetShape();
        floating_object.defineClosure<Solid, IsotropicDiffusion>(
            rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
        floating_object.generateParticles<BaseParticles, Lattice>();

        // 松弛逻辑
        return runParticleRelaxation(floating_object);
    }

    //----------------------------------------------------------------------
    //  Creating bodies with corresponding materials and particles
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(
        ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody floating_object(sph_system, makeShared<FloatingObject>("FloatingObject"));
    floating_object.defineAdaptationRatios(1.15, 1.0);
    floating_object.defineBodyLevelSetShape();

    // 定义材料，包含湿表面扩散
    floating_object.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));

    // 根据系统设置选择粒子生成方式
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        floating_object.generateParticles<BaseParticles, Reload>(floating_object.getName());
    }
    else
    {
        floating_object.generateParticles<BaseParticles, Lattice>();
    }

    // 注册额外状态变量
    BaseParticles &object_particles = floating_object.getBaseParticles();
    object_particles.registerStateVariableData<Vecd>("ViscousForceFromFluid");
    object_particles.registerStateVariableData<Vecd>("PressureForceFromFluid");
    // 确保Velocity也注册
    object_particles.registerStateVariableData<Vecd>("Velocity");

    auto viscous_force_data = object_particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    auto pressure_force_data = object_particles.getVariableDataByName<Vecd>("PressureForceFromFluid");

    // 观测体
    ObserverBody object_observer(sph_system, "ObjectObserver");
    object_observer.generateParticles<ObserverParticles>(observer_location);

    ObserverBody wetting_observer(sph_system, "WettingObserver");
    wetting_observer.generateParticles<ObserverParticles>(wetting_observer_location);

    ObserverBody leading_edge_observer(sph_system, "LeadingEdgeObserver");
    leading_edge_observer.generateParticles<ObserverParticles>(leading_edge_location);

    //----------------------------------------------------------------------
    //  Define body relation map
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation object_inner(floating_object);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &floating_object});
    ContactRelation object_contact(floating_object, {&water_block});
    ContactRelation object_observer_contact(object_observer, {&floating_object});
    ContactRelation wetting_observer_contact(wetting_observer, {&floating_object});
    ContactRelation leading_edge_observer_contact(leading_edge_observer, {&floating_object});

    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    //  Define the fluid dynamics
    //----------------------------------------------------------------------
    // 湿表面扩散模型
    GetDiffusionTimeStepSize get_diffusion_time_step(floating_object);

    // 湿表面扩散松弛
    using ObjectFluidDiffusionDirichlet =
        DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;
    ObjectFluidDiffusionDirichlet object_wetting(object_contact);

    // 初始化湿表面条件
    SimpleDynamics<WettingBodyInitialCondition<FluidBody>>
        wetting_water_initial_condition(water_block, fluid_moisture);
    SimpleDynamics<WettingBodyInitialCondition<SolidBody>>
        wetting_wall_initial_condition(wall_boundary, wall_moisture);
    SimpleDynamics<WettingBodyInitialCondition<SolidBody>>
        wetting_object_initial_condition(floating_object, object_moisture);

    // 设置物体初始速度
    Vecd initial_velocity = calculateInitialVelocity();
    SimpleDynamics<ObjectInitialVelocity> object_set_initial_velocity(floating_object, initial_velocity);

    // 重力
    Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);

    // 自由表面识别（使用带湿表面耦合的版本）
    InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex>
        free_stream_surface_indicator(water_block_inner, water_block_contact);

    // 法向计算
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> object_normal_direction(floating_object);

    // 流体动力学
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>transport_velocity_correction(water_block_inner, water_block_contact);


    // 时间步控制
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep>fluid_acoustic_time_step(water_block);

    //----------------------------------------------------------------------
    //  Algorithms of FSI
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(object_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>>pressure_force_from_fluid(object_contact);

    // 力统计（全局坐标系）
    ReducedQuantityRecording<QuantitySummation<Vecd>>write_total_viscous_force(floating_object, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>>write_total_pressure_force(floating_object, "PressureForceFromFluid");

    // 前缘点位置插值
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_leading_edge_position(leading_edge_observer_contact, "Position", "Position");
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //  Building Simbody (3D)
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);

    // 创建约束区域（使用与物体相同的形状）
    SolidBodyPartForSimbody structure_system(floating_object,
                                             createSimbodyConstrainShape3D());

    SimTK::Body::Rigid fixed_spot_info(*structure_system.body_part_mass_properties_);


    // 使用正确的系留点
    //Vecd tethering_point = object_center; // 使用物体中心作为系留点
    Vecd actual_centroid = structure_system.initial_mass_center_;
    Vecd tethering_point = actual_centroid; // 使用实际质心
    SimTK::MobilizedBody::Planar structure_mob(matter.Ground(),
                                               SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], tethering_point[2])),
                                               fixed_spot_info, SimTK::Transform(SimTKVec3(0)));

    // 质量属性
    SimTK::Body::Rigid structure_info(*structure_system.body_part_mass_properties_);

    // 重力
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(0.0, 0.0, -gravity_g), 0.0);

    // 离散力
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);

    // 初始化状态
    SimTK::State state = MBsystem.realizeTopology();


    SimTK::Vec3 mobilizer_vel(0.0, initial_velocity[0], initial_velocity[1]);
    structure_mob.setU(state, mobilizer_vel);

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);


    // 输出Simbody状态
    std::cout << "Debug: Simbody initialized successfully" << std::endl;
    std::cout << "Debug: Joint U = ("
              << structure_mob.getU(state)[0] << ", "
              << structure_mob.getU(state)[1] << ", "
              << structure_mob.getU(state)[2] << ")" << std::endl;
    std::cout << "Mass: " << structure_system.body_part_mass_properties_->getMass() << std::endl;
    std::cout << "Center of mass: " << structure_system.initial_mass_center_.transpose() << std::endl;

    //----------------------------------------------------------------------
    //  Coupling between SimBody and SPH
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_structure(structure_system, MBsystem, structure_mob, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_on_structure(structure_system, MBsystem, structure_mob, integ);

    //----------------------------------------------------------------------
    //  I/O operations and observations
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(floating_object, diffusion_species_name);

    ObservedQuantityRecording<Vecd> write_object_displacement("Position", object_observer_contact);
    ObservedQuantityRecording<Real> write_object_wetting(diffusion_species_name, wetting_observer_contact);
    ObservedQuantityRecording<Vecd> write_leading_edge_position("Position", leading_edge_observer_contact);

    SummaryOutput3D summary_output("./output/SummaryOutput3D.dat"); // 创建汇总输出

    //----------------------------------------------------------------------
    //  Prepare the simulation
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    wall_boundary_normal_direction.exec();
    object_normal_direction.exec();

    wetting_water_initial_condition.exec();
    wetting_wall_initial_condition.exec();
    wetting_object_initial_condition.exec();

    object_set_initial_velocity.exec();
    Real dt_thermal = get_diffusion_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity.exec();

    //----------------------------------------------------------------------
    //  Time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1;
    int observation_sample_interval = screen_output_interval * 1;

    Real end_time = 0.015; // 模拟结束时间
    Real output_interval = 0.0005;

    /*Real output_interval = end_time / 100.0;*/

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //  First output
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_object_displacement.writeToFile(number_of_iterations);
    write_object_wetting.writeToFile(number_of_iterations);
    write_leading_edge_position.writeToFile(number_of_iterations);

    std::cout << "Initialization complete. Starting simulation..." << std::endl;
    std::cout << "Initial velocity: " << initial_velocity.transpose() << std::endl;
    std::cout << "Object center: " << object_center.transpose() << std::endl;

    //----------------------------------------------------------------------
    //  Main loop
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;

        while (integration_time < output_interval)
        {

            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();
            fluid_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;

            while (relaxation_time < Dt)
            {
                viscous_force_from_fluid.exec();  // 计算流体对圆柱的粘性力（全局）
                pressure_force_from_fluid.exec(); // 计算流体对圆柱的压力力（全局）

                // 获取时间步长（考虑流体声学时间步和扩散时间步）
                Real dt_fluid = fluid_acoustic_time_step.exec();
                Real dt_diff = get_diffusion_time_step.exec();

                dt = SMIN(SMIN(dt_fluid, dt_diff), Dt - relaxation_time);
                // 防止dt为0或负值
                if (dt < 1e-12)
                {
                    std::cerr << "Warning: dt too small: " << dt << std::endl;
                    dt = 1e-7;
                }

                fluid_pressure_relaxation.exec(dt);
                // 在压力松弛后重新计算力，用于下一步的Simbody耦合
                 pressure_force_from_fluid.exec();
                 viscous_force_from_fluid.exec();
                fluid_density_relaxation.exec(dt);
  
                object_wetting.exec(dt);// 湿表面扩散

                 // 耦合多体动力学
                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, structure_mob,
                                                force_on_structure.exec());
                // --- 再执行约束更新，把刚体新位姿同步回 SPH 颗粒（在步进之后）
                constraint_on_structure.exec();

                // 更新前缘点位置（或其它需基于新位姿的 SPH 操作）
                interpolation_leading_edge_position.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;
            // 输出和记录
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(6)
                          << "N=" << number_of_iterations
                          << " Time = " << physical_time
                          << " Dt = " << Dt << " dt = " << dt << "\n"
                          << " integration_time = " << integration_time
                          << " output_interval = " << output_interval << "\n";

                // 执行插值以更新前缘观测点的值
                interpolation_leading_edge_position.exec();

                // 获取前缘观测点数据
                auto &leading_edge_particles = leading_edge_observer.getBaseParticles();

                Vecd leading_edge_pos = leading_edge_particles.getVariableDataByName<Vecd>("Position")[0];
                //SimTK::State &output_state = integ.updAdvancedState();
                // 获取全局力
                write_total_viscous_force.writeToFile(number_of_iterations);
                write_total_pressure_force.writeToFile(number_of_iterations);

                // 计算粒子上的总力（全局坐标系）
                BaseParticles &object_particles = floating_object.getBaseParticles();
                Vecd total_viscous_force_global(0.0, 0.0, 0.0);
                Vecd total_pressure_force_global(0.0, 0.0, 0.0);

                auto viscous_force_data = object_particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
                auto pressure_force_data = object_particles.getVariableDataByName<Vecd>("PressureForceFromFluid");

                size_t total_real_particles = object_particles.TotalRealParticles();
                for (size_t i = 0; i < total_real_particles; ++i)
                {
                    total_viscous_force_global += viscous_force_data[i];
                    total_pressure_force_global += pressure_force_data[i];
                }


                SimTK::State current_state = integ.getAdvancedState(); // 获取当前物体的旋转角度（从Simbody状态）


                MBsystem.realize(current_state, SimTK::Stage::Velocity); // 确保状态已实现

                //Real rotation_angle_z = getObjectRotationAngle(structure_mob, current_state);
                Real rotation_angle_z = initial_pitch_angle+getObjectRotationAngle(structure_mob, current_state);

                // 将全局力转换到局部坐标系
                Vecd viscous_force_local = transformGlobalForceToLocal(total_viscous_force_global, rotation_angle_z);
                Vecd pressure_force_local = transformGlobalForceToLocal(total_pressure_force_global, rotation_angle_z);



                Vecd total_force_local = viscous_force_local + pressure_force_local;

                summary_output.writeData(physical_time,
                                         leading_edge_pos,
                                         total_viscous_force_global,
                                         total_pressure_force_global,
                                         viscous_force_local,
                                         pressure_force_local,
                                         total_force_local); // 新增参数

                if (number_of_iterations % observation_sample_interval == 0)
                {
                    write_object_displacement.writeToFile(number_of_iterations);
                    write_object_wetting.writeToFile(number_of_iterations);
                    write_leading_edge_position.writeToFile(number_of_iterations);
                }
            }

            number_of_iterations++;
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }

            // 更新邻居列表和配置

                water_block.updateCellLinkedList();
                wall_boundary.updateCellLinkedList();
                floating_object.updateCellLinkedList();

                water_block_inner.updateConfiguration();
                object_inner.updateConfiguration();
                object_contact.updateConfiguration();
                water_block_complex.updateConfiguration();
                leading_edge_observer_contact.updateConfiguration();

                free_stream_surface_indicator.exec();
                interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}