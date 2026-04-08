/**
 * @file	water entry and exit.cpp
 * @brief	2D water entry and exit example with surface wetting considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach,
 *          especially when coupled with the wetting.
 * @author  Shuoguo Zhang and Xiangyu Hu
 */
//#include "2d_flow_around_cylinder.h"
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//	核心物理参数（实际工程值）
//----------------------------------------------------------------------
// 几何参数
Real actual_cylinder_length = 0.4796; /**< 圆柱总长 (m) */
Real actual_cylinder_radius = 0.025;  /**< 圆柱半径 (m) */
// 质量与转动惯量（kg·m²）
Real actual_total_mass = 1.855; /**< 圆柱总质量 (kg) */
Real Ix = 0.0009;               /**< 绕圆柱轴线（x轴）转动惯量 */
Real Iy = 0.020;                /**< 绕y轴转动惯量 */
Real Iz = 0.020;                /**< 绕z轴转动惯量（2D模拟主惯量） */
// 缩放因子：原raw_points总长度=4.796 → 实际0.4796m → 缩放因子=0.1
Real scale_factor = actual_cylinder_length / 4.796; 
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4;                                                /**< Water tank length. */
Real DH = 5;                                                /**< Water tank height. */
Real LH = 2;                                             /**< Water column height. */
Real particle_spacing_ref = 0.0025; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                       /**< Thickness of tank wall. */

Vec2d cylinder_center(0.05 * DL, LH + 0.2);               /**< Location of the cylinder center. */

// 初始速度参数
Real initial_speed =70;                         /**< Initial velocity magnitude (m/s). */
Real initial_angle =-20 * Pi / 180.0;          /**< Initial velocity angle (radians, negative = downward). */
Real initial_rotation_angle = -20 * Pi / 180.0; /**< Initial body rotation (radians). */
Real initial_angular_velocity = 0.0; // 可选：初始角速度（rad/s）（如果需要旋转入水，可设置）
Vecd InitialVelocity(initial_speed *cos(initial_angle), initial_speed *sin(initial_angle));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;     /**< Fluid density. */
Real rho0_s = 2.063;   /**< Cylinder density. */
Real gravity_g = 9.81; /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
Real mu_f = 8.9e-7;      /**< Water dynamics viscosity. */
//----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
//std::string diffusion_species_name = "Phi";
//Real diffusion_coeff = 780 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
//Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
//Real cylinder_moisture = 0.0;                                /**< cylinder moisture. */
//Real wall_moisture = 1.0;                                    /**< wall moisture. */


//----------------------------------------------------------------------
// 辅助函数：计算2D多边形的质心（面积加权法）
//----------------------------------------------------------------------
Vecd calculatePolygonCentroid(const std::vector<Vecd> &vertices)
{
    Real area = 0.0;
    Vecd centroid(0.0, 0.0);
    size_t n = vertices.size();

    for (size_t i = 0; i < n - 1; ++i)
    {
        Real x_i = vertices[i][0], y_i = vertices[i][1];
        Real x_j = vertices[i + 1][0], y_j = vertices[i + 1][1];
        Real cross = x_i * y_j - x_j * y_i;

        area += cross;
        centroid[0] += (x_i + x_j) * cross;
        centroid[1] += (y_i + y_j) * cross;
    }

    area *= 0.5;
    if (fabs(area) < 1e-9)
        return Vecd(0.0, 0.0); // 避免除零

    centroid /= (6.0 * area);
    return centroid;
}
//----------------------------------------------------------------------
// 辅助函数：获取缩放后的原始多边形点
//----------------------------------------------------------------------
std::vector<Vecd> getScaledRawPoints()
{
    std::vector<Vecd> raw_points_unscaled = {
        Vecd(2.780, 0.1055),
        Vecd(2.780, -0.1055), // 先向下
        Vecd(0.843, -0.2925),
        Vecd(-1.249, -0.2925),
        Vecd(-1.249, -0.180),
        Vecd(-2.016, -0.180),
        Vecd(-2.016, 0.180), // 再向上
        Vecd(-1.249, 0.180),
        Vecd(-1.249, 0.2925),
        Vecd(0.843, 0.2925),
        Vecd(2.780, 0.1055) // 闭合
    };

    // 缩放至实际尺寸
    std::vector<Vecd> scaled_points;
    for (const auto &p : raw_points_unscaled)
    {
        scaled_points.push_back(Vecd(p[0] * scale_factor, p[1] * scale_factor));
    }
    return scaled_points;
}

//----------------------------------------------------------------------
// 计算前段中心（入水点）的初始位置（缩放+旋转+平移）
//----------------------------------------------------------------------
Vecd calculateFrontCenterPosition()
{
    // 原始前段中心坐标
    Real front_center_raw_x = 2.780;
    Real front_center_raw_y = -0.1055;
    // 1. 缩放
    Vecd front_center_scaled(front_center_raw_x * scale_factor, front_center_raw_y * scale_factor);
    // 2. 计算原始多边形缩放后的质心
    std::vector<Vecd> raw_points = getScaledRawPoints();
    Vecd raw_centroid = calculatePolygonCentroid(raw_points);
    // 3. 计算前段中心相对质心的偏移
    Vecd front_center_relative = front_center_scaled - raw_centroid;
    // 4. 应用初始旋转
    Rotation2d rotation(initial_rotation_angle);
    Vecd front_center_rotated = rotation * front_center_relative;
    // 5. 平移到圆柱全局中心
    Vecd front_center_final = front_center_rotated + cylinder_center;
    return front_center_final;
}

//----------------------------------------------------------------------
// 坐标变换函数 - 全局力转换为圆柱自身坐标系力
//----------------------------------------------------------------------

Real transformGlobalForceToLocalX(const Vecd &global_force, Real rotation_angle)
{
    // 2D旋转矩阵的转置（全局→局部坐标系）
    Real cos_theta = std::cos(rotation_angle);
    Real sin_theta = std::sin(rotation_angle);
    return global_force[0] * cos_theta + global_force[1] * sin_theta; // 局部X方向力 = F_global_x * cosθ + F_global_y * sinθ
}
//----------------------------------------------------------------------
// 获取圆柱实时旋转角度（从Simbody State中提取）
//----------------------------------------------------------------------

Real getCylinderRotationAngle(SimTK::MobilizedBody::Planar &tethered_spot, const SimTK::State &state)
{

    SimTK::Rotation rot = tethered_spot.getBodyRotation(state);
    SimTK::Vec3 angles = rot.convertRotationToBodyFixedXYZ();
    Real angle_from_rot = angles[2]; // 绕Z轴的旋转

    return angle_from_rot; // 使用旋转矩阵的角度
}

StdVec<Vecd> observer_location = {cylinder_center}; /**< Displacement observation point. */
StdVec<Vecd> wetting_observer_location =
    {cylinder_center - Vecd(0.0, 0.01)}; /**< wetting observation point. */
// 新增：前段中心观测点初始位置
StdVec<Vecd> front_center_observer_location = {calculateFrontCenterPosition()}; /**< Front center observation point. */

//----------------------------------------------------------------------
// 汇总输出类，将多个观察量合并到一个文件
//----------------------------------------------------------------------
class SummaryOutput
{
  public:
    SummaryOutput(const std::string &filename)
        : output_file_(filename)
    {
        output_file_ << "Time[s]   "
                     << "ViscousForceGlobal_X   ViscousForceGlobal_Y   "
                     << "PressureForceGlobal_X   PressureForceGlobal_Y   "
                     << "ViscousForceLocalX   PressureForceLocalX   TotalForceLocalX   "
                     << "FrontCenterObserver_Position_X   FrontCenterObserver_Position_Y\n";
    }

    ~SummaryOutput()
    {
        output_file_.close();
    }

    // 扩展参数：同时接收全局力和局部力
    void writeData(Real time,
                   const Vecd &total_viscous_force_global,
                   const Vecd &total_pressure_force_global,
                   Real viscous_local_x,
                   Real pressure_local_x,
                   const Vecd &front_center_position)
    {
        Real total_force_local_x = viscous_local_x + pressure_local_x;

        output_file_ << std::scientific << std::setprecision(9)
                     << time << " "
                     << total_viscous_force_global[0] << " " << total_viscous_force_global[1] << " "
                     << total_pressure_force_global[0] << " " << total_pressure_force_global[1] << " "
                     << viscous_local_x << " " << pressure_local_x << " " << total_force_local_x << " "
                     << front_center_position[0] << " " << front_center_position[1] << "\n";
        output_file_.flush();
    }

  private:
    std::ofstream output_file_;
};

//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block;
    water_block.push_back(Vecd(0.0, 0.0));
    water_block.push_back(Vecd(0.0, LH));
    water_block.push_back(Vecd(DL, LH));
    water_block.push_back(Vecd(DL, 0.0));
    water_block.push_back(Vecd(0.0, 0.0));

    return water_block;
}
class WettingFluidBody : public MultiPolygonShape
{
  public:
    explicit WettingFluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

//class WettingFluidBodyInitialCondition : public LocalDynamics
//{
//  public:
//    explicit WettingFluidBodyInitialCondition(SPHBody &sph_body)
//        : LocalDynamics(sph_body),
//          pos_(particles_->getVariableDataByName<Vecd>("Position")),
//          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};
//
//    void update(size_t index_i, Real dt)
//    {
//        phi_[index_i] = fluid_moisture;
//    };
//
//  protected:
//    Vecd *pos_;
//    Real *phi_;
//};
//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall;
    outer_wall.push_back(Vecd(-BW, -BW));
    outer_wall.push_back(Vecd(-BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, -BW));
    outer_wall.push_back(Vecd(-BW, -BW));

    return outer_wall;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall;
    inner_wall.push_back(Vecd(0.0, 0.0));
    inner_wall.push_back(Vecd(0.0, DH));
    inner_wall.push_back(Vecd(DL, DH));
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));

    return inner_wall;
}
class WettingWallBody : public MultiPolygonShape
{
  public:
    explicit WettingWallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
//class WettingWallBodyInitialCondition : public LocalDynamics
//{
//  public:
//    explicit WettingWallBodyInitialCondition(SPHBody &sph_body)
//        : LocalDynamics(sph_body),
//          pos_(particles_->getVariableDataByName<Vecd>("Position")),
//          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};
//
//    void update(size_t index_i, Real dt)
//    {
//        phi_[index_i] = wall_moisture;
//    };
//
//  protected:
//    Vecd *pos_;
//    Real *phi_;
//};

//----------------------------------------------------------------------
//	Definition for cylinder body
//----------------------------------------------------------------------
class WettingCylinderBody : public MultiPolygonShape
{
  public:
    explicit WettingCylinderBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Rotation2d rotation(initial_rotation_angle);
        std::vector<Vecd> raw_points = getScaledRawPoints();
        Vecd raw_centroid = calculatePolygonCentroid(raw_points);

        // 旋转+平移逻辑
        std::vector<Vecd> polygon_points;
        for (const auto &raw_point : raw_points)
        {
            Vecd point_relative = raw_point - raw_centroid; // 原始点相对自身质心
            Vecd rotated_point = rotation * point_relative; // 绕自身质心旋转
            // 旋转后的点 直接平移到 全局目标中心（cylinder_center）
            Vecd final_point = rotated_point + cylinder_center;
            polygon_points.push_back(Vecd(final_point));
        }

        multi_polygon_.addAPolygon(polygon_points, ShapeBooleanOps::add);
    }
};

//class WettingCylinderBodyInitialCondition : public LocalDynamics
//{
//  public:
//    explicit WettingCylinderBodyInitialCondition(SPHBody &sph_body)
//        : LocalDynamics(sph_body),
//          pos_(particles_->getVariableDataByName<Vecd>("Position")),
//          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};
//
//    void update(size_t index_i, Real dt)
//    {
//        phi_[index_i] = cylinder_moisture;
//    };
//
//  protected:
//    Vecd *pos_;
//    Real *phi_;
//};

//----------------------------------------------------------------------
//	The diffusion model of wetting
//----------------------------------------------------------------------
using CylinderFluidDiffusionDirichlet =
    DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;
//------------------------------------------------------------------------------
// Constrained part for Simbody
//------------------------------------------------------------------------------
MultiPolygon createSimbodyConstrainShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    Rotation2d rotation(initial_rotation_angle);
    std::vector<Vecd> raw_points = getScaledRawPoints();
    Vecd raw_centroid = calculatePolygonCentroid(raw_points);

    std::vector<Vecd> polygon_points;
    for (const auto &raw_point : raw_points)
    {
        Vecd point_relative = raw_point - raw_centroid;
        Vecd rotated_point = rotation * point_relative;
        //直接平移到cylinder_center
        Vecd final_point = rotated_point + cylinder_center;
        polygon_points.push_back(Vecd(final_point));
    }

    multi_polygon.addAPolygon(polygon_points, ShapeBooleanOps::add);
    return multi_polygon;
};

//----------------------------------------------------------------------
//	Cylinder 初始速度设置类（基于质心的速度，消除力矩）
//----------------------------------------------------------------------
class CylinderInitialVelocity : public LocalDynamics
{
  public:
    explicit CylinderInitialVelocity(SPHBody &sph_body, const Vecd &initial_vel, Real initial_angular_vel = 0.0)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          angular_vel_(particles_->getVariableDataByName<Real>("AngularVelocity")),
          initial_vel_(initial_vel),
          initial_angular_vel_(initial_angular_vel),
          mass_(particles_->getVariableDataByName<Real>("Mass"))
    {
        calculateMassCentroid();// 计算质心
    }

    void calculateMassCentroid()
    {
        centroid_ = Vecd(0.0, 0.0);
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
        Vecd r = pos_[index_i] - centroid_; // 速度 = 质心速度 + 绕质心的旋转速度

        vel_[index_i] = initial_vel_; // 设置平动速度（所有粒子相同，保证质心速度）

        if (std::abs(initial_angular_vel_) > 1e-9) // 添加绕质心的旋转（如果角速度不为0）
        {
            Vecd tangential_vel = initial_angular_vel_ * Vecd(-r[1], r[0]);
            vel_[index_i] += tangential_vel;
        }

        angular_vel_[index_i] = initial_angular_vel_; // 设置角速度
    }

  protected:
    Vecd *pos_;
    Vecd *vel_;
    Real *angular_vel_;
    Real *mass_;
    Vecd initial_vel_;
    Real initial_angular_vel_;
    Vecd centroid_;
    Real total_mass_;
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.1. 构建SPH系统
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    //sph_system.setRunParticleRelaxation(true);
    //sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.2. 创建体（流体、壁面、圆柱）
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WettingFluidBody>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WettingWallBody>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder(sph_system, makeShared<WettingCylinderBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 1.0);
    cylinder.defineBodyLevelSetShape();

    //cylinder.defineClosure<Solid, IsotropicDiffusion>(
    //    rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff)); 
    cylinder.defineClosure<Solid>(rho0_s);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();

    BaseParticles &cylinder_particles = cylinder.getBaseParticles();
    cylinder_particles.registerStateVariableData<Vecd>("Velocity");
    cylinder_particles.registerStateVariableData<Real>("AngularVelocity");

    ObserverBody cylinder_observer(sph_system, "CylinderObserver");
    cylinder_observer.generateParticles<ObserverParticles>(observer_location);

    ObserverBody wetting_observer(sph_system, "WettingObserver");
    wetting_observer.generateParticles<ObserverParticles>(wetting_observer_location);

    // 前段中心观测体（仿照OWSC的ObserverBody）
    ObserverBody front_center_observer(sph_system, "FrontCenterObserver");
    front_center_observer.generateParticles<ObserverParticles>(front_center_observer_location);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder});
    ContactRelation wetting_observer_contact(wetting_observer, {&cylinder});
    ContactRelation front_center_observer_contact(front_center_observer, {&cylinder});//前段中心观测体与圆柱的接触关系
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(cylinder);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(cylinder);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    //GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    //CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    //SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    //SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
    //SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

 // 修正：创建圆柱初始速度设置的动力学对象
    SimpleDynamics<CylinderInitialVelocity> cylinder_set_initial_velocity(cylinder, InitialVelocity, initial_angular_velocity);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    /*InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_block_contact);*/
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_block_contact); // 移除湿润相关的复杂自由表面指示器，使用原始版本
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    /** Kernel correction matrix and transport velocity formulation. */
    //InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(DynamicsArgs(water_block_inner, 0.9), water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);//修正
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);//修正

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>> pressure_force_from_fluid(cylinder_contact);

     //定义全局坐标系总粘性力/压力力的统计
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force_global(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_global(cylinder, "PressureForceFromFluid");

    // 插值更新前段中心观测点位置
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
    interpolation_front_center_position(front_center_observer_contact, "Position", "Position");
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the fixed spot. */
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    SolidBodyPartForSimbody cylinder_constraint_area(cylinder, makeShared<MultiPolygonShape>(createSimbodyConstrainShape(cylinder), "cylinder"));
    /** Mass properties of the constrained spot. */
    Vecd actual_centroid = cylinder_constraint_area.initial_mass_center_;
    Vecd tethering_point = actual_centroid; // 使用实际计算的质心
    SimTK::MassProperties cylinder_mass_props(
        cylinder_constraint_area.body_part_mass_properties_->getMass(),       // 保留原质量
        SimTKVec3(tethering_point[0], tethering_point[1], 0.0),               // 强制质心为tethering_point
        cylinder_constraint_area.body_part_mass_properties_->getUnitInertia()); // 保留原转动惯量 

    SimTK::Body::Rigid tethered_spot_info(cylinder_mass_props);
    /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], 0.0)),
                                          fixed_spot_info, SimTK::Transform(SimTKVec3(0)));
    /** Mobility of the tethered spot. */
    Vecd displacement0 = cylinder_constraint_area.initial_mass_center_ - tethering_point;


    SimTK::MobilizedBody::Planar tethered_spot(fixed_spot,
                                               SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], 0.0)),
                                               tethered_spot_info, SimTK::Transform(SimTKVec3(0)));

    // discrete forces acting on the bodies.
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-9.81), 0.0), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));
    SimTK::State state = MBsystem.realizeTopology();

    state.updQ()[0] = displacement0[0];       // x位移（从系留点到质心）
    state.updQ()[1] = displacement0[1];       // y位移（从系留点到质心）
    state.updQ()[2] = initial_rotation_angle; // 相对于父体的旋转

    SimTK::Vec3 mobilizer_vel(
        0.0,                     // 绕Z轴角速度设为0
        InitialVelocity[0], // X方向平移速度（全局坐标系）
        InitialVelocity[1]); // Y方向平移速度（全局坐标系）             
    tethered_spot.setU(state, mobilizer_vel);// 设置初始速度（U）：通过Mobilizer的setU方法

    // 设置完Q/U后，需要让Simbody重新感知状态
    MBsystem.realize(state, SimTK::Stage::Velocity);
    MBsystem.realize(state, SimTK::Stage::Acceleration);
    MBsystem.realize(state, SimTK::Stage::Dynamics);

    SimTK::Vec3 com_vel = tethered_spot.getBodyVelocity(state)[1];

    /** Time stepping method for multibody system.*/
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH..
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    //----------------------------------------------------------------------
    //定义自身坐标系X方向受力的输出文件（仿照2d_flow_around_cylinder的RegressionTestDynamicTimeWarping）
    //----------------------------------------------------------------------
    std::ofstream local_force_file("./output/CylinderLocalXForce.dat");
    local_force_file << "# Time \t ViscousForceLocalX \t PressureForceLocalX \t TotalForceLocalX\n";

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");          // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection"); // output for debug
    RestartIO restart_io(sph_system);

    ObservedQuantityRecording<Vecd> write_cylinder_displacement("Position", cylinder_observer_contact);
    //ObservedQuantityRecording<Real> write_cylinder_wetting("Phi", wetting_observer_contact);
    ObservedQuantityRecording<Vecd> write_front_center_position("Position", front_center_observer_contact);// 记录前段中心观测点的位置

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    cylinder_normal_direction.exec();
    //wetting_water_initial_condition.exec();
    //wetting_wall_initial_condition.exec();
    //wetting_cylinder_initial_condition.exec();
    //Real dt_thermal = get_thermal_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity.exec();
    cylinder_set_initial_velocity.exec(); // 执行圆柱初始速度设置

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1;
    int observation_sample_interval = screen_output_interval * 1;
    int restart_output_interval = screen_output_interval * 500;
    Real end_time = 0.02;
    Real output_interval = end_time /20.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_cylinder_displacement.writeToFile(number_of_iterations);
    //write_cylinder_wetting.writeToFile(number_of_iterations);
    write_front_center_position.writeToFile(number_of_iterations); // 初始时刻输出前段中心位置
    SummaryOutput summary_output("./output/SummaryOutput.dat"); // 创建自定义的汇总输出对象
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;

        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();

            fluid_density_by_summation.exec();
            viscous_force.exec();
            //kernel_correction_complex.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            viscous_force_from_fluid.exec();// 计算流体对圆柱的粘性力（全局）
            pressure_force_from_fluid.exec(); // 计算流体对圆柱的压力力（全局）

            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                //dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                dt = SMIN( fluid_acoustic_time_step.exec(), Dt);//移除wetting
                fluid_pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();// 计算压力力
                viscous_force_from_fluid.exec(); // 计算粘性力
                fluid_density_relaxation.exec(dt);
                //cylinder_wetting.exec(dt);

                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();

                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                constraint_tethered_spot.exec();

                interpolation_front_center_position.exec();//更新前段中心观测点位置

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                SimTK::State &output_state = integ.updAdvancedState();

                write_total_viscous_force_global.writeToFile(number_of_iterations); // 记录结果到文件
                write_total_pressure_force_global.writeToFile(number_of_iterations);

                Real cylinder_rot_angle = getCylinderRotationAngle(tethered_spot, integ.getAdvancedState()); // 获取圆柱实时旋转角度
                //获取粒子数据
                Vec2d total_viscous_force_g(0.0, 0.0);
                Vec2d total_pressure_force_g(0.0, 0.0);
             
                BaseParticles &cylinder_particles = cylinder.getBaseParticles();
                size_t total_real_particles = cylinder_particles.TotalRealParticles();

                auto viscous_force_data = cylinder_particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
                auto pressure_force_data = cylinder_particles.getVariableDataByName<Vecd>("PressureForceFromFluid");
                for (size_t i = 0; i < total_real_particles; ++i)
                {
                    total_viscous_force_g += viscous_force_data[i];
                    total_pressure_force_g += pressure_force_data[i];
                }

                Real current_rotation_from_matrix = getCylinderRotationAngle(tethered_spot, integ.getAdvancedState());
                Real total_rotation_angle = initial_rotation_angle + current_rotation_from_matrix;
                Real angular_velocity = tethered_spot.getU(integ.getAdvancedState())[0];
                Real viscous_local_x = transformGlobalForceToLocalX(total_viscous_force_g, total_rotation_angle);
                Real pressure_local_x = transformGlobalForceToLocalX(total_pressure_force_g, total_rotation_angle);

                // 输出到文件
                local_force_file << std::fixed << std::setprecision(9)
                                 << physical_time << "\t"
                                 << viscous_local_x << "\t"
                                 << pressure_local_x << "\t"
                                 << (viscous_local_x + pressure_local_x) << "\n";
                local_force_file.flush();


                //  获取前段中心观测点位置
                auto &front_center_particles = front_center_observer.getBaseParticles();
                Vecd front_center_pos = front_center_particles.getVariableDataByName<Vecd>("Position")[0];

                // 输出到汇总文件
                summary_output.writeData(physical_time,
                                         total_viscous_force_g,  // 全局粘性力
                                         total_pressure_force_g, // 全局压力力
                                         viscous_local_x,        // 局部粘性力X
                                         pressure_local_x,       // 局部压力力X
                                         front_center_pos);      // 前段中心位置

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_cylinder_displacement.writeToFile(number_of_iterations);
                    //write_cylinder_wetting.writeToFile(number_of_iterations);
                    write_front_center_position.writeToFile(number_of_iterations);// 输出前段中心位置
                    write_total_viscous_force_global.writeToFile(number_of_iterations);// 输出全局坐标系力
                    write_total_pressure_force_global.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            cylinder.updateCellLinkedList();
            water_block_inner.updateConfiguration();// 更新观测体的网格
            cylinder_inner.updateConfiguration();
            cylinder_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            
            front_center_observer_contact.updateConfiguration();//更新前段中心观测点的配置

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
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
};