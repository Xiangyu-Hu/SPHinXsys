#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1200.0;
Real Youngs_modulus1 = 3e6;//橡胶杨氏模量
Real Youngs_modulus2 = 1e7;//轮毂杨氏模量
Real poisson = 0.49;
Real R_outer = 0.300;
Real R_inner = 0.200;
Real road_length = 0.400;
Real road_height = 0.05;
//----------------------------------------------------------------------
//	Case dependent composite material
//----------------------------------------------------------------------
class TireBodyComposite : public CompositeSolid
{
  public:
    TireBodyComposite() : CompositeSolid(rho0_s)
    {
        //add<ActiveModelSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
    };
};
//----------------------------------------------------------------------
//	Case dependent initialization material ids
//----------------------------------------------------------------------
class TireMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit TireMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body){};

    void update(size_t index_i, Real dt = 0.0)
{
    // 1. 获取粒子位置
    Vec3d pos = pos_[index_i];

    // 2. 定义圆环中心（与几何一致）
    Vec3d center(road_length * 0.5, road_height + R_outer, 0.0);

    // 3. 计算径向距离（忽略 Z）
    Real radial_dist = (Vec2d(pos[0], pos[1]) - Vec2d(center[0], center[1])).norm();

    // 4. 计算 Z 方向距离（轴向方向）
    Real z_dist = std::abs(pos[2] - center[2]);

    // 5. 轴向范围的一半（由 tire_thickness 定义）
    Real half_thickness = 0.5 * 0.2;

    // 6. 判断是否为轮毂部分（内径 + Z方向限制）
    if (radial_dist < R_inner && z_dist < half_thickness)
        material_id_[index_i] = 1; // 轮毂
    else
        material_id_[index_i] = 0; // 橡胶
}
};

