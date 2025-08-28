
#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1200.0;
Real Youngs_modulus1 = 1e5;//橡胶杨氏模量
Real Youngs_modulus2 = 1.0e6;//轮毂杨氏模量
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
         // 1. 获取该粒子的全局位置
        Vec2d pos =pos_[index_i];
        // 2. 轮胎/轮毂的圆心位置，与 Geometry 部分保持一致
        Vec2d center(road_length * 0.5, road_height + R_outer);
        // 3. 计算粒子到圆心的距离
        Real dist = (pos - center).norm();
        // 4. 小于 R_inner 的部分视作轮毂（材料索引 1），
        //    否则视作橡胶（材料索引 0）
        if (dist < R_inner)
            material_id_[index_i] = 1;
        else
            material_id_[index_i] = 0;

    };
};

