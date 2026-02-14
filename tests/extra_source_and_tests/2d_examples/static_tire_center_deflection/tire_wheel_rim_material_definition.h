
#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
// Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1200.0;
Real Youngs_modulus1 = 3e6; // rubber Young's modulus
Real Youngs_modulus2 = 1e7; // rim (hub) Young's modulus
Real poisson = 0.49;
Real R_outer = 0.300;
Real R_inner = 0.200;
Real road_length = 0.40;
Real road_height = 0.05;

//----------------------------------------------------------------------
// Case-dependent composite material
//----------------------------------------------------------------------
class TireBodyComposite : public CompositeSolid
{
  public:
    TireBodyComposite() : CompositeSolid(rho0_s)
    {
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
    };
};

//----------------------------------------------------------------------
// Case-dependent initialization of material IDs
//----------------------------------------------------------------------
class TireMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit TireMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body){};

    void update(size_t index_i, Real dt = 0.0)
    {
        // 1) Get particle global position
        Vec2d pos = pos_[index_i];
        // 2) Tire/rim center position (consistent with geometry definition)
        Vec2d center(road_length * 0.5, road_height + R_outer);
        // 3) Distance from particle to center
        Real dist = (pos - center).norm();
        // 4) dist < R_inner -> rim (material ID 1); otherwise rubber (material ID 0)
        if (dist < R_inner)
            material_id_[index_i] = 1;
        else
            material_id_[index_i] = 0;
    };
};

