/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_FLOW_AROUND_CYLINDER_H
#define FVM_FLOW_AROUND_CYLINDER_H
#include "common_weakly_compressible_FVM_classes.h" // classes for weakly compressible fluid only in FVM.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 50.0;                  /**< Channel length. */
Real DH = 30.0;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
Real DL_sponge = 2.0;            /**< Sponge region to impose inflow condition. */
Real DH_sponge = 2.0;            /**< Sponge region to impose inflow condition. */
Real cylinder_radius = 1.0;      /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string ansys_mesh_file_path = "./input/fluent_0.3.msh";
//
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon water_block(createWaterBlockShape());
        add<MultiPolygonShape>(water_block, "WaterBlock");
    }
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
  public:
    FACBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                              vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
        : BoundaryConditionSetupInFVM(inner_relation, each_boundary_type_with_all_ghosts_index,
                each_boundary_type_with_all_ghosts_eij_, each_boundary_type_contact_real_index),
                fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())){};
    virtual ~FACBoundaryConditionSetup(){};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
    }
    void applyFarFieldBoundary(size_t ghost_index) override
    {
        Vecd far_field_velocity(1.0, 0.0);
        Real far_field_density = 1.0;
        Real far_field_pressure = fluid_.getPressure(far_field_density);

        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = far_field_pressure;
        rho_[ghost_index] = far_field_density;
    }

  protected:
    Fluid &fluid_;
};

#endif // EULERIAN_FLOW_AROUND_CYLINDER_H
