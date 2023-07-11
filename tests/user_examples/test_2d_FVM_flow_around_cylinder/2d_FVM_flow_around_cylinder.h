/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_FLOW_AROUND_CYLINDER_H
#define FVM_FLOW_AROUND_CYLINDER_H
#include "common_shared_FVM_classes.h"                     // shared classes for weakly-compressible and compressible fluid in FVM.
#include "common_shared_eulerian_classes.h"                // shared eulerian classes for weakly-compressible and compressible fluid.
#include "common_weakly_compressible_FVM_classes.hpp"      // classes for weakly compressible fluid only in FVM.
#include "common_weakly_compressible_eulerian_classes.hpp" // eulerian classes for weakly compressible fluid only.
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
std::string zero_three_flow_around_cylinder_mesh_file_fullpath = "./input/fluent_0.3.msh";
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
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class WeaklyCompressibleFluidInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
    };
    void update(size_t index_i, Real dt){};

  protected:
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> &rho_, &p_;
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public fluid_dynamics::FluidDataInner
{
  public:
    FACBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                              vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index) : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
                                                                                                                                                            Vol_(particles_->Vol_), vel_(particles_->vel_), mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
                                                                                                                                                            fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), total_ghost_particles_(particles_->total_ghost_particles_),
                                                                                                                                                            real_particles_bound_(particles_->real_particles_bound_), each_boundary_type_with_all_ghosts_index_(each_boundary_type_with_all_ghosts_index),
                                                                                                                                                            each_boundary_type_with_all_ghosts_eij_(each_boundary_type_with_all_ghosts_eij_), each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index){};
    virtual ~FACBoundaryConditionSetup(){};

    void resetBoundaryConditions()
    {
        for (size_t boundary_type = 0; boundary_type < each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
        {
            if (!each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
            {
                for (size_t ghost_number = 0; ghost_number != each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                {
                    size_t ghost_index = each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                    size_t index_i = each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                    if (boundary_type == 3)
                    {
                        // non-slip wall boundary
                        vel_[ghost_index] = -vel_[index_i];
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                    }
                    if (boundary_type == 9)
                    {
                        // Far-field boundary
                        Vecd far_field_velocity(1.0, 0.0);
                        Real far_field_density = 1.0;
                        Real far_field_pressure = fluid_.getPressure(far_field_density);

                        vel_[ghost_index] = far_field_velocity;
                        p_[ghost_index] = far_field_pressure;
                        rho_[ghost_index] = far_field_density;
                    }
                }
            }
        }
    };

  protected:
    StdLargeVec<Real> &rho_, &p_, &Vol_;
    StdLargeVec<Vecd> &vel_, &mom_, &pos_;
    Fluid &fluid_;
    size_t &total_ghost_particles_;
    size_t &real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};

#endif // EULERIAN_FLOW_AROUND_CYLINDER_H
