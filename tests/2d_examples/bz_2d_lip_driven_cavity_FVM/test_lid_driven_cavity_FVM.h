/**
 * @file 	shock_tube.h
 * @brief 	This is a test to show the standard Sod shock tube case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef SHOCK_TUBE_H
#define SHOCK_TUBE_H
#include "common_shared_FVM_classes.h"                     // shared classes for weakly-compressible and compressible fluid in FVM.
#include "eulerian_fluid_dynamics.hpp"
#include "common_weakly_compressible_FVM_classes.hpp"      // classes for weakly compressible fluid only in FVM.
using namespace SPH; 
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 129.0; /**< Global reference resolution. */
Real BW = resolution_ref * 4;	 /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 400.0;						/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */

std::string lid_driven_1094_mesh_file_fullpath = "./input/fluent.msh";
std::string lid_driven_4282_mesh_file_fullpath = "./input/lid_driven_cavity_4282.msh";
std::string lid_driven_16432_mesh_file_fullpath = "./input/lid_driven_cavity_16432.msh";
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_body_shape;
		water_body_shape.push_back(Vecd(0.0, 0.0));
		water_body_shape.push_back(Vecd(0.0, DH));
		water_body_shape.push_back(Vecd(DL, DH));
		water_body_shape.push_back(Vecd(DL, 0.0));
		water_body_shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public fluid_dynamics::FluidDataInner
{
  public:
    FACBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                              vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index) 
        : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
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
                    if (boundary_type == 4)
                    {
                        Vecd wall_vel = Vecd(1.0, 0.0);
                        vel_[ghost_index] = 2.0 * wall_vel - vel_[index_i];
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
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

#endif // SHOCK_TUBE_H
