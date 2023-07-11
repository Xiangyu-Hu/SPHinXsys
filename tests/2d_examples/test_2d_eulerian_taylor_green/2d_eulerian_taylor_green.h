/**
 * @file 	eulerian_taylor_green.h
 * @brief 	This is a test to show the standard eluerian taylor green case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_TAYLOR_GREEN_H
#define EULERIAN_TAYLOR_GREEN_H
#include "common_compressible_eulerian_classes.hpp" // eulerian classes for compressible fluid only.
#include "common_shared_eulerian_classes.h"         // shared eulerian classes for weakly-compressible and compressible fluid.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;                    /**< box length. */
Real DH = 1.0;                    /**< box height. */
Real resolution_ref = 1.0 / 50.0; /**< Global reference resolution. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                  /**< Reference density of fluid. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Reference sound speed. */
Real Re = 100;                      /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
Real heat_capacity_ratio = 1.4;     /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
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
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class TaylorGreenInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit TaylorGreenInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
          rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
        particles_->registerVariable(E_, "TotalEnergy");
        particles_->registerVariable(dE_dt_, "TotalEnergyChangeRate");
        particles_->registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
        gamma_ = heat_capacity_ratio;
    };

    void update(size_t index_i, Real dt)
    {
        /** initial momentum and energy profile */
        rho_[index_i] = rho0_f;
        p_[index_i] = pow(c_f, 2) * rho_[index_i] / gamma_;
        vel_[index_i][0] = -cos(2.0 * Pi * pos_[index_i][0]) *
                           sin(2.0 * Pi * pos_[index_i][1]);
        vel_[index_i][1] = sin(2.0 * Pi * pos_[index_i][0]) *
                           cos(2.0 * Pi * pos_[index_i][1]);
        mom_[index_i] = rho_[index_i] * vel_[index_i];
        Real rho_e = p_[index_i] / (gamma_ - 1.0);
        E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real gamma_;
};

#endif // EULERIAN_TAYLOR_GREEN_H
