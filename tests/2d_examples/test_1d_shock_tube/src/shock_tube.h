/**
 * @file 	shock_tube.h
 * @brief 	This is a test to show the standard Sod shock tube case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#ifndef SHOCK_TUBE_H
#define SHOCK_TUBE_H
#include "common_compressible_eulerian_classes.hpp" // eulerian classes for compressible fluid only.
#include "common_shared_eulerian_classes.h"         // shared eulerian classes for weakly-compressible and compressible fluid.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;                           /**< Tube length. */
Real particle_spacing_ref = 1.0 / 200.0; /**< Initial reference particle spacing. */
Real DH = particle_spacing_ref * 4;      /**< Tube height. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-2.0 / 5.0 * DL, 0.0), Vec2d(3.0 / 5.0 * DL, DH));
Real rho0_l = 1.0;              /**< initial density of left state. */
Real rho0_r = 0.125;            /**< initial density of right state. */
Vecd velocity_l = Vecd::Zero(); /**< initial velocity of left state. */
Vecd velocity_r = Vecd::Zero();
;               /**< initial velocity of right state. */
Real p_l = 1.0; /**< initial pressure of left state. */
Real p_r = 0.1; /**< initial pressure of right state. */
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real heat_capacity_ratio = 1.4; /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Cases-dependent geometry
//----------------------------------------------------------------------
class WaveBlock : public MultiPolygonShape
{
  public:
    explicit WaveBlock(const std::string &body_name)
        : MultiPolygonShape(body_name)
    {
        std::vector<Vecd> waves_block_shape{
            Vecd(-2.0 / 5.0 * DL, 0.0), Vecd(-2.0 / 5.0 * DL, DH), Vecd(3.0 / 5.0 * DL, DH),
            Vecd(3.0 / 5.0 * DL, 0.0), Vecd(-2.0 / 5.0 * DL, 0.0)};
        multi_polygon_.addAPolygon(waves_block_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class ShockTubeInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit ShockTubeInitialCondition(SPHBody &sph_body)
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
        if (pos_[index_i][0] < DL / 10.0)
        {
            // initial left state pressure,momentum and energy profile
            rho_[index_i] = rho0_l;
            p_[index_i] = p_l;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i] = velocity_l;
            mom_[index_i] = rho0_l * velocity_l;
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
        if (pos_[index_i][0] > DL / 10.0)
        {
            // initial right state pressure,momentum and energy profile
            rho_[index_i] = rho0_r;
            p_[index_i] = p_r;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i] = velocity_r;
            mom_[index_i] = rho0_r * velocity_r;
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
    }

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real gamma_;
};
#endif // SHOCK_TUBE_H
