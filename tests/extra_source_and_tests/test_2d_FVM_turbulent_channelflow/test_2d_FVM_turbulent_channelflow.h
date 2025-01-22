/**
 * @file 	test_2d_FVM_turbulent_channelflow.h
 * @brief 	This is a test to show high reynolds number turbulent channel flow in FVM.
 * @author 	Yash Mandaokar, Feng Wang and Xiangyu Hu
 */

#ifndef FVM_TURBULENT_CHANNEL_FLOW_H
#define FVM_TURBULENT_CHANNEL_FLOW_H             
#include "common_weakly_compressible_FVM_classes.h"
#include "turbulence_model.hpp"
#include "rans_turbulence_dynamics.hpp"
#include "extended_eulerian_riemann_solver.cpp"

using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 120.0;                   /**< Channel length. */
Real DH = 2.0;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real U_f = 1.0;
Real rey_bulk = 20000.0;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk;       /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                               /**< Reference sound speed. */
Real C_mu = 0.09;

Real turbulent_intensity = 0.05;
Real length_scale = 0.07 * 2 * DH / pow(C_mu, 0.75);
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_file_path = "./input/Channel_mesh.msh";
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));

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
//	Initialization
//----------------------------------------------------------------------

class TurbulentChannelFlowInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
    public:
    explicit TurbulentChannelFlowInitialCondition(SPHBody &sph_body)
          : FluidInitialCondition(sph_body), C_mu_(0.09), Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
            K_(this->particles_->template registerStateVariable<Real>("TKE")),
            Eps_(this->particles_->template registerStateVariable<Real>("Dissipation")),
            mu_t_(this->particles_->template registerStateVariable<Real>("TurblunetViscosity")),
            rho_(this->particles_->template getVariableDataByName<Real>("Density")), 
            p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
            mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
            mom_(this->particles_->template getVariableDataByName<Vecd>("Momentum"))
        {};

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 0.2;
        Vecd initial_velocity(1.0, 0.0);
        vel_[index_i] = initial_velocity;
        K_[index_i] = (3.0 / 2.0) * (initial_velocity.squaredNorm()) * pow(turbulent_intensity, 2.0);
        Eps_[index_i] = pow(K_[index_i], 1.5) / length_scale;
        mu_t_[index_i] = C_mu_ * rho_[index_i] * pow(K_[index_i], 2.0) / Eps_[index_i];
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
    }
protected:
  Real C_mu_;
  Real *Vol_, *K_, *Eps_, *mu_t_, *rho_, *p_, *mass_;
  Vecd *mom_;
};
//----------------------------------------------------------------------
//	Case dependent boundary condition
//----------------------------------------------------------------------
class TurbulentChannelFlowBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
public:

    TurbulentChannelFlowBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
        :BoundaryConditionSetupInFVM(inner_relation, ghost_creation), 
        K_(this->particles_->template getVariableDataByName<Real>("TKE")),
        Eps_(this->particles_->template getVariableDataByName<Real>("Dissipation")),
        mu_t_(this->particles_->template getVariableDataByName<Real>("TurblunetViscosity")),
        fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
        C_mu_(0.09){};
    virtual ~TurbulentChannelFlowBoundaryConditionSetup(){};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = K_[index_i];
        Eps_[ghost_index] = Eps_[index_i];
        mu_t_[ghost_index] = mu_t_[index_i];
    }
    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {

        Vecd inlet_velocity(1.0, 0.0);
        vel_[ghost_index] = inlet_velocity;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i]; 
        K_[ghost_index] = (3.0 / 2.0) * (vel_[ghost_index].squaredNorm()) * pow(turbulent_intensity, 2.0);
        Eps_[ghost_index] = pow(K_[ghost_index], 1.5) / length_scale;
        mu_t_[ghost_index] = C_mu_ * rho_[ghost_index] * pow(K_[ghost_index], 2.0) / Eps_[ghost_index];
    }
    void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
    {
        if (vel_[index_i][0] >= 0.0)
        {
            vel_[ghost_index] = vel_[index_i];
            p_[ghost_index] = 0.0;
            rho_[ghost_index] = rho_[index_i];
            K_[ghost_index] = K_[index_i];
            Eps_[ghost_index] = Eps_[index_i];
            mu_t_[ghost_index] = mu_t_[index_i];
        }
        else //Incase of reverse flow
        {
            vel_[ghost_index] = vel_[index_i];
            p_[ghost_index] = 0.0;
            rho_[ghost_index] = rho_[index_i];
            K_[ghost_index] = (3.0 / 2.0) * (vel_[ghost_index].squaredNorm()) * pow(turbulent_intensity, 2.0);
            Eps_[ghost_index] = pow(K_[ghost_index], 1.5) / length_scale;
            mu_t_[ghost_index] = C_mu_ * rho_[ghost_index] * pow(K_[ghost_index], 2.0) / Eps_[ghost_index];
        }
    }
protected:
    Real *K_, *Eps_, *mu_t_;
    Fluid& fluid_;
    Real C_mu_;
};
#endif // FVM_TURBULENT_CHANNEL_FLOW_H