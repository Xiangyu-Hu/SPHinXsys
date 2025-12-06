/**
 * @file 	test_2d_FVM_turbulent_channelflow.h
 * @brief 	This is a test to show high reynolds number turbulent channel flow in FVM.
 * @author 	Yash Mandaokar, Feng Wang and Xiangyu Hu
 */

#ifndef FVM_TURBULENT_CHANNEL_FLOW_H
#define FVM_TURBULENT_CHANNEL_FLOW_H
#include "common_weakly_compressible_FVM_classes.h"
#include "extended_eulerian_riemann_solver.cpp"
#include "rans_turbulence_dynamics.hpp"
#include "turbulence_model.hpp"

using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 120.0;                 /**< Channel length. */
Real DH = 2.0;                   /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
BoundingBoxd system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real U_f = 1.0;
Real rey_bulk = 20000.0;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk; /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                            /**< Reference sound speed. */
Real C_mu = 0.09;

Real turbulent_intensity = 0.05;
Real length_scale = 0.07 * 2 * DH / pow(C_mu, 0.75);
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_file_path = "./input/Channel_Mesh_short.msh";
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
          K_(this->particles_->template registerStateVariableData<Real>("TKE")),
          Eps_(this->particles_->template registerStateVariableData<Real>("Dissipation")),
          mu_t_(this->particles_->template registerStateVariableData<Real>("TurblunetViscosity")),
          rho_(this->particles_->template getVariableDataByName<Real>("Density")),
          p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
          mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
          mom_(this->particles_->template getVariableDataByName<Vecd>("Momentum")) {};

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
        : BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
          K_(this->particles_->template getVariableDataByName<Real>("TKE")),
          Eps_(this->particles_->template getVariableDataByName<Real>("Dissipation")),
          mu_t_(this->particles_->template getVariableDataByName<Real>("TurblunetViscosity")),
          fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
          C_mu_(0.09) {};
    virtual ~TurbulentChannelFlowBoundaryConditionSetup() {};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = K_[index_i];
        Eps_[ghost_index] = Eps_[index_i];
        mu_t_[ghost_index] = mu_t_[index_i];
    }
    // Fully Developed Turbulence Profiles
    std::vector<Real> y_values = {1.95, 1.85, 1.75, 1.65,
                                  1.55, 1.45, 1.35, 1.25,
                                  1.15, 1.05, 0.95, 0.85,
                                  0.75, 0.65, 0.55, 0.45,
                                  0.35, 0.25, 0.15, 0.05};

    std::vector<Real> velocity_magnitude_profile = {0.6799317, 0.86330092, 0.9300949, 0.98939854,
                                                    1.0258716, 1.062102, 1.0885813, 1.1088117,
                                                    1.1224763, 1.1285582, 1.1284137, 1.1219281,
                                                    1.1083114, 1.0883181, 1.0620406, 1.029438,
                                                    0.98861861, 0.9312641, 0.86510468, 0.68170708};

    std::vector<Real> turbulent_kinetic_energy_profile = {0.0061183018, 0.0070746327, 0.0054861265, 0.0046617105,
                                                          0.0040161642, 0.0034167622, 0.0028720065, 0.0024708209,
                                                          0.0021943266, 0.0020446058, 0.0020473369, 0.00221448,
                                                          0.0025275873, 0.0029412457, 0.0033892714, 0.0038914206,
                                                          0.0045933602, 0.0054699243, 0.007050876, 0.0061376584};

    std::vector<Real> turbulent_dissipation_rate_profile = {0.0056344033, 0.0023513641, 0.0011376451, 0.00069459993,
                                                            0.00046991525, 0.00032542294, 0.00022683406, 0.0001672248,
                                                            0.00013027232, 0.00011116108, 0.00011085657, 0.00013127667,
                                                            0.00017384913, 0.00023951771, 0.00032356411, 0.00043955375,
                                                            0.00067517039, 0.0011331175, 0.0023402071, 0.0056611444};
    // Function to find the nearest value based on the y-coordinate
    Real nearestValue(Real y, const std::vector<Real> &y_values, const std::vector<Real> &data_values)
    {
        size_t closest_index = 0;
        Real min_distance = std::abs(y - y_values[0]);

        for (size_t i = 1; i < y_values.size(); ++i)
        {
            Real distance = std::abs(y - y_values[i]);
            if (distance < min_distance)
            {
                closest_index = i;
                min_distance = distance;
            }
        }
        return data_values[closest_index];
    }

    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {
        Real y = pos_[ghost_index][1]; // pos_[ghost_index][1] contains the y-coordinate

        // Find the nearest velocity magnitude and turbulence quantities
        Real inlet_velocity_magnitude = nearestValue(y, y_values, velocity_magnitude_profile);
        Real inlet_tke = nearestValue(y, y_values, turbulent_kinetic_energy_profile);
        Real inlet_eps = nearestValue(y, y_values, turbulent_dissipation_rate_profile);

        // Imposing the inlet profiles
        Vecd inlet_velocity(inlet_velocity_magnitude, 0.0);
        vel_[ghost_index] = inlet_velocity;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        K_[ghost_index] = inlet_tke;
        Eps_[ghost_index] = inlet_eps;
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
        else // Incase of reverse flow
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
    Fluid &fluid_;
    Real C_mu_;
};
#endif // FVM_TURBULENT_CHANNEL_FLOW_H