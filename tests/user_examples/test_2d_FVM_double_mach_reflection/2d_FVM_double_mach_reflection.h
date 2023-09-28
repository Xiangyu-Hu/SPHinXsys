/**
 * @file 	2d_FVM_double_mach_reflection.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_DOUBLE_MACH_REFLECTION_H
#define FVM_DOUBLE_MACH_REFLECTION_H
#include "common_compressible_FVM_classes.h"        // eulerian classes for compressible fluid in FVM only.
#include "common_compressible_eulerian_classes.hpp" // eulerian classes for weakly compressible fluid only.
#include "common_shared_FVM_classes.h"              // shared eulerian classes for weakly-compressible and compressible fluid in FVM.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                           /**< Computation domain length. */
Real DH = 1.0;                           /**< Computation domain height. */
Real particle_spacing_ref = 1.0 / 240.0; /**< Initial reference particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_one = 1.4;                         /**< initial density of one fluid. */
Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
Real p_one = 1.0;                            /**< initial pressure of one fluid. */
Real rho0_another = 8.0;                     /**< initial density of another. */
Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string double_mach_reflection_mesh1_fullpath = "./input/double_mach_reflection_0.05.msh";
std::string double_mach_reflection_mesh2_fullpath = "./input/double_mach_reflection_0.0125.msh";
std::string double_mach_reflection_mesh3_fullpath = "./input/double_mach_reflection_0.012.msh";
//
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> CreatComputationDomain()
{
    // geometry
    std::vector<Vecd> computation_domain;
    computation_domain.push_back(Vecd(0.0, 0.0));
    computation_domain.push_back(Vecd(0.0, DH));
    computation_domain.push_back(Vecd(DL, DH));
    computation_domain.push_back(Vecd(DL, 0.0));
    computation_domain.push_back(Vecd(0.0, 0.0));
    return computation_domain;
}
class WaveBody : public ComplexShape
{
  public:
    explicit WaveBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon wave_block(CreatComputationDomain());
        add<MultiPolygonShape>(wave_block, "WaveBlock");
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DMFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit DMFInitialCondition(SPHBody &sph_body)
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
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - 1.0 / 6.0))
        {
            /** initial left wave pressure,momentum and energy profile */
            rho_[index_i] = rho0_another;
            p_[index_i] = p_another;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_another;
            vel_[index_i][1] = v_another;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
        else
        {
            rho_[index_i] = rho0_one;
            p_[index_i] = p_one;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_one;
            vel_[index_i][1] = v_one;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
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

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class DMFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
  public:
    DMFBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                              vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
        : BoundaryConditionSetupInFVM(inner_relation, each_boundary_type_with_all_ghosts_index, 
            each_boundary_type_with_all_ghosts_eij_, each_boundary_type_contact_real_index),
            E_(*particles_->getVariableByName<Real>("TotalEnergy")){};
    virtual ~DMFBoundaryConditionSetup(){};

    // Override these methods to define the specific boundary conditions
    void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
    {
        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * e_ij) - e_ij.dot(vel_[index_i]) * e_ij;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        E_[ghost_index] = E_[index_i];
    }

    void applyGivenValueInletFlow(size_t ghost_index) override
    {
        Vecd vel_another(0.0, 0.0);
        vel_another[0] = u_another;
        vel_another[1] = v_another;
        Real p_another = 140.2 / 1.2;
        Real rho_e_another = p_another / (1.4 - 1.0);
        Real E_inlet_another = rho_e_another + 0.5 * rho0_another * vel_another.squaredNorm();

        rho_[ghost_index] = rho0_another;
        p_[ghost_index] = p_another;
        vel_[ghost_index] = vel_another;
        E_[ghost_index] = E_inlet_another;
    }

    void applyOutletBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
        E_[ghost_index] = E_[index_i];
    }

    void applyTopBoundary(size_t ghost_index, size_t index_i) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real x_1 = 1.0 / 6.0 + run_time * 10.0 / sin(3.14159 / 3.0);
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - x_1))
        {
            rho_[ghost_index] = rho0_another;
            p_[ghost_index] = p_another;
            Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
            Vecd vel_another(0.0, 0.0);
            vel_another[0] = u_another;
            vel_another[1] = v_another;
            vel_[ghost_index] = vel_another;
            E_[ghost_index] = rho_e + 0.5 * rho_[ghost_index] * vel_[ghost_index].squaredNorm();
        }
        else
        {
            rho_[ghost_index] = rho0_one;
            p_[ghost_index] = p_one;
            Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
            Vecd vel_one(0.0, 0.0);
            vel_one[0] = u_one;
            vel_one[1] = v_one;
            vel_[ghost_index] = vel_one;
            E_[ghost_index] = rho_e + 0.5 * rho_[ghost_index] * vel_[ghost_index].squaredNorm();
        }
    }

  protected:
    StdLargeVec<Real> &E_;
};

#endif // FVM_DOUBLE_MACH_REFLECTION_H
