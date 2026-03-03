/**
 * @file 	2d_FVM_double_mach_reflection.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_DOUBLE_MACH_REFLECTION_H
#define FVM_DOUBLE_MACH_REFLECTION_H
#include "common_compressible_FVM_classes.h" // classes for compressible fluid only in FVM.
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0; /**< Computation domain length. */
Real DH = 1.0; /**< Computation domain height. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
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
class DMFInitialCondition : public fluid_dynamics::CompressibleFluidInitialCondition
{
  public:
    explicit DMFInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::CompressibleFluidInitialCondition(sph_body){};
    virtual ~DMFInitialCondition(){};

    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - 1.0 / 6.0))
        {
            /** initial left wave pressure,momentum and energy profile */
            rho_[index_i] = rho0_another;
            mass_[index_i] = rho_[index_i] * Vol_[index_i];
            p_[index_i] = p_another;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_another;
            vel_[index_i][1] = v_another;
            mom_[index_i] = mass_[index_i] * vel_[index_i];
            E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
        }
        else
        {
            rho_[index_i] = rho0_one;
            mass_[index_i] = rho_[index_i] * Vol_[index_i];
            p_[index_i] = p_one;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_one;
            vel_[index_i][1] = v_one;
            mom_[index_i] = mass_[index_i] * vel_[index_i];
            E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
        }
    }

  protected:
    Real gamma_ = heat_capacity_ratio;
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class DMFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
  public:
    DMFBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
        : BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
          E_(particles_->getVariableDataByName<Real>("TotalEnergy")),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")){};
    virtual ~DMFBoundaryConditionSetup(){};

    // Override these methods to define the specific boundary conditions
    void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
    {
        rho_[ghost_index] = rho_[index_i];
        mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * e_ij) - e_ij.dot(vel_[index_i]) * e_ij;
        mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
        p_[ghost_index] = p_[index_i];
        E_[ghost_index] = E_[index_i];
    }

    void applyGivenValueInletFlow(size_t ghost_index) override
    {
        Vecd vel_another(0.0, 0.0);
        vel_another[0] = u_another;
        vel_another[1] = v_another;
        Real p_another = 140.2 / 1.2;
        Real rho_e_another = p_another / (1.4 - 1.0);
        Real mass_another = rho_e_another * Vol_[ghost_index];
        Vecd momentum_another = mass_another * vel_another;
        Real E_inlet_another = (rho_e_another + 0.5 * rho0_another * vel_another.squaredNorm()) * Vol_[ghost_index];

        rho_[ghost_index] = rho0_another;
        mass_[ghost_index] = mass_another;
        p_[ghost_index] = p_another;
        vel_[ghost_index] = vel_another;
        mom_[ghost_index] = momentum_another;
        E_[ghost_index] = E_inlet_another;
    }

    void applyOutletBoundary(size_t ghost_index, size_t index_i) override
    {
        rho_[ghost_index] = rho_[index_i];
        mass_[ghost_index] = mass_[index_i];
        vel_[ghost_index] = vel_[index_i];
        mom_[ghost_index] = mom_[index_i];
        p_[ghost_index] = p_[index_i];
        E_[ghost_index] = E_[index_i];
    }

    void applyTopBoundary(size_t ghost_index, size_t index_i) override
    {
        Real x_1 = 1.0 / 6.0 + *physical_time_ * 10.0 / sin(3.14159 / 3.0);
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - x_1))
        {
            rho_[ghost_index] = rho0_another;
            mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
            p_[ghost_index] = p_another;
            Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
            Vecd vel_another(0.0, 0.0);
            vel_another[0] = u_another;
            vel_another[1] = v_another;
            vel_[ghost_index] = vel_another;
            mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
            E_[ghost_index] = (rho_e + 0.5 * rho_[ghost_index] * vel_[ghost_index].squaredNorm()) * Vol_[ghost_index];
        }
        else
        {
            rho_[ghost_index] = rho0_one;
            mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
            p_[ghost_index] = p_one;
            Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
            Vecd vel_one(0.0, 0.0);
            vel_one[0] = u_one;
            vel_one[1] = v_one;
            vel_[ghost_index] = vel_one;
            mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
            E_[ghost_index] = (rho_e + 0.5 * rho_[ghost_index] * vel_[ghost_index].squaredNorm()) * Vol_[ghost_index];
        }
    }

  protected:
    Real *E_;
    Real *physical_time_;
};
#endif // FVM_DOUBLE_MACH_REFLECTION_H
