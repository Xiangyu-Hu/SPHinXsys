/**
 * @file 	test_3d_FVM_incompressible_channel_flow.h
 * @brief 	This is a test to show inviscid incompressible channel flow using .msh files from ICEM and FLUENT.
 * @author 	Yash Mandaokar, Zhentong Wang and Xiangyu Hu
 */

#ifndef TEST_3D_FVM_INCOMPRESSIBLE_CHANNEL_FLOW_H
#define TEST_3D_FVM_INCOMPRESSIBLE_CHANNEL_FLOW_H

#include "common_weakly_compressible_FVM_classes.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1;            /**< Computation domain length. */
Real DH = 0.6494805454; /**< Computation domain height. */
Real DW = 0.038968832;  /**< Computation domain width. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-0.3, 0.0, 0.0), Vec3d(0.469846, 0.5, 0.03));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Density. */
Real U_f = 1.0;    /**< freestream velocity. */
Real c_f = 10.0 * U_f;
Real mu_f = 0.0; /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_fullpath = "./input/Channel_ICEM.msh";
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------

class AirBody : public ComplexShape
{
  public:
    explicit AirBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wave(0.5 * DH, 0.5 * DL, 0.5 * DW);
        Transform translation_wave(halfsize_wave);
        add<GeometricShapeBox>(Transform(translation_wave), halfsize_wave);
    }
};
///----------------------------------------------------------------------
//	Initialization
//----------------------------------------------------------------------
class InvCFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit InvCFInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body),
          rho_(particles_->registerStateVariableData<Real>("Density")),
          p_(particles_->registerStateVariableData<Real>("Pressure")),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
          mom_(this->particles_->template getVariableDataByName<Vecd>("Momentum")) {};

  protected:
    Real *rho_, *p_;
    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 50 / 117.6655;
        vel_[index_i][0] = 1.0;
        vel_[index_i][1] = 0.0;
        vel_[index_i][2] = 0.0;
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
    }

  protected:
    Vecd *vel_;
    Real *Vol_, *mass_;
    Vecd *mom_;
};
///----------------------------------------------------------------------
//	InvCFBoundaryConditionSetup
//----------------------------------------------------------------------
class InvCFBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
  public:
    InvCFBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
        : BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
          fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())) {};
    virtual ~InvCFBoundaryConditionSetup() {};

    void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
    {
        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij)) + (-e_ij.dot(vel_[index_i]) * (e_ij));
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
    }
    void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
    {
        Vecd far_field_velocity(1.0, 0.0, 0.0);
        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
    }
    void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = vel_[index_i];
        p_[ghost_index] = 100.0 / 117.6655;
        rho_[ghost_index] = rho_[index_i];
    }
    virtual void applySymmetryBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
    {
        vel_[ghost_index] = (vel_[index_i] - 2 * e_ij.dot(vel_[index_i]) * e_ij);
        rho_[ghost_index] = rho_[index_i];
        p_[ghost_index] = p_[index_i];
    }

  protected:
    Fluid &fluid_;
};
#endif // TEST_3D_INCOMPRESSIBLE_CHANNEL_FLOW_H