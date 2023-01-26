/**
 * @file 	3d_channel_flow.h
 * @brief 	This is the case file for the test of fluid - structure interaction.
 * @details  We consider a flow - induced vibration of an elastic plate behind a cylinder in 3D.
 * @author 	Anastazja Broniatowska
 */

#ifndef TEST_3D_TUREK_HRON_CASE_H
#define TEST_3D_TUREK_HRON_CASE_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
int resolution = 25;				
Real resolution_ref = 0.2;				/**< Global reference resolution. */

//CHANEL
Real DL = 11.0;							/**< Channel length. */
Real DH = 4.1;
Real DW = 4.1;						/**< Channel width. */
Real BW = resolution_ref * 4.0;		/**< Boundary width, determined by specific layer of boundary particles. */
Real inflow_length = resolution_ref * 20;		
Vecd translation_fluid(DL / 2, DH / 2, DW / 2);

//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1000.0;											  /**< Density. */
Real U_f = 6.0;												  /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;										  /**< Speed of sound. */
Real Re = 100.0;											  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
// FLUID
SharedPtr<ComplexShape> createFluid()
{
	auto fluid_shape = makeShared<ComplexShape>("Fluid");
	fluid_shape->add<TriangleMeshShapeBrick>(Vec3d(DL / 2, DH / 2, DW / 2 ), resolution, translation_fluid);
	return fluid_shape;
}
//WALL
SharedPtr<ComplexShape> createWall()
{
	auto wall_shape = makeShared<ComplexShape>("Wall");
	wall_shape->add<TriangleMeshShapeBrick>(Vec3d(DL / 2 + BW , DH / 2 + BW, DW / 2 + BW ), resolution, translation_fluid);
	wall_shape->subtract<TriangleMeshShapeBrick>(Vec3d(DL / 2 + BW , DH / 2 , DW / 2), resolution, translation_fluid);
	return wall_shape;
}

struct FreeStreamVelocity
{
	Real u_ref_, t_ref_;

	template <class BoundaryConditionType>
	FreeStreamVelocity(BoundaryConditionType &boundary_condition)
		: u_ref_(U_f), t_ref_(2.0) {}

	Vecd operator()(Vecd &position, Vecd &velocity)
	{
		Vecd target_velocity = Vecd::Zero();
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		target_velocity[0] = u_ave * (DH/2 - position[1]) * (DH/2 + position[1]) / DH / DH;
		return target_velocity;
	}
};


#endif // TUREK_HRON_CASE_H
