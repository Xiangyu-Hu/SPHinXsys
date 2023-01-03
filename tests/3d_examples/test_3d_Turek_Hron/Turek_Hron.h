/**
 * @file 	Turek_Hron.h
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
Real DL = 11.0;							/**< Channel length. */
Real DH = 4.1;							/**< Channel height. */
int resolution = 10;				
Real resolution_ref = 0.2;				/**< Global reference resolution. */
Real BW = resolution_ref * 4.0;			/**< Boundary width, determined by specific layer of boundary particles. */
Real inflow_length = resolution_ref * 20.0;
Vecd insert_cylinder_center(2.0, 2.0, DH /2);	/**< Location of the cylinder center. */
Real insert_cylinder_radius = 0.5;		/**< Radius of the cylinder. */
/** Beam related parameters. */
Real bh = 0.4 * insert_cylinder_radius; /**< Height of the beam. */
Real bl = 7.0 * insert_cylinder_radius; /**< Length of the beam. */
Vecd translation_fluid(DL / 2, DH / 2, DH / 2);
Vecd translation_plate(insert_cylinder_center[0] + (bl + insert_cylinder_radius) /2, insert_cylinder_center[1], insert_cylinder_center[2]);

//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;											  /**< Density. */
Real U_f = 1.0;												  /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;										  /**< Speed of sound. */
Real Re = 100.0;											  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;	/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;

//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
// FLUID
SharedPtr<ComplexShape> createFluid()
{
	auto fluid_shape = makeShared<ComplexShape>("Fluid");
	fluid_shape->add<TriangleMeshShapeBrick>(Vec3d(DL / 2, DH / 2, DH / 2 ), resolution, translation_fluid);
	fluid_shape->subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0,0,1.0), insert_cylinder_radius, DH / 2, resolution, insert_cylinder_center);
	fluid_shape->subtract<TriangleMeshShapeBrick>(Vec3d(bl/2, bh/2 , DH / 2 ), resolution, translation_plate);
	return fluid_shape;
}
//WALL
SharedPtr<ComplexShape> createWall()
{
	auto wall_shape = makeShared<ComplexShape>("Wall");
	wall_shape->add<TriangleMeshShapeBrick>(Vec3d(DL / 2 + BW , DH / 2 + BW, DH / 2 + BW ), resolution, translation_fluid);
	wall_shape->subtract<TriangleMeshShapeBrick>(Vec3d(DL / 2 + BW , DH / 2, DH / 2), resolution, translation_fluid);
	return wall_shape;
}

//inserted body
SharedPtr<ComplexShape> createInsertBody()
{
	auto insert_body = makeShared<ComplexShape>("Body");
	insert_body->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0,0,1.0), insert_cylinder_radius, DH / 2, resolution, insert_cylinder_center);
	insert_body->add<TriangleMeshShapeBrick>(Vec3d(bl/2, bh/2 , DH / 2 ), resolution, translation_plate);
	return insert_body;
}

/** create the beam base as constrain shape. */
SharedPtr<ComplexShape> createBeamBaseShape()
{
	auto constraint = makeShared<ComplexShape>("Constraint");
	constraint->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0,0,1.0), insert_cylinder_radius, DH / 2, resolution, insert_cylinder_center);
	constraint->subtract<TriangleMeshShapeBrick>(Vec3d(bl/2, bh/2 , DH / 2 ), resolution, translation_plate);
	return constraint;
}

struct FreeStreamVelocity
{
	Real u_ref_, t_ref_;

	template <class BoundaryConditionType>
	FreeStreamVelocity(BoundaryConditionType &boundary_condition)
		: u_ref_(U_f), t_ref_(4.0) {}

	Vecd operator()(Vecd &position, Vecd &velocity)
	{
		Vecd target_velocity = Vecd::Zero();
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		target_velocity[0] = 6.0 * u_ave * (DH/2 - position[1]) * (DH/2 + position[1]) / DH / DH;
		return target_velocity;
	}
};

#endif // TUREK_HRON_CASE_H
