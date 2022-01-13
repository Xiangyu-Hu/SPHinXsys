/**
* @file 	2d_cylinder_flow.h
* @brief 	This is the case file for the test of flow passing by a cylinder.
* @details  We consider a flow passing by a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;							/**< Channel length. */
Real DH = 10.0;							/**< Channel height. */
Real resolution_ref = 0.2;				/**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 10.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0;	/**< Sponge region to impose freestream condition. */
Vec2d insert_circle_center(4.0, 5.0);	/**< Location of the cylinder center. */
Real insert_circle_radius = 0.75;		/**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;											  /**< Density. */
Real U_f = 1.0;												  /**< freestream velocity. */
Real c_f = 10.0 * U_f;										  /**< Speed of sound. */
Real Re = 100.0;											  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	SPH bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	return water_block_shape;
}
/** create a water block buffer shape. */
MultiPolygon createBufferShape()
{
	std::vector<Vecd> buffer_shape;
	buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	buffer_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	buffer_shape.push_back(Vecd(DL, DH + DH_sponge));
	buffer_shape.push_back(Vecd(DL, DH));
	buffer_shape.push_back(Vecd(0.0, DH));
	buffer_shape.push_back(Vecd(0.0, 0.0));
	buffer_shape.push_back(Vecd(DL, 0.0));
	buffer_shape.push_back(Vecd(DL, -DH_sponge));
	buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Define case dependent bodies material, constraint and boundary conditions.
//----------------------------------------------------------------------
/** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Define parametrization for this case.
//----------------------------------------------------------------------
class ParameterizedWaterMaterial : public BaseParameterization<WeaklyCompressibleFluid>
{
public:
	ParameterizedWaterMaterial(ParameterizationIO &parameterization_io, Real rho0, Real c0, Real mu)
		: BaseParameterization<WeaklyCompressibleFluid>(parameterization_io, rho0, c0, mu)
	{
		getAParameter("WaterMaterial", "Viscosity", mu_);
	}
};
/** Definition of the cylinder. */
class Cylinder : public SolidBody
{
public:
	Cylinder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
//----------------------------------------------------------------------
//	Case dependent flow boundary condition.
//----------------------------------------------------------------------
class FreeStreamCondition : public fluid_dynamics::FlowRelaxationBuffer
{
	Real u_ave_, u_ref_, t_ref;

public:
	FreeStreamCondition(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: fluid_dynamics::FlowRelaxationBuffer(fluid_body, constrained_region),
		u_ave_(0), u_ref_(U_f), t_ref(2.0) {}
	Vecd getTargetVelocity(Vecd &position, Vecd &velocity)
	{
		return Vecd(u_ave_, 0.0);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
	}
};
//----------------------------------------------------------------------
//	Observer particle generator.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** the measureing particles */
		Vec2d point_coordinate_1(3.0, 5.0);
		Vec2d point_coordinate_2(4.0, 5.0);
		Vec2d point_coordinate_3(5.0, 5.0);

		positions_volumes_.push_back(std::make_pair(point_coordinate_1, 0.0));
		positions_volumes_.push_back(std::make_pair(point_coordinate_2, 0.0));
		positions_volumes_.push_back(std::make_pair(point_coordinate_3, 0.0));
	}
};
