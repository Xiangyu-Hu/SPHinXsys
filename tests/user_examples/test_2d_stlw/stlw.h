/**
 * @file 	 stlw.h
 * @brief 	 This is the case file for 2D still water.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 50.0; 				/**< TOTAL SIMULATION TIME*/
Real DL = 3.0;									/**< Tank length. */
Real DH = 2.5;									/**< Tank height. */
Real WH = 2.0;									/**< Water block height. */
Real particle_spacing_ref = 0.05;	
Real BW = particle_spacing_ref * 4.0;			/**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-DL -BW, -DH -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;					   	/**< Reference density of fluid. */
Real gravity_g = 9.81;				       	/**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); 	/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                     	/**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Water block
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(Vec2d(-DL/2,-DH/2));
		water_block_shape.push_back(Vec2d(-DL/2,0));
		water_block_shape.push_back(Vec2d(DL/2,0));
		water_block_shape.push_back(Vec2d(DL/2,-DH/2));
		water_block_shape.push_back(Vec2d(-DL/2,-DH/2));

		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);

	}
};
//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-DL/2, -DH/2)+Vec2d(-BW,-BW));
		outer_wall_shape.push_back(Vec2d(-DL/2, DH/2)+Vec2d(-BW,0));
		outer_wall_shape.push_back(Vec2d(DL/2, DH/2)+Vec2d(+BW,0));
		outer_wall_shape.push_back(Vec2d(DL/2, -DH/2)+Vec2d(+BW,-BW));
		outer_wall_shape.push_back(Vec2d(-DL/2, -DH/2)+Vec2d(-BW,-BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(-DL/2, -DH/2));
		inner_wall_shape.push_back(Vec2d(-DL/2, DH/2));
		inner_wall_shape.push_back(Vec2d(DL/2, DH/2));
		inner_wall_shape.push_back(Vec2d(DL/2, -DH/2));
		inner_wall_shape.push_back(Vec2d(-DL/2, -DH/2));


		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

	}
};
//----------------------------------------------------------------------
//	create mesuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
MultiPolygon createFreeSurfaceGauge()
{	
		/** Geometry definition. */
		std::vector<Vecd> point;
		point.push_back(Vecd(DL/3-h,0));
		point.push_back(Vecd(DL/3-h,DH));
		point.push_back(Vecd(DL/3+h,DH));
		point.push_back(Vecd(DL/3+h,DH));
		point.push_back(Vecd(DL/3-h,0));

		MultiPolygon multi_polygon_;

		multi_polygon_.addAPolygon(point, ShapeBooleanOps::add);

		return multi_polygon_;
}