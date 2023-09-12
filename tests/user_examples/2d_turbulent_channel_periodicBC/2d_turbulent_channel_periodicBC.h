/**
 * @file 	two_phase_dambreak_static_confinement.h
 * @brief 	Numerical parameters and body definition for 2D two-phase dambreak flow.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 120;						  /**< Reference length. */
Real DH = 2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.1;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW));
/** Observation locations*/
Real x_observe = 0.90 * DL;
Real x_observe_start = 0.90 * DL;
Real observe_spacing_x = 0.02 * DL;
int num_observer_points_x = 1;
int num_observer_points = 20;
Real observe_spacing = DH / num_observer_points;
StdVec<Vecd> observation_locations;
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */             
Real U_f = 1.0;	   /**< Characteristic velocity. */
//Real Re = 40000.0;					/**< Reynolds number. */
Real Re = 500.0;

Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
Real gravity_g = 3.0 * mu_f * U_f / rho0_f / (DH/2.0)/ (DH/2.0) ; /**< Gravity force of fluid. */
Real c_f = 10.0 * U_f;


//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------

/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(-DL_sponge, 0.0),
	Vecd(-DL_sponge, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(-DL_sponge, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-DL_sponge - 2.0 * BW, -BW), //1
	Vecd(-DL_sponge - 2.0 * BW, DH + BW), //2
	Vecd(DL + 2.0 * BW , DH + BW), //3
	Vecd(DL + 2.0 * BW , -BW), //4
	Vecd(-DL_sponge - 2.0 * BW, -BW), //1
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(-DL_sponge - 3.0 * BW, 0.0), //1
	Vecd(-DL_sponge - 3.0 * BW, DH), //2
	Vecd(DL + 3.0 * BW  , DH), //3
	Vecd(DL + 3.0 * BW , 0.0), //4
	Vecd(-DL_sponge - 3.0 * BW, 0.0), //1
};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};

class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};



