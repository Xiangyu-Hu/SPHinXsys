#include "sphinxsys.h" // SPHinXsys Library.
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"
#include "k-epsilon_turbulent_model.cpp"
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Global parameters on the turbulent properties
//----------------------------------------------------------------------
Real resolution_ref = 0.1;		       	                  /**< Initial reference particle spacing. */
//Real y_p_theo = 0.05;                                   /**< Turbulent: Theoretical distance from the first particle P to wall  */
Real y_p_theo = 0.0;                                      /**< Turbulent: Theoretical distance from the first particle P to wall  */
//Real offset_dist_ref = y_p_theo - 0.5 * resolution_ref; /**< Turbulent: offset distance for keeping y+ unchanged */
Real offset_dist_ref = 0.0;  
/** If this offset value is zero, that means the distance to dummy interface will be real y_p */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;						    		          /**< Reference length. */
Real DH = 16.0 - 2.0 * offset_dist_ref;                   /**< Reference and the height of main channel. */
Real BW = resolution_ref * 4;		                      /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;                     /**< Reference size of the emitter buffer to impose inflow condition. */
Vec2d insert_circle_center(10.0, 0.5 * DH);   /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;              /**< Radius of the cylinder. */
Real characteristic_length = insert_circle_radius * 2.0; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
StdVec<int> id_exclude;
// Observation locations
Vec2d point_coordinate_1(3.0, 5.0);
Vec2d point_coordinate_2(4.0, 5.0);
Vec2d point_coordinate_3(5.0, 5.0);
StdVec<Vecd> observation_locations = { point_coordinate_1, point_coordinate_2, point_coordinate_3 };
//-------------------------------------------------------
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
/** Observation locations, but for channel flow, the cell-based monitoring approach is used, parameters are defined in corresponding .cpp file*/
Real x_observe = 0.90 * DL;
Real x_observe_start = 0.90 * DL;
Real observe_spacing_x = 0.02 * DL;
int num_observer_points_x = 1;

/**Cell method to monitor Center Line data*/
/**Input mannually*/
Real observe_x_ratio = 1.1;
Real observe_x_spacing = resolution_ref * observe_x_ratio;
Real x_start_f = 0.0;
Real x_start_b = insert_circle_center[0] + insert_circle_radius;
int num_observer_points_f = std::round((insert_circle_center[0]- insert_circle_radius) / observe_x_spacing);//** Build observers in front of the cylinder *
int num_observer_points_b = std::round((DL- insert_circle_center[0]- insert_circle_radius) / observe_x_spacing);//** Build observers behind the cylinder *
int num_observer_points = num_observer_points_f + num_observer_points_b;
StdVec<Real> monitor_bound_y = { 7.85 ,8.15 };
/**Input mannually*/
StdVec<Real> monitor_bound_x_f, monitor_bound_x_b;

Real observe_spacing = DH / num_observer_points;
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
//Real Re = 30000.0;					/**< Reynolds number. */
Real Re = 100.0;
Real mu_f = rho0_f * U_f * (2.0 * (insert_circle_radius + 2.0 * offset_dist_ref)) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the emitter block. */
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(-DL_sponge, 0.0),
	Vecd(-DL_sponge, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(-DL_sponge, 0.0),
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
		multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
	}
};
class Cylinder : public MultiPolygonShape
{
public:
	explicit Cylinder(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
	Real u_ref_, t_ref_;

	template <class BoundaryConditionType>
	FreeStreamVelocity(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(2.0) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = Vecd::Zero();
		Real run_time = GlobalStaticVariables::physical_time_;
		target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		return target_velocity;
	}
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
	Real u_ref_, t_ref_;
	AlignedBoxShape& aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(2.0),
		aligned_box_(boundary_condition.getAlignedBox()),
		halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = velocity;
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		if (aligned_box_.checkInBounds(0, position))
		{
			/* Fully-developed velocity inlet */
			//target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
			
			/* Uniform velocity inlet */
			target_velocity[0] = u_ave;
			
			/* Fix velocity in Y direction */
			//target_velocity[1] = 0.0;
		}
		return target_velocity;
	}
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
	Real du_ave_dt_, u_ref_, t_ref_;

public:
	explicit TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

	virtual Vecd InducedAcceleration(const Vecd& position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
	}
};
//----------------------------------------------------------------------
//	Define turbulent inflow boundary condition
//----------------------------------------------------------------------

