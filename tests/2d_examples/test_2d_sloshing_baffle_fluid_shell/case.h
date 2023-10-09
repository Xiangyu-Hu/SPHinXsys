/**
* @file 	case.h
* @brief 	Fluid-shell interaction in sloshing flow. 
* @details 	Here, the first fluid-shell interaction test is presented. 
* @author 	Chi Zhang
*/

#include "sphinxsys.h"
#define PI 3.1415926
 /**
* @brief Namespace cite here.
*/
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 1.0; 								           /**< Tank length. */
Real DH = 0.7; 							              /**< Tank height. */
Real L_W = 1.0;                                      /**< water width. */
Real L_H = 0.15;                                    /**< water depth. */
Real Gate_x = 0.5 * L_W;						   /**< Width of the gate. */
Real Gate_width = 0.006;						  /**< Width of the gate. */
Real Gate_height = 0.174;						 /**< Height of the gate. */
Real particle_spacing_ref = Gate_width; 	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0; 		   /**< Extending width for BCs. */
Real thickness = Gate_width;
int particle_number_mid_surface = (0.026 + Gate_height + BW) / particle_spacing_ref;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// parameters for liquid sloshing Case:Xue Mian
Real a_0 = 0.01; 					      /**< amplitude of the sloshing. */
Real w_0 = 1.2 * 4.142814038; 			 /**< frequency of the sloshing  0.583*8.9556. */
Real k_0 = a_0 * w_0 * w_0;             /**< parameter of sloshing x = k_0 * sin(2*pi*f_0*t). */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;						/**< Reference density of fluid. */
Real gravity_g = 9.81; 					    /**< Value of gravity. */
Real U_max = 2.0*sqrt(gravity_g*L_H);		/**< Characteristic velocity. */
Real c_f = 10.0* U_max;					    /**< Reference sound speed. */
Real mu_f = 1.0e-6;
/**
 * @brief Material properties of the elastic gate.
 */
Real rho0_s = 1250;					/**< Reference density of gate. */
Real Youngs_modulus = 30.0e6;		/**< reference Youngs modulus. */
Real poisson = 0.47; 				/**< Poisson ratio. */
Real physical_viscosity = 1000.0; 	/**< physical damping, here we choose the same value as numerical viscosity. */
/** create a water block shape */
StdVec<Vecd> water_block_shape{
	Vecd(0.0, 0.0), Vecd(0.0, L_H), Vecd(L_W, L_H), Vecd(L_W, 0.0), Vecd(0.0, 0.0)};
/** Baffle base shape */
StdVec<Vecd> baffle_base_shap{
	Vecd(Gate_x - 0.0262, 0.0),
	Vecd(Gate_x - 0.0262, 0.006),
	Vecd(Gate_x - 0.0262 + 0.0075, 0.006),
	Vecd(Gate_x - 0.0262 + 0.0075, 0.0125),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.0125),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.026),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022, 0.026),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022, 0.0125),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008, 0.0125),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008, 0.006),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008 + 0.0075, 0.006),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008 + 0.0075, 0.0),
	Vecd(Gate_x - 0.0262, 0.0)
};
/** Inner wall shape */
StdVec<Vecd> inner_wall_shape{
    Vecd(0.0, 0.0),
	Vecd(0.0, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(0.0, 0.0)
};
/** Outer wall shape */
StdVec<Vecd> outer_wall_shape{
	Vecd(-BW, -BW),
	Vecd(-BW, DH + BW),
	Vecd(DL + BW, DH + BW),
	Vecd(DL + BW, -BW),
	Vecd(-BW, -BW)
};
/** Baffle shape. */
StdVec<Vecd> baffle_shape{
    Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002, -BW),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002, 0.026 + Gate_height),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006, 0.026 + Gate_height),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006, -BW),
	Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002, -BW)
};
/* Case-dependent geometries.*/
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(baffle_shape, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(baffle_base_shap, ShapeBooleanOps::sub);
	}
};
/* Case-dependent geometries.*/
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(baffle_base_shap, ShapeBooleanOps::add);
	}
};
/** Particle generator and constraint boundary for shell baffle. */
class ShellBaffleParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit ShellBaffleParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_mid_surface; i++)
		{
			Real x = Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 - 0.4 * particle_spacing_ref; //0.5, 0.75, 1.5
			Real y = -BW + Real(i) * particle_spacing_ref;
			initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
			Vec2d normal_direction = Vec2d(1.0, 0);
			initializeSurfaceProperties(normal_direction, thickness);
		}
	}
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][1] < 0.026)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
/**
 * @brief 	define external force.
 */
class VariableGravity : public Gravity
{
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g)) {};
	void UpdateAcceleration()
	{

		if (GlobalStaticVariables::physical_time_ <= 1.0)
		{
			global_acceleration_[0] = 0.0;
		}

		if (GlobalStaticVariables::physical_time_ > 1.0)
		{
			global_acceleration_[0]
				= k_0 * cos(w_0 * (GlobalStaticVariables::physical_time_ - 1.0));
		}

	}
};
/**
 * @brief  Setup for gate and fluid boservers. 
 */
Real h = 1.37 * particle_spacing_ref;
Real probe_x1 = DL - 0.032;
MultiPolygon CreateWaveProbeShape1()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x1 - h, 0.0));
	pnts.push_back(Vecd(probe_x1 - h, DH));
	pnts.push_back(Vecd(probe_x1 + h, DH));
	pnts.push_back(Vecd(probe_x1 + h, 0.0));
	pnts.push_back(Vecd(probe_x1 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}

Real probe_x2 = 0.032;
MultiPolygon  CreateWaveProbeShape2()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x2 - h, 0.0));
	pnts.push_back(Vecd(probe_x2 - h, DH));
	pnts.push_back(Vecd(probe_x2 + h, DH));
	pnts.push_back(Vecd(probe_x2 + h, 0.0));
	pnts.push_back(Vecd(probe_x2 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}

Real probe_x3 = DL / 2.0;
MultiPolygon  CreateWaveProbeShape3()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x3 - h, 0.0));
	pnts.push_back(Vecd(probe_x3 - h, DH));
	pnts.push_back(Vecd(probe_x3 + h, DH));
	pnts.push_back(Vecd(probe_x3 + h, 0.0));
	pnts.push_back(Vecd(probe_x3 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}
/**
 * @brief Define the observer body.
 */
StdVec<Vecd> baffle_disp_probe_location{Vecd(Gate_x, 0.195), Vecd(Gate_x, 0.145), Vecd(Gate_x, 0.095), Vecd(Gate_x, 0.05)};
StdVec<Vecd> baffle_pressure_probe_location{Vecd(Gate_x - Gate_width / 2, 0.19), Vecd(Gate_x + Gate_width / 2, 0.19)};    
StdVec<Vecd> fluid_pressure_probe_location{Vecd(DL, 0.045), Vecd(DL, 0.095), Vecd(DL, 0.165), Vecd(DL, 0.22)};