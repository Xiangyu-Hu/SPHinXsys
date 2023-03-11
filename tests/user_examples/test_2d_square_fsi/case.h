/**
 * @file 	case.h
 * @brief 	This is the case study the flow induced vibration of a flexible thin structure
 *			attached to a square. 
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

Real DL = 19.5;									/**< Channel length. */
Real DH = 8;									/**< Channel height. */
Vec2d insert_circle_center(5.0, 0.5 * DH);		/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;				/**< Radius of the cylinder. */
Real bh = 0.1 * insert_circle_radius;			/**< Height of the beam. */
Real bl = 8.0 * insert_circle_radius;			/**< Length of the beam. */
Real particle_spacing_ref = bh;					/**< Global reference resolution. */
Real DL_sponge = particle_spacing_ref * 20.0; 	/**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0;			/**< Boundary width, determined by specific layer of boundary particles. */
Real thickness = bh;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
/** Observation locations for flow field. */
StdVec<Vecd> observation_locations = {Vec2d(3.0, 5.0), Vec2d(4.0, 5.0), Vec2d(5.0, 5.0)};
/**
 * @brief Material properties of the water and the flexible strucutre. 
 */
// Real rho0_f = 1.18e-3;
// Real U_f = 51.3;
// Real c_f = 10.0 * U_f;
// Real Re = 332.6;
// Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; 

// Real rho0_s = 0.1; 
// Real poisson = 0.45; 
// Real Youngs_modulus =2.5e6;
// Real physical_viscosity = 0.1 * Youngs_modulus; 
Real rho0_f = 1.0;
Real U_f = 1.0;
Real c_f = 10.0 * U_f;
Real Re = 100.0;
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; 

Real rho0_s = 10.0; 
Real poisson = 0.45; 
Real Ae = 5.0e3;	
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
Real physical_viscosity = 0.1 * Youngs_modulus; 
/** Geometry for bodies. */
std::vector<Vecd> water_block_shape{ Vecd(-DL_sponge, 0.0),Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
/** create a filament shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
/** create a gate shape */
std::vector<Vecd> filament_shape{BLB, BLT, BRT, BRB, BLB};
/** create a square shape */
Vec2d SLB(insert_circle_center[0]- insert_circle_radius, insert_circle_center[1] - insert_circle_radius);
Vec2d SLT(insert_circle_center[0]- insert_circle_radius, insert_circle_center[1] + insert_circle_radius);
Vec2d SRB(insert_circle_center[0] + insert_circle_radius, insert_circle_center[1] - insert_circle_radius);
Vec2d SRT(insert_circle_center[0] + insert_circle_radius, insert_circle_center[1] + insert_circle_radius);
std::vector<Vecd> square_shape{SLB, SLT, SRT, SRB, SLB};

Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
/** Water block shape. */
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(filament_shape, ShapeBooleanOps::sub);
		multi_polygon_.addAPolygon(square_shape, ShapeBooleanOps::sub);
	}
};
/* Case-dependent geometries.*/
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(square_shape, ShapeBooleanOps::add);
	}
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int( bl / particle_spacing_ref) + 1;
class ShellBaffleParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit ShellBaffleParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_mid_surface; i++)
		{
			Real x = insert_circle_center[0] + insert_circle_radius + Real(i) * particle_spacing_ref - 0.5 * particle_spacing_ref;
			Real y = insert_circle_center[1] + 0.5 * particle_spacing_ref;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
            Vec2d normal_direction = Vec2d(0, 1.0);
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
		if (base_particles_.pos_[index_i][0] <= insert_circle_center[0] + insert_circle_radius)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
/**
 * Free-stream velocity
*/
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
		target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		return target_velocity;
	}
};
/** Define time dependent acceleration in x-direction. */
class TimeDependentAcceleration : public Gravity
{
	Real du_ave_dt_, u_ref_, t_ref_;

public:
	explicit TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
	}
};