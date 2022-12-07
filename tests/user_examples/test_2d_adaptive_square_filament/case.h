/**
 * @file 	case.h
 * @brief 	This is the case study the flow induced vibration of a flexible thin structure
 *			attached to a square. 
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

Real DL = 19.5;								/**< Channel length. */
Real DH = 8;								/**< Channel height. */
Vec2d insert_circle_center(5.0, 0.5 * DH);	/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;			/**< Radius of the cylinder. */
Real bh = 0.1 * insert_circle_radius;		/**< Height of the beam. */
Real bl = 8.0 * insert_circle_radius;		/**< Length of the beam. */
Real particle_spacing_ref = 2.0 * bh;				/**< Global reference resolution. */
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
std::vector<Vecd> initial_refinement_region{ Vecd(-DL_sponge - BW, 2.0), Vecd(-DL_sponge - BW, 6.0), Vecd(DL  + BW, 6.0), Vecd(DL  + BW, 2.0), Vecd(-DL_sponge  - BW, 2.0)};
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
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon outer_boundary(water_block_shape);
		add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
		MultiPolygon square(square_shape);
		subtract<MultiPolygonShape>(square);
		MultiPolygon filament(filament_shape);
		subtract<MultiPolygonShape>(filament);
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
Real spacing_ref = 0.5 * particle_spacing_ref;
int particle_number_mid_surface = int( bl / spacing_ref) + 1;
class ShellBaffleParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit ShellBaffleParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_mid_surface; i++)
		{
			Real x = insert_circle_center[0] + insert_circle_radius + Real(i) * spacing_ref - 0.5 * spacing_ref;
			Real y = insert_circle_center[1] + 0.5 * spacing_ref;
            initializePositionAndVolumetricMeasure(Vecd(x, y), spacing_ref);
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
/** Define emitter buffer inflow boundary condition. */
class EmitterBufferInflowCondition : public fluid_dynamics::InflowVelocityCondition
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowCondition(BodyAlignedBoxByCell &aligned_box_part)
		: InflowVelocityCondition(aligned_box_part),
		u_ave_(0), u_ref_(U_f), t_ref_(2.0) {}

	Vecd getPrescribedVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];

		if (position[0] < 0.0)
		{
			u = u_ave_;
			v = 0.0;
		}
		return Vecd(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
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
