/**
 * @file 	freestream_flow_around_cylinder.h
 * @brief 	This is the case file for the test of free-stream boundary condition.
 * @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 19.5;                                    /**< Channel length. */
Real DH = 8.0;                                     /**< Channel height. */
Real particle_spacing_ref = 0.1;                   /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0;      /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;              /**< Sponge region to impose injection. */
Vec2d insert_square_center(5.0, 0.5 * DH);         /**< Location of the square center. */
Real insert_square_halflength = 0.5;               /**< Halflength of the square. */
Real particle_spacing_beam = particle_spacing_ref; /**< Initial beam particle spacing. */
Real bh = particle_spacing_beam;                   /**< Height of the beam. */
Real bl = 8.0 * insert_square_halflength;          /**< Length of the beam. */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                                /**< Density. */
Real U_f = 1.0;                                                   /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                            /**< Speed of sound. */
Real Re = 100.0;                                                  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_square_halflength) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 10.0;
Real poisson = 0.45;
Real Ae = 5e3;
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
Real physical_viscosity = 0.1 * Youngs_modulus;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
//	water block shape
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
/** create a square shape */
Vec2d SLB(insert_square_center[0] - insert_square_halflength, insert_square_center[1] - insert_square_halflength);
Vec2d SLT(insert_square_center[0] - insert_square_halflength, insert_square_center[1] + insert_square_halflength);
Vec2d SRB(insert_square_center[0] + insert_square_halflength, insert_square_center[1] - insert_square_halflength);
Vec2d SRT(insert_square_center[0] + insert_square_halflength, insert_square_center[1] + insert_square_halflength);
std::vector<Vecd> square_shape{SLB, SLT, SRT, SRB, SLB};
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_square_center[0] + insert_square_halflength, insert_square_center[1] - hbh);
Vec2d BLT(insert_square_center[0] + insert_square_halflength, insert_square_center[1] + hbh);
Vec2d BRB(insert_square_center[0] + insert_square_halflength + bl, insert_square_center[1] - hbh);
Vec2d BRT(insert_square_center[0] + insert_square_halflength + bl, insert_square_center[1] + hbh);
std::vector<Vecd> beam_shape{BLB, BLT, BRT, BRB, BLB};
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(square_shape, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::sub);
    }
};
class Square : public MultiPolygonShape
{
  public:
    explicit Square(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(square_shape, ShapeBooleanOps::add);
    }
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int(bl / particle_spacing_beam);
class BeamParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit BeamParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real x = insert_square_center[0] + insert_square_halflength + (Real(i) + 0.5) * particle_spacing_beam;
            Real y = insert_square_center[1];
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_beam);
            Vec2d normal_direction = Vec2d(0, 1.0);
            initializeSurfaceProperties(normal_direction, bh);
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
        if (base_particles_.pos_[index_i][0] <= insert_square_center[0] + insert_square_halflength + 2 * particle_spacing_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
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
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real t_ref_, u_ref_, du_ave_dt_;

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
