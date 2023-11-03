/**
 * @file 	case.h
 * @brief 	This is the case study the flow in the channel.
 * @author  Chi Zhang & Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_), inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
        {
            contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
        inner_relation.base_particles_.registerVariable(dW_ijV_je_ij_ttl, "TotalKernelGrad");
        inner_relation.base_particles_.registerVariable(number_of_inner_neighbor, "InnerNeighborNumber");
        inner_relation.base_particles_.registerVariable(number_of_contact_neighbor, "ContactNeighborNumber");
    }

    inline void exec()
    {
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                int N_inner_number = 0;
                int N_contact_number = 0;
                Real W_ijV_j_ttl_i = 0;
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }
};

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real DH = 12.5 * scale;                 /**< Channel height. */
Real DL = DH * 9.0;                     /**< Channel length. */
Real resolution_ref = DH / Real(25);    /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */
Real valve_x = 0.5 * DL;                /**<Position of the valve>*/
Real resolution_valve = resolution_ref;
Real valve_thickness = resolution_valve;
Real valve_length = DH - 0.5 * resolution_valve;
Real base_length = 4 * resolution_valve;
/** Inflow shape*/
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1100.0;               // blood density
const Real mu_f = 3.6e-3;                 // blood viscosity
const Real Re = 2000;                     /**< Reynolds number. */
const Real U_f = Re * mu_f / rho0_f / DH; // average velocity at peak systole
const Real U_max = 1.5 * U_f;
const Real c_f = 10.0 * U_max; /**< Speed of sound. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 1000.0;
const Real Youngs_modulus = 3.6e6; // use softer material for development
const Real poisson = 0.49;
const Real shape_constant = 0.4;
const Real physical_viscosity = 1e3 * shape_constant / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * valve_thickness; /**< physical damping. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.emplace_back(-DL_sponge, 0.0);
    water_block_shape.emplace_back(-DL_sponge, DH);
    water_block_shape.emplace_back(DL, DH);
    water_block_shape.emplace_back(DL, 0.0);
    water_block_shape.emplace_back(-DL_sponge, 0.0);

    return water_block_shape;
}
/** create a valve shape */
std::vector<Vecd> createValveShape()
{
    std::vector<Vecd> valve_shape;
    valve_shape.emplace_back(valve_x - 0.5 * valve_thickness, DH);
    valve_shape.emplace_back(valve_x + 0.5 * valve_thickness, DH);
    valve_shape.emplace_back(valve_x + 0.5 * valve_thickness, DH - valve_length);
    valve_shape.emplace_back(valve_x - 0.5 * valve_thickness, DH - valve_length);
    valve_shape.emplace_back(valve_x - 0.5 * valve_thickness, DH);

    return valve_shape;
}
/** Inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.emplace_back(-DL_sponge - BW, 0.0);
    inner_wall_shape.emplace_back(-DL_sponge - BW, DH);
    inner_wall_shape.emplace_back(DL + BW, DH);
    inner_wall_shape.emplace_back(DL + BW, 0.0);
    inner_wall_shape.emplace_back(-DL_sponge - BW, 0.0);

    return inner_wall_shape;
}
/** Outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.emplace_back(-DL_sponge - BW, -BW);
    outer_wall_shape.emplace_back(-DL_sponge - BW, DH + BW);
    outer_wall_shape.emplace_back(DL + BW, DH + BW);
    outer_wall_shape.emplace_back(DL + BW, -BW);
    outer_wall_shape.emplace_back(-DL_sponge - BW, -BW);

    return outer_wall_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createValveShape(), ShapeBooleanOps::sub);
    }
};
/* Case-dependent geometries.*/
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int(valve_length / resolution_ref);
class ValveParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ValveParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real y = DH - (Real(i) + 0.5) * resolution_ref;
            initializePositionAndVolumetricMeasure(Vecd(valve_x, y), resolution_valve);
            Vecd normal_direction(1.0, 0.0);
            initializeSurfaceProperties(normal_direction, valve_thickness);
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
        if (base_particles_.pos_[index_i][1] > DH - base_length)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(1.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = 0.5 * u_ref_ * (1.0 - cos(2 * Pi * run_time / t_ref_));
        if (aligned_box_.checkInBounds(0, position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};
