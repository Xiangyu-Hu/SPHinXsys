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

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 1.0;                                     /**< Tank length. */
Real DH = 0.7;                                     /**< Tank height. */
Real L_W = DL;                                     /**< water width. */
Real L_H = 0.15;                                   /**< water depth. */
Real Gate_x = 0.5 * L_W;                           /**< Width of the gate. */
Real Gate_width = 0.006;                           /**< Width of the gate. */
Real Gate_height = 0.18;                           /**< Height of the gate. */
Real particle_spacing_gate = Gate_width;           /**< Initial reference particle spacing of gate. */
Real particle_spacing_ref = particle_spacing_gate; /**< Initial reference particle spacing of water. */
Real BW = particle_spacing_ref * 4.0;              /**< Extending width for BCs. */
int particle_number_mid_surface = (0.026 + Gate_height + BW) / particle_spacing_gate;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// parameters for liquid sloshing Case:Xue Mian
Real a_0 = 0.01;              /**< amplitude of the sloshing. */
Real w_0 = 1.2 * 4.142814038; /**< frequency of the sloshing  0.583*8.9556. */
Real k_0 = a_0 * w_0 * w_0;   /**< parameter of sloshing x = k_0 * sin(2*pi*f_0*t). */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;                     /**< Reference density of fluid. */
Real gravity_g = 9.81;                    /**< Value of gravity. */
Real U_max = 2.0 * sqrt(gravity_g * L_H); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                  /**< Reference sound speed. */
Real mu_f = 1.0e-6;
/**
 * @brief Material properties of the elastic gate.
 */
Real rho0_s = 1250;           /**< Reference density of gate. */
Real Youngs_modulus = 30.0e6; /**< reference Youngs modulus. */
Real poisson = 0.47;          /**< Poisson ratio. */
Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * Gate_width;
/** create a water block shape */
double Gate_x_middle = (Gate_x - 0.0262 + Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008 + 0.0075) * 0.5;
StdVec<Vecd> water_block_shape{
    Vecd(0.0, 0.0), Vecd(0.0, L_H), Vecd(Gate_x_middle - 0.5 * particle_spacing_ref, L_H), Vecd(Gate_x_middle - 0.5 * particle_spacing_ref, 0.0), Vecd(0.0, 0.0)};
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
    Vecd(Gate_x - 0.0262, 0.0)};
/** Inner wall shape */
StdVec<Vecd> inner_wall_shape{
    Vecd(0.0, 0.0),
    Vecd(0.0, DH),
    Vecd(DL, DH),
    Vecd(DL, 0.0),
    Vecd(0.0, 0.0)};
/** Outer wall shape */
StdVec<Vecd> outer_wall_shape{
    Vecd(-BW, -BW),
    Vecd(-BW, DH + BW),
    Vecd(DL + BW, DH + BW),
    Vecd(DL + BW, -BW),
    Vecd(-BW, -BW)};
/* Case-dependent geometries.*/
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
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
            Real x = Gate_x_middle; // 0.5, 0.75, 1.5
            Real y = -BW + Real(i) * particle_spacing_gate;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_gate);
            Vec2d normal_direction = Vec2d(1.0, 0);
            initializeSurfaceProperties(normal_direction, Gate_width);
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
    VariableGravity() : Gravity(Vecd(0.0, -gravity_g)){};
    void UpdateAcceleration()
    {

        if (GlobalStaticVariables::physical_time_ <= 1.0)
        {
            global_acceleration_[0] = 0.0;
        }

        if (GlobalStaticVariables::physical_time_ > 1.0)
        {
            global_acceleration_[0] = k_0 * cos(w_0 * (GlobalStaticVariables::physical_time_ - 1.0));
        }
    }
};

/**
 * @brief Define the observer body.
 */
StdVec<Vecd> baffle_disp_probe_location{Vecd(Gate_x, 0.195), Vecd(Gate_x, 0.145), Vecd(Gate_x, 0.095), Vecd(Gate_x, 0.05)};
StdVec<Vecd> baffle_pressure_probe_location{Vecd(Gate_x - Gate_width / 2, 0.19), Vecd(Gate_x + Gate_width / 2, 0.19)};
StdVec<Vecd> fluid_pressure_probe_location{Vecd(DL, 0.045), Vecd(DL, 0.095), Vecd(DL, 0.165), Vecd(DL, 0.22)};