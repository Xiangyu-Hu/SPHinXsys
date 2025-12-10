/**
 * @file 	2d_eulerian_supersonic_flow_around_cylinder.h
 * @brief 	This is the compressible test for Eulerian supersonic flow around a cylinder
            with the FVM boundary algorithm in the SPHinXsys.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "eulerian_ghost_boundary.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real particle_spacing_ref = 1.0 / 7.0; /**< Initial reference particle spacing. */
Real calculation_circle_radius = 11.0;
Vecd insert_circle_center(7.0, 0.0);
Vecd calculation_circle_center(calculation_circle_radius, 0.0);
Real insert_circle_radius = 1.0;
Real calculation_circle_radius_with_BW = 11.0 + 4.0 * particle_spacing_ref;
BoundingBoxd system_domain_bounds(Vec2d(0.0, -calculation_circle_radius_with_BW), Vec2d(calculation_circle_radius_with_BW, calculation_circle_radius_with_BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho_reference = 1.0;
Real rho_farfield = 1.0;                     /**< initial density of one fluid. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */
Real p_farfield = 1.0 / heat_capacity_ratio; /**< initial pressure of one fluid. */
Real Mach_number = 2.0;
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> Creatrightsqureshape()
{
    // geometry
    std::vector<Vecd> right_square_shape;
    right_square_shape.push_back(Vecd(calculation_circle_center[0], -calculation_circle_radius));
    right_square_shape.push_back(Vecd(calculation_circle_center[0], calculation_circle_radius));
    right_square_shape.push_back(Vecd(calculation_circle_center[0] + calculation_circle_radius, calculation_circle_radius));
    right_square_shape.push_back(Vecd(calculation_circle_center[0] + calculation_circle_radius, -calculation_circle_radius));
    right_square_shape.push_back(Vecd(calculation_circle_center[0], -calculation_circle_radius));
    return right_square_shape;
}
class FluidBlock : public MultiPolygonShape
{
  public:
    explicit FluidBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(calculation_circle_center, calculation_circle_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(Creatrightsqureshape(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class SupersonicFlowInitialCondition : public fluid_dynamics::CompressibleFluidInitialCondition
{
  public:
    explicit SupersonicFlowInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::CompressibleFluidInitialCondition(sph_body){};
    void update(size_t index_i, Real dt)
    {
        p_[index_i] = p_farfield;
        rho_[index_i] = rho_farfield;
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        Real sound_speed = sqrt(p_[index_i] * gamma_ / rho_[index_i]);
        vel_[index_i][0] = Mach_number * sound_speed;
        vel_[index_i][1] = 0.0;
        mom_[index_i] = mass_[index_i] * vel_[index_i];
        Real rho_e = p_[index_i] / (gamma_ - 1.0);
        E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    Real gamma_ = heat_capacity_ratio;
};

//----------------------------------------------------------------------
//	SupersonicFlowBoundaryConditionSetup
//----------------------------------------------------------------------
class SupersonicFlowBoundaryConditionSetup : public GhostBoundaryConditionSetupInESPH
{
  public:
    SupersonicFlowBoundaryConditionSetup(BaseInnerRelation &inner_relation, GhostCreationInESPH &ghost_creation)
        : GhostBoundaryConditionSetupInESPH(inner_relation, ghost_creation),
          p_(particles_->getVariableDataByName<Real>("Pressure")),
          E_(particles_->getVariableDataByName<Real>("TotalEnergy"))
    {
        setupBoundaryTypes();
    };
    virtual ~SupersonicFlowBoundaryConditionSetup(){};

    // Here, we follow the FVM labelling different boundary conditions with different number
    void setupBoundaryTypes() override
    {
        for (size_t ghost_number = 0; ghost_number != real_and_ghost_particle_data_.size(); ++ghost_number)
        {
            size_t index_i = real_and_ghost_particle_data_[ghost_number].real_index_;
            // around cylinder
            if ((pos_[index_i] - insert_circle_center).norm() <= insert_circle_radius + 5.0 * particle_spacing_ref)
            {
                boundary_type_[index_i] = 3;
            }
            // outer boundary
            else if ((pos_[index_i] - insert_circle_center).norm() > insert_circle_radius + 5.0 * particle_spacing_ref)
            {
                boundary_type_[index_i] = 9;
            }
        }
    };

    // For inner boundaries, we set the refective boundary conditions as boundary_type 3.
    void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
    {
        rho_[ghost_index] = rho_[index_i];
        p_[ghost_index] = p_[index_i];
        Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij)) + (-e_ij.dot(vel_[index_i]) * (e_ij));
        mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
        mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
        E_[ghost_index] = rho_e * Vol_[ghost_index] + 0.5 * mass_[ghost_index] * vel_[ghost_index].squaredNorm();
    }

    // For outer boundaries, we set the non-refective far-field boundary conditions as boundary_type 9.
    void applyFarFieldBoundary(size_t ghost_index, size_t index_i) override
    {
        Vecd normal_direction_index_i = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd velocity_farfield = Vecd::Zero();
        Real sound_speed = sqrt(p_farfield * heat_capacity_ratio / rho_farfield);
        velocity_farfield[0] = Mach_number * sound_speed;
        Real velocity_farfield_normal = velocity_farfield.dot(normal_direction_index_i);
        Real velocity_boundary_normal = vel_[index_i].dot(normal_direction_index_i);

        // judge the inflow boundary condition
        if (normal_direction_index_i[0] <= 0.0 || fabs(normal_direction_index_i[1]) > fabs(normal_direction_index_i[0]))
        {
            // supersonic inflow condition
            if (fabs(velocity_boundary_normal) >= sound_speed)
            {
                vel_[ghost_index] = velocity_farfield;
                p_[ghost_index] = p_farfield;
                rho_[ghost_index] = rho_farfield;
                mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
                mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
                Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
                E_[ghost_index] = rho_e * Vol_[ghost_index] + 0.5 * mass_[ghost_index] * vel_[ghost_index].squaredNorm();
            }
            // subsonic inflow condition
            if (fabs(velocity_boundary_normal) < sound_speed)
            {
                Real inner_weight_summation = W0_ * Vol_[index_i];
                Real rho_summation = 0.0;
                Real p_summation = 0.0;
                Real vel_normal_summation = 0.0;
                size_t total_inner_neighbor_particles = 0;
                Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real W_ij = inner_neighborhood.W_ij_[n];
                    inner_weight_summation += W_ij * Vol_[index_j];
                    rho_summation += rho_[index_j];
                    vel_normal_summation += vel_[index_j].dot(normal_direction_index_i);
                    p_summation += p_[index_j];
                    total_inner_neighbor_particles += 1;
                }
                Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
                Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
                Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

                p_[ghost_index] = p_average * inner_weight_summation + p_farfield * (1.0 - inner_weight_summation);
                rho_[ghost_index] = rho_average * inner_weight_summation + rho_farfield * (1.0 - inner_weight_summation);
                Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
                vel_[ghost_index] = vel_normal * normal_direction_index_i + (velocity_farfield - velocity_farfield_normal * normal_direction_index_i);
                mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
                mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
                Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
                E_[ghost_index] = rho_e * Vol_[ghost_index] + 0.5 * mass_[ghost_index] * vel_[ghost_index].squaredNorm();
            }
        }
        else
        {
            // supersonic outflow condition
            if (fabs(velocity_boundary_normal) >= sound_speed)
            {
                vel_[ghost_index] = vel_[index_i];
                p_[ghost_index] = p_[index_i];
                rho_[ghost_index] = rho_[index_i];
                mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
                mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
                Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
                E_[ghost_index] = rho_e * Vol_[ghost_index] + 0.5 * mass_[ghost_index] * vel_[ghost_index].squaredNorm();
            }
            // subsonic outflow condition
            if (fabs(velocity_boundary_normal) < sound_speed)
            {
                Real inner_weight_summation = W0_ * Vol_[index_i];
                Real rho_summation = 0.0;
                Real p_summation = 0.0;
                Real vel_normal_summation(0.0);
                Vecd vel_tangential_summation = Vecd::Zero();
                size_t total_inner_neighbor_particles = 0;
                Neighborhood &inner_neighborhood = inner_configuration_[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real W_ij = inner_neighborhood.W_ij_[n];
                    inner_weight_summation += W_ij * Vol_[index_j];
                    rho_summation += rho_[index_j];
                    vel_normal_summation += vel_[index_j].dot(normal_direction_index_i);
                    vel_tangential_summation += vel_[index_j] - vel_[index_j].dot(normal_direction_index_i) * normal_direction_index_i;
                    p_summation += p_[index_j];
                    total_inner_neighbor_particles += 1;
                }
                Real rho_average = rho_summation / (total_inner_neighbor_particles + TinyReal);
                Real vel_normal_average = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
                Vecd vel_tangential_average = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);
                Real p_average = p_summation / (total_inner_neighbor_particles + TinyReal);

                p_[ghost_index] = p_average * inner_weight_summation + p_farfield * (1.0 - inner_weight_summation);
                rho_[ghost_index] = rho_average * inner_weight_summation + rho_farfield * (1.0 - inner_weight_summation);
                Real vel_normal = vel_normal_average * inner_weight_summation + velocity_farfield_normal * (1.0 - inner_weight_summation);
                vel_[ghost_index] = vel_normal * normal_direction_index_i + vel_tangential_average;
                mass_[ghost_index] = rho_[ghost_index] * Vol_[ghost_index];
                mom_[ghost_index] = mass_[ghost_index] * vel_[ghost_index];
                Real rho_e = p_[ghost_index] / (heat_capacity_ratio - 1.0);
                E_[ghost_index] = rho_e * Vol_[ghost_index] + 0.5 * mass_[ghost_index] * vel_[ghost_index].squaredNorm();
            }
        }
    }

  protected:
    Real *p_, *E_;
};
