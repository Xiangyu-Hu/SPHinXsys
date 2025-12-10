/**
 * @file 	diffusion_RobinBC.h
 * @brief 	This is the head files used by diffusion_RobinBC.cpp.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#ifndef DIFFUSION_ROBIN_BC_H
#define DIFFUSION_ROBIN_BC_H

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 100.0;
Real BW = resolution_ref * 2.0;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
std::string diffusion_species_name = "Phi";
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 100.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real convection = 100.0;
Real T_infinity = 400.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
    std::vector<Vecd> thermalDomainShape;
    thermalDomainShape.push_back(Vecd(0.0, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, H));
    thermalDomainShape.push_back(Vecd(L, H));
    thermalDomainShape.push_back(Vecd(L, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, 0.0));

    return thermalDomainShape;
}

std::vector<Vecd> left_temperature_region{
    Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW), Vecd(0.4 * L, H + BW),
    Vecd(0.4 * L, H), Vecd(0.3 * L, H)};

std::vector<Vecd> right_temperature_region{
    Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW), Vecd(0.7 * L, H + BW),
    Vecd(0.7 * L, H), Vecd(0.6 * L, H)};

std::vector<Vecd> convection_region{
    Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0), Vecd(0.55 * L, 0),
    Vecd(0.55 * L, -BW), Vecd(0.45 * L, -BW)};

//----------------------------------------------------------------------
// Define extra classes which are used in the main program.
// These classes are defined under the namespace of SPH.
//----------------------------------------------------------------------
namespace SPH
{
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
    }
};

class DirichletWallBoundary : public MultiPolygonShape
{
  public:
    explicit DirichletWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(left_temperature_region, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(right_temperature_region, ShapeBooleanOps::add);
    }
};

class RobinWallBoundary : public MultiPolygonShape
{
  public:
    explicit RobinWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(convection_region, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature;
    };

  protected:
    Real *phi_;
};

class DirichletWallBoundaryInitialCondition : public LocalDynamics
{
  public:
    explicit DirichletWallBoundaryInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;

        if (pos_[index_i][1] > H && pos_[index_i][0] > 0.3 * L && pos_[index_i][0] < 0.4 * L)
        {
            phi_[index_i] = left_temperature;
        }
        if (pos_[index_i][1] > H && pos_[index_i][0] > 0.6 * L && pos_[index_i][0] < 0.7 * L)
        {
            phi_[index_i] = right_temperature;
        }
    }

  protected:
    Vecd *pos_;
    Real *phi_;
};

class RobinBoundaryDefinition : public LocalDynamics
{
  public:
    explicit RobinBoundaryDefinition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)),
          phi_convection_(particles_->template getVariableDataByName<Real>(diffusion_species_name + "Convection")),
          phi_infinity_(*(this->particles_->template getSingularVariableByName<Real>(diffusion_species_name + "Infinity")->Data())) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;
        if (pos_[index_i][1] < 0 && pos_[index_i][0] > 0.45 * L && pos_[index_i][0] < 0.55 * L)
        {
            phi_convection_[index_i] = convection;
            phi_infinity_ = T_infinity;
        }
    }

  protected:
    Vecd *pos_;
    Real *phi_, *phi_convection_;
    Real &phi_infinity_;
};

using DiffusionBodyRelaxation = DiffusionBodyRelaxationComplex<
    BaseDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet, Robin>;

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the middle line. */
    size_t number_of_observation_points = 5;
    Real range_of_measure = L;
    Real start_of_measure = 0;

    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(
            0.5 * L, range_of_measure * Real(i) / Real(number_of_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
} // namespace SPH
#endif // DIFFUSION_ROBIN_BC_H