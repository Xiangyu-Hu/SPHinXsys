/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar using updated Lagragian SPH.
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu
 */
#pragma once
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 0.00391; /**< X-direction domain. */
Real PW = 0.02346; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 12.0;
Real SL = particle_spacing_ref * 4.0;
Real inner_circle_radius = PL;
Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 2700.0; /**< Reference density. */
Real poisson = 0.3;   /**< Poisson ratio. */
Real Youngs_modulus = 78.2e9;
Real yield_stress = 0.29e9;
Real vel_0 = 373.0;
Real U_max = vel_0;
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
/** Define the wall. */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
        Vecd translation_holder(0.0, 0.0, -0.5 * SL);
        add<TriangleMeshShapeBrick>(halfsize_holder, resolution, translation_holder);
    }
};
/** Define the body. */
class Column : public ComplexShape
{
  public:
    explicit Column(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(0.0, 0.0, 0.5 * PW + particle_spacing_ref);
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1.0), inner_circle_radius,
                                       0.5 * PW, resolution, translation_column);
    }
};
/**
 * application dependent initial condition
 */
class InitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][2] = -vel_0;
    }
};
/** Contact force. */
class DynamicContactForceWithWall : public LocalDynamics,
                                    public DataDelegateContact
{
  public:
    explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0)
        : LocalDynamics(solid_body_contact_relation.getSPHBody()),
          DataDelegateContact(solid_body_contact_relation),
          continuum_(DynamicCast<GeneralContinuum>(this, sph_body_->getBaseMaterial())),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          penalty_strength_(penalty_strength)
    {
        impedance_ = continuum_.ReferenceDensity() * sqrt(continuum_.ContactStiffness());
        reference_pressure_ = continuum_.ReferenceDensity() * continuum_.ContactStiffness();
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_n_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));
        }
    };
    virtual ~DynamicContactForceWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real particle_spacing_j1 = 1.0 / contact_bodies_[k]->getSPHAdaptation().ReferenceSpacing();
            Real particle_spacing_ratio2 =
                1.0 / (getSPHAdaptation().ReferenceSpacing() * particle_spacing_j1);
            particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

            Vecd *n_k = contact_n_[k];
            Vecd *vel_n_k = contact_vel_[k];
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Vecd n_k_j = n_k[index_j];
                Real impedance_p = 0.5 * impedance_ * (vel_[index_i] - vel_n_k[index_j]).dot(-n_k_j);
                Real overlap = contact_neighborhood.r_ij_[n] * n_k_j.dot(e_ij);
                Real delta = 2.0 * overlap * particle_spacing_j1;
                Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
                Real penalty_p = penalty_strength_ * beta * fabs(overlap) * reference_pressure_;
                // force due to pressure
                force -= 2.0 * (impedance_p + penalty_p) * e_ij.dot(n_k_j) *
                         n_k_j * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            }
        }
        force_prior_[index_i] += force * Vol_[index_i];
    };

  protected:
    GeneralContinuum &continuum_;
    Real *Vol_;
    Vecd *vel_, *force_prior_; // note that prior force directly used here
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_n_;
    Real penalty_strength_;
    Real impedance_, reference_pressure_;
};
