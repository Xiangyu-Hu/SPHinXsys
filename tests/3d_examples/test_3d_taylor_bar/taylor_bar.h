/**
 * @file 	taylor_bar.h
 * @brief 	This is the case setup for plastic taylor bar.
 * @author 	Xiaojing Tang, Dong Wu and Xiangyu Hu
 */
#pragma once

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 0.0032; /**< X-direction domain. */
Real PW = 0.0324; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 5.0;
/** YOU can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real SL = particle_spacing_ref * 4.0; /**< Length of the holder is one layer particle. */
Real inner_circle_radius = PL;

Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; /**< Reference density. */
Real poisson = 0.35;  /**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;

class WallShape : public ComplexShape
{
  public:
    explicit WallShape(const std::string &shape_name) : ComplexShape(shape_name)
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
        Vecd translation_column(0.0, 0.0, 0.6 * PW);
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1.0), inner_circle_radius,
                                       0.5 * PW, resolution, translation_column);
    }
};

/**
 * application dependent initial condition
 */

class InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][2] = -227.0;
    }
};

/**
 * @class DynamicContactForceWithWall
 * @brief Computing the contact force with a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class DynamicContactForceWithWall : public LocalDynamics,
                                    public DataDelegateContact
{
  public:
    explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0)
        : LocalDynamics(solid_body_contact_relation.getSPHBody()),
          DataDelegateContact(solid_body_contact_relation),
          solid_(DynamicCast<Solid>(this, sph_body_->getBaseMaterial())),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          penalty_strength_(penalty_strength)
    {
        impedance_ = sqrt(solid_.ReferenceDensity() * solid_.ContactStiffness());
        reference_pressure_ = solid_.ContactStiffness();
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
    Solid &solid_;
    Real *Vol_;
    Vecd *vel_, *force_prior_; // note that prior force directly used here
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_n_;
    Real penalty_strength_;
    Real impedance_, reference_pressure_;
};
