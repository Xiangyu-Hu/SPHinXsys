#include "sphinxsys.h"
using namespace SPH;

//------------Coupling class----------------------
class CouplingPart : public BodyPartByParticle,
                     public DataDelegateContact
{
  public:
    explicit CouplingPart(BaseContactRelation &contact_relation)
        : BodyPartByParticle(contact_relation.getSPHBody(), "CoupledParticles"),
          DataDelegateContact(contact_relation)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&CouplingPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
        std::cout << "Number of coupled particles: " << body_part_particles_.size() << std::endl;
    };
    void tagManually(size_t index_i)
    {
        if (std::any_of(contact_configuration_.begin(), contact_configuration_.end(), [&](auto &k)
                        { return (*k)[index_i].current_size_ > 0; }))
        {
            body_part_particles_.emplace_back(index_i);
        }
    }
    void reset_ids()
    {
        body_part_particles_.clear();
        TaggingParticleMethod tagging_particle_method = std::bind(&CouplingPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
        std::cout << "Number of coupled particles: " << body_part_particles_.size() << std::endl;
    }
};

class SolidRigidCouplingStressFirst : public BaseLocalDynamics<BodyPartByParticle>,
                                      public DataDelegateContact
{
  private:
    Real rho0_;
    Real *mass_;
    Vecd *vel_;
    Matd *stress_PK1_B_;
    StdVec<Real *> contact_Vol_;
    StdVec<Matd *> contact_stress_PK1_B_;

    Vecd *coupling_force_;

  public:
    explicit SolidRigidCouplingStressFirst(CouplingPart &coupling_part)
        : BaseLocalDynamics<BodyPartByParticle>(coupling_part),
          DataDelegateContact(coupling_part.getBodyRelation()),
          rho0_(particles_->getBaseMaterial().ReferenceDensity()),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          vel_(particles_->registerStateVariable<Vecd>("Velocity")),
          stress_PK1_B_(particles_->registerStateVariable<Matd>("StressPK1OnParticle")),
          coupling_force_(particles_->registerStateVariable<Vecd>("CouplingForce"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_stress_PK1_B_.push_back(contact_particles_[k]->registerStateVariable<Matd>("StressPK1OnParticle"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            const Matd *stress_PK1_B_k = contact_stress_PK1_B_[k];
            const auto *Vol_k = contact_Vol_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                // In Integration1stHalf, a pair-wise dissipation term is added to the stress
                // However, tests have shown that this term leads to stress oscillations in solid-solid coupling
                // The cause is not yet clear, but the term is omitted here to avoid the oscillations
                force += mass_[index_i] / rho0_ * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] *
                         (stress_PK1_B_[index_i] + stress_PK1_B_k[index_j]) * e_ij;
            }
        }
        coupling_force_[index_i] = force;
    }

    inline void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] += coupling_force_[index_i] * dt / mass_[index_i];
    }
};

class SolidRigidCouplingStressSecond : public BaseLocalDynamics<BodyPartByParticle>,
                                       public DataDelegateContact
{
  private:
    Vecd *vel_;
    Matd *B_;
    Matd *F_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_Vol_;

    Matd *coupling_dF_dt_;

  public:
    explicit SolidRigidCouplingStressSecond(CouplingPart &coupling_part)
        : BaseLocalDynamics<BodyPartByParticle>(coupling_part),
          DataDelegateContact(coupling_part.getBodyRelation()),
          vel_(particles_->registerStateVariable<Vecd>("Velocity")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          F_(particles_->registerStateVariable<Matd>("DeformationGradient")),
          coupling_dF_dt_(particles_->registerStateVariable<Matd>("CouplingDeformationGradientChangeRate"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->registerStateVariable<Vecd>("Velocity"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Matd deformation_gradient_change_rate = Matd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            const auto *Vol_k = contact_Vol_[k];
            const auto *vel_k = contact_vel_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];

                SPH::Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                deformation_gradient_change_rate -= (vel_[index_i] - vel_k[index_j]) * gradW_ij.transpose();
            }
        }
        deformation_gradient_change_rate *= B_[index_i];
        coupling_dF_dt_[index_i] = deformation_gradient_change_rate;
    }

    inline void update(size_t index_i, Real dt = 0.0)
    {
        F_[index_i] += 0.5 * coupling_dF_dt_[index_i] * dt;
    }
};

class RigidSolidCoupling : public BaseForcePrior<BodyPartByParticle>,
                           public DataDelegateContact
{
  private:
    Real rho0_;
    Real *mass_;
    Matd *stress_PK1_B_;
    StdVec<Real *> contact_Vol_;
    StdVec<Matd *> contact_stress_PK1_B_;

    Vecd *coupling_force_;

  public:
    explicit RigidSolidCoupling(CouplingPart &coupling_part)
        : BaseForcePrior<BodyPartByParticle>(coupling_part, "CouplingForce"),
          DataDelegateContact(coupling_part.getBodyRelation()),
          rho0_(particles_->getBaseMaterial().ReferenceDensity()),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          stress_PK1_B_(particles_->registerStateVariable<Matd>("StressPK1OnParticle")),
          coupling_force_(particles_->registerStateVariable<Vecd>("CouplingForce"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_stress_PK1_B_.push_back(contact_particles_[k]->registerStateVariable<Matd>("StressPK1OnParticle"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            const auto *stress_PK1_B_k = contact_stress_PK1_B_[k];
            const auto *Vol_k = contact_Vol_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                // In Integration1stHalf, a pair-wise dissipation term is added to the stress
                // However, tests have shown that this term leads to stress oscillations in solid-solid coupling
                // The cause is not yet clear, but the term is omitted here to avoid the oscillations
                force += mass_[index_i] / rho0_ * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] *
                         (stress_PK1_B_[index_i] + stress_PK1_B_k[index_j]) * e_ij;
            }
        }
        coupling_force_[index_i] = force;
    }
};

class RigidBodyPseudoDeformationGradient : public BaseLocalDynamics<BodyPartByParticle>,
                                           public DataDelegateInner,
                                           public DataDelegateContact
{
  private:
    Vecd *vel_;
    Real *Vol_;
    Matd *B_;
    Matd *F_;
    Matd *dF_dt_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_Vol_;

  public:
    explicit RigidBodyPseudoDeformationGradient(CouplingPart &coupling_part, BaseInnerRelation &inner_relation)
        : BaseLocalDynamics<BodyPartByParticle>(coupling_part),
          DataDelegateInner(inner_relation),
          DataDelegateContact(coupling_part.getBodyRelation()),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          F_(particles_->registerStateVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
          dF_dt_(particles_->registerStateVariable<Matd>("DeformationRate"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Velocity"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &vel_n_i = vel_[index_i];

        Matd deformation_gradient_change_rate = Matd::Zero();

        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            deformation_gradient_change_rate -= (vel_n_i - vel_[index_j]) * gradW_ijV_j.transpose();
        }

        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            const auto *Vol_k = contact_Vol_[k];
            const auto *vel_k = contact_vel_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];

                Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                deformation_gradient_change_rate -= (vel_n_i - vel_k[index_j]) * gradW_ijV_j.transpose();
            }
        }

        dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
    };

    inline void update(size_t index_i, Real dt = 0.0)
    {
        F_[index_i] += 0.5 * dF_dt_[index_i] * dt;
    }
};

class RigidBodyPseudoPK1Stress : public BaseLocalDynamics<BodyPartByParticle>
{
  private:
    ElasticSolid &elastic_solid_;
    Matd *stress_PK1_B_;
    Matd *B_;
    Matd *F_;
    Matd *dF_dt_;

  public:
    explicit RigidBodyPseudoPK1Stress(CouplingPart &coupling_part, ElasticSolid &elastic_solid)
        : BaseLocalDynamics<BodyPartByParticle>(coupling_part),
          elastic_solid_(elastic_solid),
          stress_PK1_B_(particles_->registerStateVariable<Matd>("StressPK1OnParticle")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          F_(particles_->registerStateVariable<Matd>("DeformationGradient")),
          dF_dt_(particles_->registerStateVariable<Matd>("DeformationRate")) {}

    inline void update(size_t index_i, Real dt = 0.0)
    {
        F_[index_i] += dF_dt_[index_i] * 0.5 * dt;
        stress_PK1_B_[index_i] = elastic_solid_.StressPK1(F_[index_i], index_i) * B_[index_i].transpose();
    }
};
//------------function to setup this test----------------------
BoundingBox union_bounding_box(const BoundingBox &a, const BoundingBox &b)
{
    BoundingBox out = a;
    out.first_[0] = std::min(a.first_[0], b.first_[0]);
    out.first_[1] = std::min(a.first_[1], b.first_[1]);
    out.first_[2] = std::min(a.first_[2], b.first_[2]);
    out.second_[0] = std::max(a.second_[0], b.second_[0]);
    out.second_[1] = std::max(a.second_[1], b.second_[1]);
    out.second_[2] = std::max(a.second_[2], b.second_[2]);
    return out;
}

class FixPart : public BodyPartByParticle
{
  private:
    std::function<bool(Vec3d &)> contains_;

  public:
    FixPart(SPHBody &body, const std::string &body_part_name, std::function<bool(Vec3d &)> contains)
        : BodyPartByParticle(body, body_part_name),
          contains_(std::move(contains))
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (contains_(pos_[index_i]))
            body_part_particles_.push_back(index_i);
    };
};

Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    // the physical viscosity is defined in the paper of prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant (0.4 for Beam)
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}