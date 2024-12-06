#include "elastic_solid.h"
#include "base_particles.hpp"

#ifdef max
#undef max
#endif

namespace SPH
{
//=================================================================================================//
void ElasticSolid::setSoundSpeeds()
{
    c0_ = sqrt(K0_ / rho0_);
    ct0_ = sqrt(E0_ / rho0_);
    cs0_ = sqrt(G0_ / rho0_);
};
//=================================================================================================//
Real ElasticSolid::PairNumericalDamping(Real dE_dt_ij, Real smoothing_length)
{
    return 0.5 * rho0_ * c0_ * dE_dt_ij * smoothing_length;
}
//=================================================================================================//
Matd ElasticSolid::DeviatoricKirchhoff(const Matd &deviatoric_be)
{
    return G0_ * deviatoric_be;
}
//=================================================================================================//
Vecd *ElasticSolid::AverageVelocity(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("AverageVelocity");
}
//=================================================================================================//
Vecd *ElasticSolid::AverageAcceleration(BaseParticles *base_particles)
{
    return base_particles->registerStateVariable<Vecd>("AverageAcceleration");
}
//=================================================================================================//
DiscreteVariable<Vecd> *ElasticSolid::AverageVelocityVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableOnly<Vecd>("AverageVelocity");
}
//=================================================================================================//
DiscreteVariable<Vecd> *ElasticSolid::AverageAccelerationVariable(BaseParticles *base_particles)
{
    return base_particles->registerStateVariableOnly<Vecd>("AverageAcceleration");
}
//=================================================================================================//
LinearElasticSolid::
    LinearElasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio) : ElasticSolid(rho0)
{
    material_type_name_ = "LinearElasticSolid";
    E0_ = youngs_modulus;
    nu_ = poisson_ratio;
    G0_ = getShearModulus(youngs_modulus, poisson_ratio);
    K0_ = getBulkModulus(youngs_modulus, poisson_ratio);
    lambda0_ = getLambda(youngs_modulus, poisson_ratio);
    setSoundSpeeds();
    setContactStiffness(c0_);
}
//=================================================================================================//
Real LinearElasticSolid::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
{
    return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Real LinearElasticSolid::getShearModulus(Real youngs_modulus, Real poisson_ratio)
{
    return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
}
//=================================================================================================//
Real LinearElasticSolid::getLambda(Real youngs_modulus, Real poisson_ratio)
{
    return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Matd LinearElasticSolid::StressPK1(Matd &F, size_t index_i)
{
    return F * StressPK2(F, index_i);
}
//=================================================================================================//
Matd LinearElasticSolid::StressPK2(Matd &F, size_t index_i)
{
    Matd strain = 0.5 * (F.transpose() + F) - Matd::Identity();
    return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain;
}
//=================================================================================================//
Matd LinearElasticSolid::StressCauchy(Matd &almansi_strain, size_t index_i)
{
    return lambda0_ * almansi_strain.trace() * Matd::Identity() + 2.0 * G0_ * almansi_strain;
}
//=================================================================================================//
Real LinearElasticSolid::VolumetricKirchhoff(Real J)
{
    return K0_ * J * (J - 1);
}
//=================================================================================================//
Matd SaintVenantKirchhoffSolid::StressPK2(Matd &F, size_t index_i)
{
    Matd strain = 0.5 * (F.transpose() * F - Matd::Identity());
    return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain;
}
//=================================================================================================//
Matd NeoHookeanSolid::StressPK2(Matd &F, size_t index_i)
{
    // This formulation allows negative determinant of F. Please refer Eq. (12) in
    // Smith et al. (2018) Stable Neo-Hookean Flesh Simulation.
    // ACM Transactions on Graphics, Vol. 37, No. 2, Article 12.
    Matd right_cauchy = F.transpose() * F;
    Real J = F.determinant();
    return G0_ * Matd::Identity() + (lambda0_ * (J - 1.0) - G0_) * J * right_cauchy.inverse();
}
//=================================================================================================//
Matd NeoHookeanSolid::StressCauchy(Matd &almansi_strain, size_t index_i)
{
    Matd B = (-2.0 * almansi_strain + Matd::Identity()).inverse();
    Real J = sqrt(B.determinant());
    Matd cauchy_stress = 0.5 * K0_ * (J - 1.0 / J) * Matd::Identity() +
                         G0_ * pow(J, -2.0 * OneOverDimensions - 1.0) *
                             (B - OneOverDimensions * B.trace() * Matd::Identity());
    return cauchy_stress;
}
//=================================================================================================//
Real NeoHookeanSolid::VolumetricKirchhoff(Real J)
{
    return 0.5 * K0_ * (J * J - 1);
}
//=================================================================================================//
Matd NeoHookeanSolidIncompressible::StressPK2(Matd &F, size_t index_i)
{
    Matd right_cauchy = F.transpose() * F;
    Real I_1 = right_cauchy.trace();       // first strain invariant
    Real I_3 = right_cauchy.determinant(); // first strain invariant
    return G0_ * pow(I_3, -1.0 / 3.0) * (Matd::Identity() - 1.0 / 3.0 * I_1 * right_cauchy.inverse());
}
//=================================================================================================//
Matd NeoHookeanSolidIncompressible::
    StressCauchy(Matd &almansi_strain, size_t index_i)
{
    // TODO: implement
    return {};
}
//=================================================================================================//
Real NeoHookeanSolidIncompressible::VolumetricKirchhoff(Real J)
{
    return 0.5 * K0_ * (J * J - 1);
}
//=================================================================================================//
OrthotropicSolid::OrthotropicSolid(Real rho_0, std::array<Vecd, Dimensions> a, std::array<Real, Dimensions> E,
                                   std::array<Real, Dimensions> G, std::array<Real, Dimensions> poisson)
    // since Poisson ratio can be larger than 0.5 for orthotropic materials,
    // the sound speed needs to be reset after parameters are computed
    : LinearElasticSolid(rho_0, *std::max_element(E.cbegin(), E.cend()),
                         *std::max_element(poisson.cbegin(), poisson.cend())),
      a_(a), E_(E), G_(G), poisson_(poisson)
{
    // parameters for derived class
    material_type_name_ = "OrthotropicSolid";
    CalculateA0();
    CalculateAllMu();
    CalculateAllLambda();
    // we take the max. K in three pinciple directions to approximate the maximum of
    // the Bulk modulus --> for time step size calculation
    reset_sound_speeds();
};
//=================================================================================================//
Matd OrthotropicSolid::StressPK2(Matd &F, size_t index_i)
{
    Matd strain = 0.5 * (F.transpose() * F - Matd::Identity());
    Matd stress_PK2 = Matd::Zero();
    for (int i = 0; i < Dimensions; i++)
    {
        // outer sum (a{1-3})
        Matd Summa2 = Matd::Zero();
        for (int j = 0; j < Dimensions; j++)
        {
            // inner sum (b{1-3})
            Summa2 += Lambda_(i, j) * (CalculateBiDotProduct(A_[i], strain) * A_[j] +
                                       CalculateBiDotProduct(A_[j], strain) * A_[i]);
        }
        stress_PK2 += Mu_[i] * (A_[i] * strain + strain * A_[i]) + 0.5 * Summa2;
    }
    return stress_PK2;
}
//=================================================================================================//
Matd StressCauchy(Matd &almansi_strain, size_t particle_index_i)
{
    // TODO: implement
    throw(std::runtime_error("Cauchy stress of orthotropic material not yet implemented"));
}
//=================================================================================================//
Real OrthotropicSolid::VolumetricKirchhoff(Real J)
{
    return K0_ * J * (J - 1);
}
//=================================================================================================//
void OrthotropicSolid::CalculateA0()
{
    for (int i = 0; i < Dimensions; ++i)
        A_[i] = a_[i] * a_[i].transpose();
}
//=================================================================================================//
void OrthotropicSolid::reset_sound_speeds()
{
    std::array<Real, Dimensions> K;
    for (int i = 0; i < Dimensions; ++i)
        K[i] = (Lambda_.row(i).sum() + 2.0 * Mu_[i]) / 3.0;
    K0_ = *std::max_element(K.cbegin(), K.cend());
    G0_ = *std::max_element(G_.cbegin(), G_.cend());
    lambda0_ = Lambda_.maxCoeff();
    setSoundSpeeds();
    setContactStiffness(c0_);
}
//=================================================================================================//
Matd FeneNeoHookeanSolid::StressPK2(Matd &F, size_t index_i)
{
    Matd right_cauchy = F.transpose() * F;
    Matd strain = 0.5 * (right_cauchy - Matd::Identity());
    Real J = F.determinant();
    return G0_ / (1.0 - 2.0 * strain.trace() / j1_m_) * Matd::Identity() +
           (lambda0_ * (J - 1.0) - G0_) * J * right_cauchy.inverse();
}
//=================================================================================================//
Real Muscle::getShearModulus(const Real (&a0)[4], const Real (&b0)[4])
{
    // This is only the background material property.
    // The previous version seems not correct because it leads to
    // that shear modulus is even bigger than bulk modulus.
    return a0[0];
}
//=================================================================================================//
Real Muscle::getPoissonRatio(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4])
{
    Real shear_modulus = getShearModulus(a0, b0);
    return 0.5 * (3.0 * bulk_modulus - 2.0 * shear_modulus) /
           (3.0 * bulk_modulus + shear_modulus);
}
//=================================================================================================//
Real Muscle::getYoungsModulus(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4])
{
    return 3.0 * bulk_modulus * (1.0 - 2.0 * getPoissonRatio(bulk_modulus, a0, b0));
}
//=================================================================================================//
Matd Muscle::StressPK2(Matd &F, size_t i)
{
    Matd right_cauchy = F.transpose() * F;
    Real I_ff_1 = (right_cauchy * f0_).transpose() * f0_ - 1.0;
    Real I_ss_1 = (right_cauchy * s0_).transpose() * s0_ - 1.0;
    Real I_fs = (right_cauchy * f0_).transpose() * s0_;
    Real I_1_1 = right_cauchy.trace() - Real(Dimensions);
    Real J = F.determinant();
    return a0_[0] * exp(b0_[0] * I_1_1) * Matd::Identity() +
           (lambda0_ * (J - 1.0) - a0_[0]) * J * right_cauchy.inverse() +
           2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * f0f0_ +
           2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * s0s0_ +
           a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * f0s0_;
}
//=================================================================================================//
Real Muscle::VolumetricKirchhoff(Real J)
{
    return K0_ * J * (J - 1);
}
//=================================================================================================//
Matd LocallyOrthotropicMuscle::StressPK2(Matd &F, size_t i)
{
    Matd right_cauchy = F.transpose() * F;
    Real I_ff_1 = (right_cauchy * local_f0_[i]).transpose() * local_f0_[i] - 1.0;
    Real I_ss_1 = (right_cauchy * local_s0_[i]).transpose() * local_s0_[i] - 1.0;
    Real I_fs = (right_cauchy * local_f0_[i]).transpose() * local_s0_[i];
    Real I_1_1 = right_cauchy.trace() - Real(Dimensions);
    Real J = F.determinant();
    return a0_[0] * exp(b0_[0] * I_1_1) * Matd::Identity() +
           (lambda0_ * (J - 1.0) - a0_[0]) * J * right_cauchy.inverse() +
           2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * local_f0f0_[i] +
           2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * local_s0s0_[i] +
           a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * local_f0s0_[i];
}
//=================================================================================================//
void LocallyOrthotropicMuscle::registerLocalParameters(BaseParticles *base_particles)
{
    Muscle::registerLocalParameters(base_particles);
    local_f0_ = base_particles->registerStateVariable<Vecd>("Fiber");
    local_s0_ = base_particles->registerStateVariable<Vecd>("Sheet");
}
void LocallyOrthotropicMuscle::registerLocalParametersFromReload(BaseParticles *base_particles)
{
    Muscle::registerLocalParametersFromReload(base_particles);
    local_f0_ = base_particles->registerStateVariableFromReload<Vecd>("Fiber");
    local_s0_ = base_particles->registerStateVariableFromReload<Vecd>("Sheet");
}
//=================================================================================================//
void LocallyOrthotropicMuscle::initializeLocalParameters(BaseParticles *base_particles)
{
    Muscle::initializeLocalParameters(base_particles);
    local_f0f0_ = base_particles->registerStateVariable<Matd>(
        "FiberFiberTensor", [&](size_t i) -> Matd
        { return local_f0_[i] * local_f0_[i].transpose(); });
    local_s0s0_ = base_particles->registerStateVariable<Matd>(
        "SheetSheetTensor", [&](size_t i) -> Matd
        { return local_s0_[i] * local_s0_[i].transpose(); });
    local_f0s0_ = base_particles->registerStateVariable<Matd>(
        "FiberSheetTensor", [&](size_t i) -> Matd
        { return local_f0_[i] * local_s0_[i].transpose() +
                 local_s0_[i] * local_f0_[i].transpose(); });
}
//=================================================================================================//
} // namespace SPH
