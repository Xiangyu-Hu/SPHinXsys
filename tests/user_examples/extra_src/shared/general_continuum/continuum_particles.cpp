#include "continuum_particles.h"
#include "general_continuum.h"
namespace SPH
{
    ContinuumParticles::
        ContinuumParticles(SPHBody &sph_body, GeneralContinuum *continuum)
        : BaseParticles(sph_body, continuum), continuum_(*continuum) {}
    //=================================================================================================//
    void ContinuumParticles::initializeOtherVariables()
    {
        BaseParticles::initializeOtherVariables();

        registerVariable(acc_shear_, "AccelerationByShear");
        registerVariable(stress_tensor_, "StressTensor");
        registerVariable(stress_tensor_rate_, "StressTensorRate");
        registerVariable(shear_stress_, "ShearStress");
        registerVariable(shear_stress_rate_, "ShearStressRate");
        registerVariable(von_mises_stress_, "VonMisesStress");
        registerVariable(von_mises_strain_, "VonMisesStrain");
        registerVariable(velocity_gradient_, "VelocityGradient");
        registerVariable(strain_tensor_, "StrainTensor");
        registerVariable(strain_tensor_rate_, "StrainTensorRate");
        //----------------------------------------------------------------------
        //		register sortable particle data
        //----------------------------------------------------------------------
        registerSortableVariable<Vecd>("AccelerationByShear");
        registerSortableVariable<Matd>("StressTensor");
        registerSortableVariable<Matd>("StressTensorRate");
        registerSortableVariable<Matd>("ShearStress");
        registerSortableVariable<Matd>("ShearStressRate");
        registerSortableVariable<Real>("VonMisesStress");
        registerSortableVariable<Real>("VonMisesStrain");
        registerSortableVariable<Matd>("VelocityGradient");
        registerSortableVariable<Matd>("StrainTensor");
        registerSortableVariable<Matd>("StrainTensorRate");
        //----------------------------------------------------------------------
        registerVariable(pos0_, "InitialPosition", [&](size_t i) -> Vecd
                         { return pos_[i]; });
        registerVariable(n_, "NormalDirection");
        registerVariable(n0_, "InitialNormalDirection",
                         [&](size_t i) -> Vecd
                         { return n_[i]; });
    }
    //=================================================================================================//
    PlasticContinuumParticles::
        PlasticContinuumParticles(SPHBody& sph_body, PlasticContinuum* plastic_continuum)
        : ContinuumParticles(sph_body, plastic_continuum), plastic_continuum_(*plastic_continuum) {}
    //=================================================================================================//
    void PlasticContinuumParticles::initializeOtherVariables()
    {
        ContinuumParticles::initializeOtherVariables();

        registerVariable(elastic_strain_tensor_3D_, "ElasticStrainTensor3D");
        registerVariable(elastic_strain_rate_3D_, "ElasticStrainRate3D");
        registerVariable(strain_tensor_3D_, "StrainTensor3D");
        registerVariable(stress_tensor_3D_, "StressTensor3D");
        registerVariable(strain_rate_3D_, "StrainRate3D");
        registerVariable(stress_rate_3D_, "StressRate3D");
        registerVariable(shear_stress_3D_, "ShearStress3D");
        registerVariable(shear_strain_3D_, "ShearStrain3D");
        registerVariable(shear_stress_rate_3D_, "ShearStressRate3D");
        registerVariable(shear_strain_rate_3D_, "ShearStrainRate3D");
        registerVariable(vertical_stress_, "VerticalStress");
        registerVariable(acc_deviatoric_plastic_strain_, "AccDeviatoricPlasticStrain");
        //----------------------------------------------------------------------
        //		register sortable particle data
        //----------------------------------------------------------------------
        registerSortableVariable<Mat3d>("ElasticStrainTensor3D");
        registerSortableVariable<Mat3d>("ElasticStrainRate3D");
        registerSortableVariable<Mat3d>("StrainTensor3D");
        registerSortableVariable<Mat3d>("StressTensor3D");
        registerSortableVariable<Mat3d>("StrainRate3D");
        registerSortableVariable<Mat3d>("StressRate3D");

        registerSortableVariable<Mat3d>("ShearStress3D");
        registerSortableVariable<Mat3d>("ShearStrain3D");
        registerSortableVariable<Mat3d>("ShearStressRate3D");
        registerSortableVariable<Mat3d>("ShearStrainRate3D");

        registerSortableVariable<Real>("VerticalStress");
        registerSortableVariable<Real>("AccDeviatoricPlasticStrain");
    }
    //=================================================================================================//
    Real  PlasticContinuumParticles::getDeviatoricPlasticStrain(Mat3d& strain_tensor)
    {
        Mat3d deviatoric_strain_tensor = strain_tensor - (1 / (Real)Dimensions) * strain_tensor.trace() * Mat3d::Identity();
        Real sum = (deviatoric_strain_tensor.cwiseProduct(deviatoric_strain_tensor)).sum();
        return sqrt(sum * 2 / 3);
    }

} // namespace SPH