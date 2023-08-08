#ifndef CONTINUUM_PARTICLES_H
#define CONTINUUM_PARTICLES_H
#include "base_particles.hpp"
#include "general_continuum.h"
namespace SPH
{
    class GeneralContinuum;

    class ContinuumParticles : public BaseParticles
    {
      public:
        StdLargeVec<Matd> strain_tensor_;
        StdLargeVec<Matd> strain_tensor_rate_;
        StdLargeVec<Vecd> acc_shear_;

        StdLargeVec<Matd> stress_tensor_;
        StdLargeVec<Matd> stress_tensor_rate_;

        StdLargeVec<Matd> shear_stress_;
        StdLargeVec<Matd> shear_stress_rate_;
        StdLargeVec<Matd> velocity_gradient_;

        StdLargeVec<Real> von_mises_stress_;
        StdLargeVec<Real> von_mises_strain_;

        StdLargeVec<Vecd> pos0_; /**< initial position */
        StdLargeVec<Vecd> n_;    /**<  current normal direction */
        StdLargeVec<Vecd> n0_;   /**<  initial normal direction */

        GeneralContinuum &continuum_;

        ContinuumParticles(SPHBody &sph_body, GeneralContinuum *continuum);
        virtual ~ContinuumParticles(){};

        virtual void initializeOtherVariables() override;
        virtual ContinuumParticles *ThisObjectPtr() override { return this; };
    };

    class PlasticContinuumParticles : public ContinuumParticles
    {
    public:
        StdLargeVec<Mat3d> elastic_strain_tensor_3D_;
        StdLargeVec<Mat3d> elastic_strain_rate_3D_;

        StdLargeVec<Mat3d> strain_tensor_3D_;
        StdLargeVec<Mat3d> stress_tensor_3D_;
        StdLargeVec<Mat3d> strain_rate_3D_;
        StdLargeVec<Mat3d> stress_rate_3D_;

        StdLargeVec<Mat3d> shear_stress_3D_;
        StdLargeVec<Mat3d> shear_strain_3D_;
        StdLargeVec<Mat3d> shear_stress_rate_3D_;
        StdLargeVec<Mat3d> shear_strain_rate_3D_;

        StdLargeVec<Real> vertical_stress_;
        StdLargeVec<Real> acc_deviatoric_plastic_strain_;

        Real getDeviatoricPlasticStrain(Mat3d& strain_tensor);

        PlasticContinuum& plastic_continuum_;

        PlasticContinuumParticles(SPHBody& sph_body, PlasticContinuum* plastic_continuum);
        virtual ~PlasticContinuumParticles() {};

        virtual void initializeOtherVariables() override;
        virtual ContinuumParticles* ThisObjectPtr() override { return this; };
    };
} // namespace SPH

#endif // CONTINUUM_PARTICLES_H