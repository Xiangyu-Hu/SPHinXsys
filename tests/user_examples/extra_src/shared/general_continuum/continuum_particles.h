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
} // namespace SPH

#endif // CONTINUUM_PARTICLES_H