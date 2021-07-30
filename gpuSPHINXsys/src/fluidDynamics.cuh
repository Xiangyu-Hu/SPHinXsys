/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module collects all the functions required for
 * the fluid dynamics related calculations
 */

#ifndef FLUIDDYNAMICS_CUH
#define FLUIDDYNAMICS_CUH

#include"System.h"
#include"ParticleData.cuh"
#include<thrust/device_vector.h>
#include<memory>

namespace gpu{
class fluidDynamics{

private:
    cudaStream_t stream;

    shared_ptr<ParticleData> pd;
    shared_ptr<CellList> nl;
    shared_ptr<System> sys;

    real c_f, rho0_f;
    real3 bodyForce;
    Box  box;

public:
    struct Parameters{
        cudaStream_t stream;
        real c_f, rho0_f;
        real3 bodyForce;
        Box  box;
    };

    fluidDynamics(shared_ptr<ParticleData> pd,
                  shared_ptr<CellList> nl,
                  shared_ptr<System> sys,
                  Parameters par);
    ~fluidDynamics();

    enum densityOption{RiemannDensity, Continuity};
    enum forceOption{RiemannForce, Artificial};

    template<int densityOpt>
    void calcDensity(real Dt, real h);
    void calcPressureBC(real h);
    template<int forceOpt>
    void calcForce(real h, real physicalTime);
    void calcInitNumDensity(real h);
    void densitySumFreeSurface(real h);
    void densitySumLightPhase(real h);

};

}//namespace gpu

#include"fluidDynamics.cu"
#endif

