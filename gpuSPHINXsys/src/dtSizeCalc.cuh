/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module collects all the functions required for
 * the calculation of the dynamic timeStep sizes
 */

#ifndef DTSIZECALC_CUH
#define DTSIZECALC_CUH

#include"System.h"
#include"ParticleData.cuh"
#include"CellList.cuh"
#include<thrust/device_vector.h>
#include<memory>

namespace gpu{
class dtSizeCalc{

private:
    cudaStream_t stream;
    real dtMin = 1.e-5f;

    enum reducedPar{Acc, Vel}; //reduced Parameter: acc or vel

    thrust::device_vector<char> tempStorage; //temporary storage - cub::Reduce
    thrust::device_vector<real> velMag, AccMag;

    shared_ptr<ParticleData> pd;
    shared_ptr<CellList> nl;
    shared_ptr<System> sys;
    
public:
    struct Parameters{
        cudaStream_t stream;
    };

    dtSizeCalc(shared_ptr<ParticleData> pd,
               shared_ptr<CellList> nl,
               shared_ptr<System> sys,
               Parameters par);
    ~dtSizeCalc();

    template<int reductionOpt>
    real maxAccVel();
    real calcDt(real h, real c_f);
    real calcDtAdv(real h, real U_f);
    real calcDtAcs(real h, real c_f);

};
}//namespace gpu

#include"dtSizeCalc.cu"
#endif
  
