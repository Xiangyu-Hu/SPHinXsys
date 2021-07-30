/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module collects all the functions required for
 * the calculation of the dynamic timeStep sizes
 */

#include"dtSizeCalc.cuh"

namespace gpu{

dtSizeCalc::dtSizeCalc(shared_ptr<ParticleData> pd,
                       shared_ptr<CellList> nl,
                       shared_ptr<System> sys,
                       Parameters par):
    pd(pd), nl(nl), sys(sys),
    stream(par.stream){
    printf("|dtSizeCalc| \tis called with dtMin: %f \n", dtMin);
}

dtSizeCalc::~dtSizeCalc(){
    printf("|dtSizeCalc| \tcall ended! \n");
}

namespace dtSizeCalc_ns{

//function to get the magnitude of acc and vel
__global__ void magAccVel(int N,
                          const int* __restrict__ groupIndex,
                          real3* __restrict__ vel,
                          real4* __restrict__ force,
                          real* __restrict__ velMag,
                          real* __restrict__ AccMag){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=N) return;
    const real3  veli = vel[groupIndex[i]];
    const real3  acci = make_real3(force[groupIndex[i]]);
    velMag[groupIndex[i]] = sqrtf(dot(veli, veli));
    AccMag[groupIndex[i]] = sqrtf(dot(acci, acci));
}

}//dtSizeCalc_ns

//function to find the max acceleration and velocity in the systems
template<int reductionOpt>
real dtSizeCalc::maxAccVel()
{
    int np = pd->getNumParticles();
    int Nthreads=128;
    int Nblocks=np/Nthreads + ((np%Nthreads)?1:0);
    auto groupIndex = nl->getGroupIndexIterator();
#ifndef _TRANSPORT_VELOCITY_
    auto   vel = pd->getVel(access::location::gpu, access::mode::read);
#else
    auto   vel = pd->getVel_tv(access::location::gpu, access::mode::read);
#endif
    auto force = pd->getForce(access::location::gpu, access::mode::read);

    velMag.resize(np);
    AccMag.resize(np);

    auto velMag_ptr = thrust::raw_pointer_cast(velMag.data());
    auto AccMag_ptr = thrust::raw_pointer_cast(AccMag.data());

    //TODO: this coud also be improved by templates when force and vel are of the same type (e.g. real3)
    dtSizeCalc_ns::magAccVel<<<Nblocks, Nthreads, 0, stream>>>(np,
                                                               groupIndex,
                                                               vel.raw(),
                                                               force.raw(),
                                                               velMag_ptr,
                                                               AccMag_ptr);

    //find the Max of reducedPar(acc or vel) among all particles
    real *maxPar;
    cudaMalloc(&maxPar, sizeof(real));
    {
        size_t newSize = 0;
        if(reductionOpt == reducedPar::Vel)
            cub::DeviceReduce::Max(nullptr, newSize, velMag_ptr, maxPar, np);
        else if(reductionOpt == reducedPar::Acc)
            cub::DeviceReduce::Max(nullptr, newSize, AccMag_ptr, maxPar, np);
        else
            throw std::runtime_error("|dtSizeCalc| \treducedPar is not valid!");


        if(newSize > tempStorage.size()){
            tempStorage.resize(newSize);
        }
    }
    size_t size = tempStorage.size();
    if(reductionOpt == reducedPar::Vel)
        cub::DeviceReduce::Max((void*)thrust::raw_pointer_cast(tempStorage.data()), size, velMag_ptr, maxPar, np);
    else if(reductionOpt == reducedPar::Acc)
        cub::DeviceReduce::Max((void*)thrust::raw_pointer_cast(tempStorage.data()), size, AccMag_ptr, maxPar, np);


    real max = 0;
    CudaSafeCall(cudaMemcpy(&max, maxPar, sizeof(real), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(maxPar));
    return max;
}

//calculate the advection time step size
//Eq. 8 in doi.org/10.1016/j.jcp.2019.109135
real dtSizeCalc::calcDtAdv(real h, real U_f){
    real cflAdv = 0.25;
    real velMax = 0.0;
    velMax = maxAccVel<reducedPar::Vel>();
    real Umax = std::max(U_f, velMax);
    //dt1 based on advection
    const real dt1 = h / (Umax + 1e-6f);
    //dt2 based on viscous terms (TODO)
    const real dt2 = std::numeric_limits<int>::max();// h*h/kinViscosity;
    //new value of the dynamic time step.
    real dtSize = cflAdv * std::min(dt1, dt2);
    if(dtSize<real(dtMin))
        dtSize=real(dtMin);
    return dtSize;
}

//calculate the acoustic time step size
//Eq. 9 in doi.org/10.1016/j.jcp.2019.109135
real dtSizeCalc::calcDtAcs(real h, real c_f){
    real cflAcs = 0.6;
    real velMax = 0.0;
    velMax = maxAccVel<reducedPar::Vel>();
    //new value of the dynamic time step.
    real dtSize = cflAcs * h / (c_f + velMax + 1e-6f);
    if(dtSize<real(dtMin))
        dtSize=real(dtMin);
    return dtSize;
}

//calculate time dynamic step size
real dtSizeCalc::calcDt(real h, real c_f){
    real CFL = 0.25;
    real accMax = 0.0;
    real velMax = 0.0;
    velMax = maxAccVel<reducedPar::Vel>();
    accMax = maxAccVel<reducedPar::Acc>();
    //dt1 based on force per unit mass.
    const real dtF = (accMax)?sqrtf(h/accMax):std::numeric_limits<int>::max();
    //dt2 based on advection
    const real dtAd = h/(std::max(c_f,velMax*10.f));
    //new value of the dynamic time step.
    real dtSize=CFL*std::min(dtF,dtAd);
    if(dtSize<real(dtMin))
        dtSize=real(dtMin);
    return dtSize;
}

}//namespace gpu
