/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module integrates the particles dynamics using
 * different time marching algorithms
 */

#ifndef TIMESTEPPING_CUH
#define TIMESTEPPING_CUH

#include"System.h"
#include"ParticleData.cuh"
#include"CellList.cuh"
#include<memory>

namespace gpu{
class timeStepping{

private:
    shared_ptr<ParticleData> pd;
    shared_ptr<CellList> nl;
    shared_ptr<System> sys;

    real U_f, h, c_f, rho0_f;
    real3 bodyForce;
    std::tuple<float, float, float, std::string> probe1, probe2;
    Box  box;
    
    cudaStream_t stream;

    real dt = 1.e-5f;           //acoustic time step size
    real Dt = 1.e-5f;           //advection time step size
    real End_time = 10.0;       //End point of the simulation
    int  screen_interval = 100; //screen info output interval
    int  output_interval = 5;  //number of outputs per second
    int  output_counter = 0;    //output counter
    int  iter_counter = 0;      //number of iterations
    real physical_time = 0.0;   //physical simulation time

    template<int halfSteps>
    void integration(real currentDt);

public:
    struct Parameters{
        real U_f, h, c_f, rho0_f;
        real3 bodyForce;
        std::tuple<float, float, float, std::string> probe1, probe2;
        Box  box;
    };

    enum updateOption{firstHalf, secondHalf};

    timeStepping(shared_ptr<ParticleData> pd,
                 shared_ptr<CellList> nl,
                 shared_ptr<System> sys,
                 Parameters par);

    ~timeStepping();

    void velVerletIntg();
    void dualCriteriaIntg();

    /// TODO: this funcs will be needed when no pd, via cudaMemcpy
    /*template<class InOut>void outputToFile(InOut InOut_ptr, int outputOpt, int outputCount);*/

};

}//namespace gpu

#include"timeStepping.cu"
#endif
  
