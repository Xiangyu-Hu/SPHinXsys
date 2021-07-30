/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This is the main interface to communicate between
 * Host and Device
 */

#include"gpuSPHinxsys.cuh"
#include"System.h"
#include"ParticleData.cuh"
#include"ParticleGroup.cuh"
#include"CellList.cuh"
#include"timeStepping.cuh"
#include"dtSizeCalc.cuh"


using gpu::access;
using real = gpu::real;
using real3 = gpu::real3;
using real4 = gpu::real4;
using timeStepping = gpu::timeStepping;


gpuSPHinxsys::gpuSPHinxsys(Parameters par):
    particle_size(par.particle_size),
    box_size(par.box_size),
    body_force(par.body_force),
    probe1(par.probe1), probe2(par.probe2),
    U_f(par.U_f), c_f(par.c_f), rho0_f(par.rho0_f),
    c_s(par.c_s), rho0_s(par.rho0_s), rho0_g(par.rho0_g)
{
    std::cout << "\t*************************************************" << "\n";
    std::cout << "\t       gpuSPHinXsys: An SPH solver on GPUs       " << "\n";
    std::cout << "\t*************************************************" << "\n";
    //    /*checkCudaErrors*/(cudaMalloc((void **)&devicePar->VelMax, sizeof(real)));
}

gpuSPHinxsys::~gpuSPHinxsys()
{
    std::cout << "\t*************************************************" << "\n";
    std::cout << "\t             gpuSPHinXsys call ended!            " << "\n";
    std::cout << "\t*************************************************" << "\n";
    //    if(devicePar->VelMax != NULL) cudaFree(devicePar->VelMax);
}

struct GPU{
    std::shared_ptr<gpu::ParticleData> pd;
    std::shared_ptr<gpu::System> sys;
    std::shared_ptr<gpu::CellList> nl; //neighborlist
};

//initialize a GPU system
GPU initializeGPU(size_t np){
    auto sys = std::make_shared<gpu::System>();
    auto  pd = std::make_shared<gpu::ParticleData>(np, sys);
    auto  pg = std::make_shared<gpu::ParticleGroup>(pd, sys, "All");
    auto  nl = std::make_shared<gpu::CellList>(pd, pg, sys);
    return {pd, sys, nl};
}

//Copy particle position_type to device and allocate density and mass
void copyPositionsToDevice(const GPU & gpu_state,
                           const std::vector<float> &fluidParticles,
                           const std::vector<float> &wallParticles,
                           const std::vector<float> &thirdBodyParticles,
                           float rho0_fluid,
                           float rho0_gas,
                           int DIM,
                           float dp){
    auto  pos_type = gpu_state.pd->getPos(access::location::cpu, access::mode::write);
    auto       rho = gpu_state.pd->getRho(access::location::cpu, access::mode::write);
    auto      rho0 = gpu_state.pd->getRho0(access::location::cpu, access::mode::write);
    auto      mass = gpu_state.pd->getMass(access::location::cpu, access::mode::write);
    auto       vol = gpu_state.pd->getVol(access::location::cpu, access::mode::write);
    size_t npFluid = fluidParticles.size()/3;
    size_t  npWall = wallParticles.size()/3;
    size_t  npThirdBody = thirdBodyParticles.size()/3;

    for(int i = 0; i<npFluid; i++){
        pos_type[i].x = fluidParticles[3*i];
        pos_type[i].y = fluidParticles[3*i+1];
        pos_type[i].z = fluidParticles[3*i+2];
        pos_type[i].w = LIQUID; //fluid particles have type 0
        rho[i] = rho0_fluid;
        rho0[i] = rho0_fluid;
        mass[i] = std::pow(dp, DIM)*rho[i];
        vol[i] = mass[i]/rho[i];
        // check NaN in fluid pos
        if (std::isnan(pos_type[i].x) or
            std::isnan(pos_type[i].y) or
            std::isnan(pos_type[i].z)){
            std::cout << "ERROR in Fluid position copied to Device!" <<"\n";
            exit(1);
        }
    }
    for(int i = 0; i<npWall; i++){
        pos_type[i+npFluid].x = wallParticles[3*i];
        pos_type[i+npFluid].y = wallParticles[3*i+1];
        pos_type[i+npFluid].z = wallParticles[3*i+2];
        pos_type[i+npFluid].w = WALL; //wall particles have type 1
        rho[i+npFluid] = rho0_fluid;
        rho0[i+npFluid] = rho0_fluid;
        mass[i+npFluid] = std::pow(dp, DIM)*rho[i+npFluid];
        vol[i+npFluid] = mass[i+npFluid]/rho[i+npFluid];
        // check NaN in wall pos
        if (std::isnan(pos_type[i+npFluid].x) or
            std::isnan(pos_type[i+npFluid].y) or
            std::isnan(pos_type[i+npFluid].z)){
            std::cout << "ERROR in Wall position copied to Device!" <<"\n";
            exit(1);
        }
    }
    for(int i = 0; i<npThirdBody; i++){
        pos_type[i+npFluid+npWall].x = thirdBodyParticles[3*i];
        pos_type[i+npFluid+npWall].y = thirdBodyParticles[3*i+1];
        pos_type[i+npFluid+npWall].z = thirdBodyParticles[3*i+2];
        pos_type[i+npFluid+npWall].w = THIRDBODY; //3rd body particles have type 2
        rho[i+npFluid+npWall] = rho0_gas;
        rho0[i+npFluid+npWall] = rho0_gas;
        mass[i+npFluid+npWall] = std::pow(dp, DIM)*rho[i+npFluid+npWall];
        vol[i+npFluid+npWall] = mass[i+npFluid+npWall]/rho[i+npFluid+npWall];
        // check NaN in THIRDBODY pos
        if (std::isnan(pos_type[i+npFluid+npWall].x) or
            std::isnan(pos_type[i+npFluid+npWall].y) or
            std::isnan(pos_type[i+npFluid+npWall].z)){
            std::cout << "ERROR in Thirdbody position copied to Device!" <<"\n";
            exit(1);
        }
    }
}


// set initial velocity and pressure to zero
void initialize_properties(const GPU & gpu_state){
    auto vel = gpu_state.pd->getVel(access::location::cpu, access::mode::write);
    auto   p = gpu_state.pd->getPressure(access::location::cpu, access::mode::write);
    std::fill(vel.begin(), vel.end(), real3());
    std::fill(p.begin(), p.end(), real());
#ifdef _TRANSPORT_VELOCITY_
    auto vel_tv = gpu_state.pd->getVel_tv(access::location::cpu, access::mode::write);
    auto   F_Pb = gpu_state.pd->getF_Pb(access::location::cpu, access::mode::write);
    std::fill(vel_tv.begin(), vel_tv.end(), real3());
    std::fill(F_Pb.begin(), F_Pb.end(), real3());
#endif
}

void gpuSPHinxsys::call_gpuSPHinxsys(std::vector<float> &fluidParticles,
                                     std::vector<float> &wallParticles,
                                     std::vector<float> &thirdBodyParticles){

    int     npFluid = fluidParticles.size()/3;
    int      npWall = wallParticles.size()/3;
    int npThirdBody = thirdBodyParticles.size()/3;
    int          np = npFluid + npWall + npThirdBody;

    printf("Copied to device: npFluid: %d npWall: %d npThirdBody: %d np: %d \n", npFluid, npWall, npThirdBody, np);

    //copy the box size
    real3 boxSize = gpu::make_real3(std::get<0>(box_size),
                                    std::get<1>(box_size),
                                    std::get<2>(box_size));
    //copy the external body force
    real3 bodyForceTmp = gpu::make_real3(std::get<0>(body_force),
                                          std::get<1>(body_force),
                                          std::get<2>(body_force));

    int DIM = boxSize.z > real(0.0)?3:2;

    auto gpu_state = initializeGPU(np);

    copyPositionsToDevice(gpu_state, fluidParticles, wallParticles,
                          thirdBodyParticles, rho0_f, rho0_g, DIM, particle_size);

    //creat a pointer to the timeStepping module
    timeStepping::Parameters parTS;
    parTS.U_f = U_f;
    parTS.box = gpu::Box(boxSize);
    parTS.h = 1.3*particle_size;
    parTS.c_f = c_f;
    parTS.rho0_f = rho0_f;
    parTS.bodyForce = bodyForceTmp;
    parTS.probe1 = probe1;
    parTS.probe2 = probe2;
    auto timeStepping_ptr = std::make_shared<timeStepping>(gpu_state.pd, gpu_state.nl, gpu_state.sys, parTS);

    initialize_properties(gpu_state);

    gpu_state.pd->sortParticles();

    //the time integration functions to Run the simulation
    timeStepping_ptr->velVerletIntg(); //velocity Verlet scheme
//    timeStepping_ptr->dualCriteriaIntg(); //dual criteria scheme
//    timeStepping_ptr->solidDynaIntg(); //integration scheme for solid dynamics

    gpu_state.sys->finish();
}
