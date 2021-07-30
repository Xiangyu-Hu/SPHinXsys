/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module integrates the particles dynamics using
 * different time marching algorithms
 */

#include"timeStepping.cuh"
#include"dtSizeCalc.cuh"
#include"fluidDynamics.cuh"
#include"inOut.cuh"
#include <chrono>


namespace gpu{

timeStepping::timeStepping(shared_ptr<ParticleData> pd,
                           shared_ptr<CellList> nl,
                           shared_ptr<System> sys,
                           timeStepping::Parameters par):
    pd(pd), nl(nl), sys(sys), U_f(par.U_f),
    box(par.box), h(par.h), c_f(par.c_f), probe1(par.probe1),
    probe2(par.probe2), rho0_f(par.rho0_f), bodyForce(par.bodyForce){
    printf("|timeStepping| \tis called with c_f = %.1f, rho0_f = %.1f \n", c_f, rho0_f);
    CudaSafeCall(cudaStreamCreate(&stream));
}

timeStepping::~timeStepping(){
    cudaStreamDestroy(stream);
    printf("|timeStepping| \tcall ended! \n");
}

namespace timeStepping_ns{
//Kernel for time stepping
template<int halfSteps>
__global__ void integration_ker(real4* __restrict__ pos,
                                real3* __restrict__ vel,
                                real4* __restrict__ force,
                                real3* __restrict__ F_Pb,
                                real3* __restrict__ vel_tv,
                                const int* __restrict__ groupIndex,
                                int np,
                                real dt){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=np) return;

    vel[groupIndex[i]] += make_real3(force[groupIndex[i]])*dt*real(0.5);

    //updat positions at the first half step
    if(halfSteps==timeStepping::firstHalf){
#ifndef _TRANSPORT_VELOCITY_
        real3 newPos = make_real3(pos[groupIndex[i]]) + vel[groupIndex[i]]*dt;
        pos[groupIndex[i]] = make_real4(newPos, pos[groupIndex[i]].w);
#else
        vel_tv[groupIndex[i]] = vel[groupIndex[i]] + F_Pb[groupIndex[i]]*dt*real(0.5);
        real3 newPos = make_real3(pos[groupIndex[i]]) + vel_tv[groupIndex[i]]*dt;
        pos[groupIndex[i]] = make_real4(newPos, pos[groupIndex[i]].w);
        F_Pb[groupIndex[i]] = make_real3(0);
#endif
        //Reset force
        force[groupIndex[i]] = make_real4(0);
    }
}

}//timeStepping_ns

//function for time integration
template<int halfSteps>
void timeStepping::integration(real currentDt){
    int np = pd->getNumParticles();
    int Nthreads=128;
    int Nblocks=np/Nthreads + ((np%Nthreads)?1:0);
    auto groupIndex = nl->getGroupIndexIterator();
    auto   pos = pd->getPos(access::location::gpu, access::mode::readwrite);
    auto   vel = pd->getVel(access::location::gpu, access::mode::readwrite);
    auto force = pd->getForce(access::location::gpu, access::mode::readwrite);
#ifdef _TRANSPORT_VELOCITY_
    auto   F_Pb = pd->getF_Pb(access::location::gpu, access::mode::readwrite).raw();
    auto vel_tv = pd->getVel_tv(access::location::gpu, access::mode::readwrite).raw();
#else
    auto   F_Pb = nullptr;
    auto vel_tv = nullptr;
#endif
    timeStepping_ns::integration_ker<halfSteps><<<Nblocks, Nthreads, 0, stream>>>(pos.raw(),
                                                                                  vel.raw(),
                                                                                  force.raw(),
                                                                                  F_Pb,
                                                                                  vel_tv,
                                                                                  groupIndex,
                                                                                  np,
                                                                                  currentDt);
}

//velocity verlet integration
//called from the main inerface to run the simulation
void timeStepping::velVerletIntg(){
    dtSizeCalc::Parameters parDT;
    parDT.stream = stream;
    auto dtSizeCalc_ptr = std::make_shared<dtSizeCalc>(pd, nl, sys, parDT);

    fluidDynamics::Parameters parFD;
    parFD.stream = stream;
    parFD.box = box;
    parFD.rho0_f = rho0_f;
    parFD.bodyForce = bodyForce;
    parFD.c_f = c_f;
    auto fluidDynamics_ptr = std::make_shared<fluidDynamics>(pd, nl, sys, parFD);

    auto inOut_ptr = std::make_shared<inOut>(pd, nl, box);
    //the observer's coordinate
    real3 probe1Point = make_real3(std::get<0>(probe1),
                                   std::get<1>(probe1),
                                   std::get<2>(probe1));
    real3 probe2Point = make_real3(std::get<0>(probe2),
                                   std::get<1>(probe2),
                                   std::get<2>(probe2));

    //output the initial configuration
    inOut_ptr->outputToFile<inOut::VTU, inOut::Fluid>(0);
    inOut_ptr->outputToFile<inOut::VTU, inOut::Wall>(0);
    inOut_ptr->outputToFile<inOut::VTU, inOut::thirdBody>(0);

    real rcut = Kernel::getCutOff(h);
    //TODO this updateNeighbourList can be much improved
    nl->updateNeighbourList(box, rcut, stream);
    //calculate initial number density: sigma0
    fluidDynamics_ptr->calcInitNumDensity(h);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> interval;

    //computation loop starts
    while (physical_time < End_time){

        nl->updateNeighbourList(box, rcut, stream);

        fluidDynamics_ptr->densitySumFreeSurface(h);

        integration<updateOption::firstHalf>(dt);

//        fluidDynamics_ptr->densitySumFreeSurface(h);

//        fluidDynamics_ptr->densitySumLightPhase(h);

        fluidDynamics_ptr->calcDensity<fluidDynamics::densityOption::RiemannDensity>(dt, h);
        fluidDynamics_ptr->calcPressureBC(h);
        fluidDynamics_ptr->calcForce<fluidDynamics::forceOption::RiemannForce>(h, physical_time);

        integration<updateOption::secondHalf>(dt);

        dt = dtSizeCalc_ptr->calcDt(h, c_f);

        physical_time += dt;

        if (iter_counter % screen_interval == 0){
            printf("Step: %d \tTime: %0.3f \tdt: %f \n", iter_counter, physical_time, dt);
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        //write results into a file
        if (output_counter < physical_time*output_interval){
            printf("|I/O| \tWriting output to disk ... file No. %d \n", output_counter+1);
            //write the simulation results into a file
            inOut_ptr->outputToFile<inOut::VTU, inOut::Fluid>(output_counter+1);
            inOut_ptr->outputToFile<inOut::VTU, inOut::thirdBody>(output_counter+1);
            //write the proble signals for a givnen obsever point into a file
            inOut_ptr->probeSignalToFile(h, probe1Point, std::get<3>(probe1), physical_time, output_counter+1);
            inOut_ptr->probeSignalToFile(h, probe2Point, std::get<3>(probe2), physical_time, output_counter+1);
            output_counter++;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        interval += t3 - t2;

        //resorting particles (slightly improves the performance)
        if(iter_counter%500 == 0){
            pd->sortParticles();
        }
        iter_counter++;
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> tt = t4 - t1 - interval;
    printf("Total wall clock time for computation: %.3f seconds \n", tt.count());
    printf("Total number of Iterations: %d \n", iter_counter);
}

//Dual-Criteria time integration scheme
//called from the main inerface to run the simulation
void timeStepping::dualCriteriaIntg(){
    dtSizeCalc::Parameters parDT;
    parDT.stream = stream;
    auto dtSizeCalc_ptr = std::make_shared<dtSizeCalc>(pd, nl, sys, parDT);

    fluidDynamics::Parameters parFD;
    parFD.stream = stream;
    parFD.box = box;
    parFD.rho0_f = rho0_f;
    parFD.bodyForce = bodyForce;
    parFD.c_f = c_f;
    auto fluidDynamics_ptr = std::make_shared<fluidDynamics>(pd, nl, sys, parFD);

    auto inOut_ptr = std::make_shared<inOut>(pd, nl, box);
    //output the initial configuration
    inOut_ptr->outputToFile<inOut::VTU, inOut::Fluid>(0);
    inOut_ptr->outputToFile<inOut::VTU, inOut::Wall>(0);

    real rcut = Kernel::getCutOff(h);
    //TODO this updateNeighbourList can be much improved
    nl->updateNeighbourList(box, rcut, stream);
    //calculate initial number density: sigma0
    fluidDynamics_ptr->calcInitNumDensity(h);

    //computation loop starts
    while (physical_time < End_time){

        nl->updateNeighbourList(box, rcut, stream);

        Dt = dtSizeCalc_ptr->calcDtAdv(h, U_f);
        fluidDynamics_ptr->densitySumFreeSurface(h);

        real relaxation_time = 0.0;
        while (relaxation_time < Dt){


            integration<updateOption::firstHalf>(dt);

            fluidDynamics_ptr->calcDensity<fluidDynamics::densityOption::RiemannDensity>(dt, h);
            fluidDynamics_ptr->calcPressureBC(h);
            fluidDynamics_ptr->calcForce<fluidDynamics::forceOption::RiemannForce>(h, physical_time);

            integration<updateOption::secondHalf>(dt);

            dt = dtSizeCalc_ptr->calcDtAcs(h, c_f);

            relaxation_time += dt;
            physical_time += dt;
        }

        if (iter_counter % screen_interval == 0){
            printf("Step: %d \tTime: %0.3f \tDt: %f \tdt: %f \n", iter_counter, physical_time, Dt, dt);
        }

        //write results to a file
        if (output_counter < physical_time*output_interval){
            printf("|I/O| \tWriting output to disk ... file No. %d \n", output_counter+1);
            inOut_ptr->outputToFile<inOut::VTU, inOut::Fluid>(output_counter+1);
            output_counter++;
        }

        //resorting particles
        if(iter_counter%500 == 0){
            pd->sortParticles();
        }
        iter_counter++;
    }
}

}//namspace gpu
