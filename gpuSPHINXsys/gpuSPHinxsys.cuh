/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This is the main interface to communicate between
 * Host and Device
 */

#ifndef GPUSPHINXSYS_CUH
#define GPUSPHINXSYS_CUH

#include<iostream>
#include<fstream>
#include<vector>
#include<tuple>

class gpuSPHinxsys
{
public:

    struct Parameters{
        float particle_size;
        std::tuple<float, float, float> box_size, body_force;
        std::tuple<float, float, float, std::string> probe1, probe2;
        float U_f, c_f, c_s, rho0_f, rho0_g, rho0_s;
    };

    gpuSPHinxsys(Parameters par);
    ~gpuSPHinxsys();

    void call_gpuSPHinxsys(std::vector<float> &fluidParticles,
                           std::vector<float> &wallParticles,
                           std::vector<float> &imprtdBodyParticles);

private:

    float particle_size;
    std::tuple<float, float, float> box_size, body_force;
    std::tuple<float, float, float, std::string> probe1, probe2;
    float U_f, c_f, c_s, rho0_f, rho0_g, rho0_s;

};

#endif
