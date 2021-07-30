/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module writes the particle data on an external file
 * for post-processing porpuses
 * Two output formats are available, VTU as well as PLT
 */

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#ifndef INOUT_CUH
#define INOUT_CUH

#include"ParticleData.cuh"
#include<memory>
#include<iostream>
#include<fstream>
#include<string>
#include<experimental/filesystem>
namespace fs = std::experimental::filesystem;

namespace gpu{
class inOut{

private:

    shared_ptr<ParticleData> pd;
    shared_ptr<CellList> nl;
    Box  box;

    std::string output_folder = "./output_gpu";
    std::string observer_folder = "./observer_data";
    std::string bodyName;
    size_t npBody;

    std::vector<real4> posTypeFluid, posTypeWall, posTypeThirdBody, posTypeBody;
    std::vector<real3> velFluid, velWall, velThirdBody, velBody;
    std::vector<real> rhoFluid, rhoWall, rhoThirdBody, rhoBody, pressFluid, pressWall, pressThirdBody, pressBody;

public:

    inOut(shared_ptr<ParticleData> pd, shared_ptr<CellList> nl, Box box);
    ~inOut();

    enum outputOption{VTU, PLT};
    enum bodyOption{Fluid, Wall, thirdBody};

    template<int outputOpt, int bodyOpt>
    void outputToFile(int outputCount);

    void probeSignalToFile(real h, real3 probe, std::string probeID, real physicalTime, int counter);

};

}//namespace gpu

#include"inOut.cu"
#endif

