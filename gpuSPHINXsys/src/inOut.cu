/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module writes the particle data on an external file
 * for post-processing porpuses
 * Two output formats are available, VTU as well as PLT
 */

#include"inOut.cuh"
#include"Kernel.cuh"

using Kernel = gpu::KernelFunction::Wendland_C4;

namespace gpu{

inOut::inOut(shared_ptr<ParticleData> pd,
             shared_ptr<CellList> nl,
             Box box):
    pd(pd), nl(nl){
    printf("|inOut| \tis called  \n");
    if (fs::exists(output_folder) || fs::exists(observer_folder)){
        fs::remove_all(output_folder);
        fs::remove_all(observer_folder);
    }
    if (!fs::exists(output_folder) || !fs::exists(observer_folder)){
        fs::create_directory(output_folder);
        fs::create_directory(observer_folder);
        printf("|inOut| \toutput folders created!  \n");
    }
}

inOut::~inOut(){
    printf("|inOut| \tcall ended! \n");
}

namespace inOut_ns{

//kernel to get the physical signal at the probe
template<class NeighbourContainer>
__global__ void calcPressureSignal(NeighbourContainer ni,
                                   Box box,
                                   real h,
                                   Kernel kernel,
                                   real3 probe,
                                   real* __restrict__ probeSignal,
                                   real* __restrict__ mass,
                                   real* __restrict__ rho,
                                   real* __restrict__ p){
    real sum0 = real();
    real sum1 = real();
    auto   it = ni.begin(probe);

    while(it){
        auto neigh = *it++;
        const real4 posj = neigh.getPos();
        // include only Fluid neighboring particles
        if(posj.w == WALL) continue;
        const real3   rj = make_real3(neigh.getPos());
        const real  rhoj = rho[neigh.getGroupIndex()];
        const real    pj = p[neigh.getGroupIndex()];
        const real massj = mass[neigh.getGroupIndex()];
        const real3  rij = box.apply_pbc(probe-rj);
        sum0 += kernel(rij, h, box.boxSize.z)*massj/rhoj;
        sum1 += pj*kernel(rij, h, box.boxSize.z)*massj/rhoj;
    }
    *probeSignal = sum1/fmax(sum0, Eps);
}

}//inOut_ns


template<int outputOpt, int bodyOpt>
inline void inOut::outputToFile(int outputCount){

    const int         np = pd->getNumParticles();
    const real4 *posType = pd->getPos(access::location::cpu, access::mode::read).raw();
    const real3     *vel = pd->getVel(access::location::cpu, access::mode::read).raw();
    const real      *rho = pd->getRho(access::location::cpu, access::mode::read).raw();
    const real    *press = pd->getPressure(access::location::cpu, access::mode::read).raw();

    //seperate particles of different types
    if(bodyOpt==Fluid){
        posTypeFluid.clear();
        velFluid.clear();
        rhoFluid.clear();
        pressFluid.clear();
        for (size_t i = 0; i < np; i++){
            if (posType[i].w==LIQUID){
                posTypeFluid.push_back(posType[i]);
                velFluid.push_back(vel[i]);
                rhoFluid.push_back(rho[i]);
                pressFluid.push_back(press[i]);
            }
        }
        posTypeBody = posTypeFluid;
        velBody = velFluid;
        rhoBody = rhoFluid;
        pressBody = pressFluid;
        bodyName = "Fluid_";
        npBody = posTypeBody.size();
    }else if(bodyOpt==Wall){
        posTypeWall.clear();
        velWall.clear();
        rhoWall.clear();
        pressWall.clear();
        for (size_t i = 0; i < np; i++){
            if (posType[i].w==WALL){
                posTypeWall.push_back(posType[i]);
                velWall.push_back(vel[i]);
                rhoWall.push_back(rho[i]);
                pressWall.push_back(press[i]);
            }
        }
        posTypeBody = posTypeWall;
        velBody = velWall;
        rhoBody = rhoWall;
        pressBody = pressWall;
        bodyName = "Wall_";
        npBody = posTypeBody.size();
    }else if(bodyOpt==thirdBody){
        posTypeThirdBody.clear();
        velThirdBody.clear();
        rhoThirdBody.clear();
        pressThirdBody.clear();
        for (size_t i = 0; i < np; i++){
            if (posType[i].w==THIRDBODY){
                posTypeThirdBody.push_back(posType[i]);
                velThirdBody.push_back(vel[i]);
                rhoThirdBody.push_back(rho[i]);
                pressThirdBody.push_back(press[i]);
            }
        }
        posTypeBody = posTypeThirdBody;
        velBody = velThirdBody;
        rhoBody = rhoThirdBody;
        pressBody = pressThirdBody;
        bodyName = "ThirdBody_";
        npBody = posTypeBody.size();
    }else{
        throw std::runtime_error("|inOut| \tBody option to output is not valid!");
    }

    // write in VTU format to use in Paraview
    if(outputOpt==VTU)
    {
        std::string filefullpath = output_folder + "/" + bodyName + std::to_string(outputCount) + ".vtu";
        std::ofstream out(filefullpath.c_str(), std::ios::trunc);

        //beginning of the XML file
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << " <UnstructuredGrid>\n";
        out << "  <Piece Name =\"" << 0 << "\" NumberOfPoints=\"" << npBody << "\" NumberOfCells=\"0\">\n";

        //write position of particles
        out << "   <Points>\n";
        out << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            out << posTypeBody[i].x << " " << posTypeBody[i].y << " " << posTypeBody[i].z << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";
        out << "   </Points>\n";

        //Particles data set
        out << "   <PointData  Vectors=\"vector\">\n";
        //wrtie density
        out << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            out << rhoBody[i] << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";

        //wrtie type
        out << "    <DataArray Name=\"Type\" type=\"Int32\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            out << posTypeBody[i].w << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";

        //wrtie id
        out << "    <DataArray Name=\"Id\" type=\"Int32\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            //            out << id[i] << " ";
            out << i << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";

        //write pressure
        out << "    <DataArray Name=\"Pressure\" type=\"Float32\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            out << pressBody[i] << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";

        //write velocity
        out << "    <DataArray Name=\"Velocity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        out << "    ";
        for (size_t i = 0; i < npBody; i++) {
            out << velBody[i].x << " " << velBody[i].y << " " << velBody[i].z << " ";
        }
        out << std::endl;
        out << "    </DataArray>\n";

        //Particles data set ended
        out << "   </PointData>\n";

        //cells connectivity
        out << "   <Cells>\n";
        out << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
        out << "    </DataArray>\n";
        out << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
        out << "    </DataArray>\n";
        out << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
        out << "    </DataArray>\n";
        out << "   </Cells>\n";

        out << "  </Piece>\n";
        out << " </UnstructuredGrid>\n";
        out << "</VTKFile>\n";

        out.close();
    }

    // write in PLT format to use in TecPlot
    else if (outputOpt==PLT) {
        std::string filefullpath = output_folder + "/" + bodyName + std::to_string(outputCount) + ".plt";
        std::ofstream out(filefullpath.c_str(), std::ios::trunc);

        out<<"VARIABLES = \"x\",\"y\",\"z\",\"type\",\"vx\",\"vy\",\"vz\",\"rho\",\"p\"\n";
        for (int i = 0; i < npBody; i++){
            out << posTypeBody[i].x << " "
                << posTypeBody[i].y << " "
                << posTypeBody[i].z << " "
                << posTypeBody[i].w << " "
                << velBody[i].x << " " << velBody[i].y << " " << velBody[i].z << " "
                << rhoBody[i] << "  "
                << pressBody[i] << "\n";
        }
        out.close();
    }

    else
    {
        throw std::runtime_error("|inOut| \tOutput format is not valid!");
    }
}

//function to output the probe signal at a given probe point into a plt file
inline void inOut::probeSignalToFile(real h, real3 probe, std::string probeID, real physicalTime, int counter){

    Kernel kernel;

    real *probeSignalTmp;
    cudaMalloc(&probeSignalTmp, sizeof(real));

    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();

    real *mass = nullptr;
    if(pd->isMassAllocated())
        mass = pd->getMass(access::location::gpu, access::mode::read).raw();
    auto density = pd->getRho(access::location::gpu, access::mode::read).raw();
    auto pressure = pd->getPressure(access::location::gpu, access::mode::read).raw();

    inOut_ns::calcPressureSignal<<<1, 1>>>(ni, box, h, kernel, probe,
                                           probeSignalTmp, mass, density, pressure);

    real probeSignal = real(0.0);
    CudaSafeCall(cudaMemcpy(&probeSignal, probeSignalTmp, sizeof(real), cudaMemcpyDeviceToHost));
    CudaSafeCall(cudaFree(probeSignalTmp));

    std::string filefullpath = observer_folder + "/" + probeID + ".plt";
    std::ofstream out(filefullpath.c_str(), std::ios::app);

    if (counter == 1){
        out<<"VARIABLES=\"time\",\"probeSignal\"\n";
    }
    out <<physicalTime<<"  "<<probeSignal<<"\n";
}


}//namespace gpu
