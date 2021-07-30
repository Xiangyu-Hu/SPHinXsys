/*
 * Massoud Rezavand 2019.
 * Technical University of Munich
 * gpuSPHinxsys - An SPH solver for CUDA enabled GPUs
 *
 * This module collects all the functions required for
 * the fluid dynamics related calculations
 */

#include"fluidDynamics.cuh"
#include"Kernel.cuh"

using Kernel = gpu::KernelFunction::Wendland_C4;

namespace gpu{

fluidDynamics::fluidDynamics(shared_ptr<ParticleData> pd,
                             shared_ptr<CellList> nl,
                             shared_ptr<System> sys,
                             Parameters par):
    pd(pd), nl(nl), sys(sys),
    c_f(par.c_f), rho0_f(par.rho0_f), bodyForce(par.bodyForce),
    stream(par.stream), box(par.box){
    printf("|fluidDynamics| \tis called with c_f = %.1f, rho0_f = %.1f \n", c_f, rho0_f);
}

fluidDynamics::~fluidDynamics(){
    printf("|fluidDynamics| \tcall ended! \n");
}

namespace fluidDynamics_ns{

//Kernel to calculate drho/dt by Riemann Solvers (Continuity Eq.)
template<class NeighbourContainer, class Kernel>
__global__ void calcDensityRiemann_ker(NeighbourContainer ni,
                                       Kernel kernel,
                                       const real4* __restrict__ sortPos,
                                       const int* __restrict__ groupIndex,
                                       int np, Box box,
                                       real3* __restrict__ vel,
                                       real* __restrict__ rho,
                                       real* __restrict__ mass,
                                       real* __restrict__ p, real h, real dt,
                                       real c_f,
                                       const real* __restrict__ rho0,
                                       real* __restrict__ vol){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on LIQUID particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w != WALL){
        //Set ni to provide iterators for particle i
        ni.set(i);

        const real3   ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        const real3 veli = vel[groupIndex[i]];
        const real  rhoi = rho[groupIndex[i]];
        const real    pi = p[groupIndex[i]];
        real        drho = real();

        auto it = ni.begin(); //Iterator to the first neighbour of particle i

        while(it){
            auto neigh = *it++;
            if(neigh.getGroupIndex() == groupIndex[i]) continue; //skip if same particle

            const real3   rj = make_real3(neigh.getPos());
            const int  typej = neigh.getPos().w;
            const real3 velj = vel[neigh.getGroupIndex()];
            const real  rhoj = rho[neigh.getGroupIndex()];
            const real massj = mass[neigh.getGroupIndex()];
            const real  volj = vol[neigh.getGroupIndex()];
            const real3  rij = box.apply_pbc(ri-rj);

            //low dissipation Riemann problem
            real r2 = dot(rij, rij);
            real dist = sqrtf(r2);
            real3 _rij = rj - ri;
            real pj = p[neigh.getGroupIndex()];
            real3 e_ij = _rij*1/(dist + 1.0e-15);
            real ul = dot(e_ij, veli);
            real ur = dot(e_ij, velj);
            real v_star = (rhoi*ul+rhoj*ur+(pi-pj)/c_f)/(rhoi+rhoj);
            real aw = kernel.gradient(rij, h, box.boxSize.z);
            //only volume of wall particles into account
            if (typej == WALL)
                drho += 2.0*rhoi*volj*(v_star-ul)*aw*dist;
            else
                drho += 2.0*rhoi*massj/rhoj*(v_star-ul)*aw*dist;
        }
        rho[groupIndex[i]] += drho*dt;
        //get the volume according to rho
        vol[groupIndex[i]] = mass[groupIndex[i]]/rho[groupIndex[i]];
        // pressure calculation via the linear EoS
        p[groupIndex[i]] = c_f*c_f*(rho[groupIndex[i]] - rho0[groupIndex[i]]);
//        printf("rho = %f  and rho0 = %f \n",rho[groupIndex[i]], rho0[groupIndex[i]] );

    }
}

//Kernel to calculate drho/dt using Artificial viscosity (Continuity Eq.)
template<class NeighbourContainer, class Kernel>
__global__ void calcDensityArtificial_ker(NeighbourContainer ni,
                                          Kernel kernel,
                                          const real4* __restrict__ sortPos,
                                          const int* __restrict__ groupIndex,
                                          int np, Box box,
                                          real3* __restrict__ vel,
                                          real* __restrict__ rho,
                                          real* __restrict__ mass,
                                          real* __restrict__ p,
                                          real h, real dt,
                                          real c_f,
                                          const real* __restrict__ rho0){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on LIQUID and gas particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w != WALL){
        //Set ni to provide iterators for particle i
        ni.set(i);

        const real3   ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        const real3 veli = vel[groupIndex[i]];
        const real  rhoi = rho[groupIndex[i]];
        real        drho = real();

        auto it = ni.begin(); //Iterator to the first neighbour of particle i

        while(it){
            auto neigh = *it++;
            if(neigh.getGroupIndex() == groupIndex[i]) continue; //skip if same particle

            const real3    rj = make_real3(neigh.getPos());
            const real3  velj = vel[neigh.getGroupIndex()];
            const real   rhoj = rho[neigh.getGroupIndex()];
            const real  massj = mass[neigh.getGroupIndex()];

            const real3   rij = box.apply_pbc(ri-rj);
            //Artificial viscosity
            const real3 velij = veli - velj;
            //TODO
            //in this AV implemetaion we should use the Cubic kernel which indludes rij in there
            const real3 kernel_grad = /*rij**/rij*kernel.gradient(rij, h, box.boxSize.z);
            drho += rhoi*massj/rhoj*dot(kernel_grad, velij);
        }
        rho[groupIndex[i]] += drho*dt;
        // pressure calculation via the linear EoS
        p[groupIndex[i]] = c_f*c_f*(rho[groupIndex[i]] - rho0[groupIndex[i]]);
    }
}

//Kernel to calculate initail number density sigma0
template<class NeighbourContainer, class Kernel>
__global__ void calcInitNumDensity_ker(NeighbourContainer ni,
                                       Kernel kernel,
                                       const real4* __restrict__ sortPos,
                                       const int* __restrict__ groupIndex,
                                       real* __restrict__ sigma0,
                                       int np, Box box,
                                       real h, real rho0_f){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    ni.set(i);
    const real3 ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
    real  sum0 = real();
    auto it = ni.begin();
    while(it){
        auto neigh = *it++;
        const real3  rj = make_real3(neigh.getPos());
        const real3 rij = box.apply_pbc(ri-rj);
        sum0 += kernel(rij, h, box.boxSize.z);
    }
    sigma0[groupIndex[i]] = sum0;
}

//Kernel to update density using summation for free surface cases
template<class NeighbourContainer, class Kernel>
__global__ void densitySumFreeSurface_ker(NeighbourContainer ni,
                                          Kernel kernel,
                                          const real4* __restrict__ sortPos,
                                          const int* __restrict__ groupIndex,
                                          real* __restrict__ sigma0,
                                          real* __restrict__ rho,
                                          int np, Box box,
                                          real h,
                                          const real* __restrict__ rho0){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w != WALL){
        ni.set(i);
        const real3 ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        const real sigma0i = sigma0[groupIndex[i]];
        real rhoi = rho[groupIndex[i]];
        real  sigma = real();
        auto it = ni.begin();
        while(it){
            auto neigh = *it++;
            const real4 posj = neigh.getPos();
            // include only Fluid neighboring particles
            if((posi.w ==LIQUID and posj.w == THIRDBODY) /*or
                    (posi.w == THIRDBODY and posj.w ==LIQUID) or
                        (posi.w == THIRDBODY and posj.w ==THIRDBODY)*/) continue;
            const real3  rj = make_real3(neigh.getPos());
            const real3 rij = box.apply_pbc(ri-rj);
            sigma += kernel(rij, h, box.boxSize.z);
        }
        real rhoSum = sigma * rho0[groupIndex[i]] / sigma0i;
        rho[groupIndex[i]] = rhoSum + fmax(0.0f, (rhoi - rhoSum)) * rho0[groupIndex[i]] / rhoi;
    }
}

//Kernel to update density using summation for the lighter phase
template<class NeighbourContainer, class Kernel>
__global__ void densitySumLightPhase_ker(NeighbourContainer ni,
                                         Kernel kernel,
                                         const real4* __restrict__ sortPos,
                                         const int* __restrict__ groupIndex,
                                         real* __restrict__ sigma0,
                                         real* __restrict__ rho,
                                         int np, Box box, real h,
                                         const real* __restrict__ rho0,
                                         real* __restrict__ p,
                                         const real c_f,
                                         real* __restrict__ vol,
                                         const real* __restrict__ mass){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on gas particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w == THIRDBODY){
        ni.set(i);
        const real3 ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));

        const real   rho0i = rho0[groupIndex[i]];

        real  sum = real();
        auto it = ni.begin();
        while(it){
            auto neigh = *it++;
            const real3   rj = make_real3(neigh.getPos());
            const real3  rij = box.apply_pbc(ri-rj);

            sum += kernel(rij, h, box.boxSize.z);
            /*or in a total Lagrangian form: */
            //const real rho0j = rho0[neigh.getGroupIndex()];
            //sum += kernel(rij, h, box.boxSize.z)*2.f*rho0i/(rho0i+rho0j);
        }
        rho[groupIndex[i]] = sum * mass[groupIndex[i]];
        /*or in a total Lagrangian form: */
        //const real sigma0i = sigma0[groupIndex[i]];
        //rho[groupIndex[i]] = sum * rho0i / sigma0i;

        //get the volume according to rho
        vol[groupIndex[i]] = mass[groupIndex[i]]/rho[groupIndex[i]];
        // pressure calculation via the linear EoS
        p[groupIndex[i]] = c_f*c_f*(rho[groupIndex[i]] - rho0i);
    }

}


//Kernel to calculate pressure for wall particles
template<class NeighbourContainer, class Kernel>
__global__ void calcPressureBC_ker(NeighbourContainer ni,
                                   Kernel kernel,
                                   const real4* __restrict__ sortPos,
                                   const int* __restrict__ groupIndex,
                                   int np, Box box,
                                   real* __restrict__ rho,
                                   real* __restrict__ p,
                                   real h, real c_f,
                                   const real* __restrict__ rho0,
                                   real3 bodyForce){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on Wall particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w == WALL){
        //Set ni to provide iterators for particle i
        ni.set(i);

        const real3 ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        real  sum0 = real();
        real  sum1 = real();
        real3 sum2 = real3();
        real3   aw = real3();
        //gravity
        //in case coordinates are needed ow. directly can be used
        const real3 body_force = bodyForce;

        auto it = ni.begin(); //Iterator to the first neighbour of particle i

        while(it){
            auto neigh = *it++;
            if(neigh.getGroupIndex() == groupIndex[i]) continue; //skip if same particle
            const real4 posj = neigh.getPos();
            // include only Fluid neighboring particles
            if(posj.w == WALL) continue;

            const real3  rj = make_real3(neigh.getPos());
            const real rhoj = rho[neigh.getGroupIndex()];
            const real   pj = p[neigh.getGroupIndex()];
            const real3 rij = box.apply_pbc(ri-rj);
            const real  wij = kernel(rij, h, box.boxSize.z);
//            printf("rhoj gas = %f  \trho0j gas = %f \n", rho[neigh.getGroupIndex()], rho0[neigh.getGroupIndex()]);

            //fraction devided by rho_j to get pressure from the lighter fluid
            sum0 += wij/rhoj;
            sum1 += pj*wij/rhoj;
            sum2 += rij*wij*rhoj/rhoj;

        }
        aw = body_force;// - aw; //for later developments
        real tmp = real();
        tmp = dot(aw, sum2);
        p[groupIndex[i]] = (sum1+tmp)/(sum0+1.e-20);
        //get density for wall particles
        rho[groupIndex[i]] = p[groupIndex[i]]/(c_f*c_f) + rho0[groupIndex[i]];
//        printf("rho gas = %f  \trho0 gas = %f \n", rho[groupIndex[i]], rho0[groupIndex[i]]);
    }
}

//Kernel to calculate pressure and viscosity related forces via RiemannSolvers
template<class NeighbourContainer, class Kernel>
__global__ void calcForceRiemann_ker(NeighbourContainer ni,
                                     Kernel kernel,
                                     const real4* __restrict__ sortPos,
                                     const int* __restrict__ groupIndex,
                                     int np, Box box,
                                     real4* __restrict__ force,
                                     real3* __restrict__ vel,
                                     real* __restrict__ rho,
                                     real* __restrict__ mass,
                                     real* __restrict__ p,
                                     real3* __restrict__ vel_tv,
                                     real3* __restrict__ F_Pb, real P_b,
                                     real h, real c_f, real3 bodyForce,
                                     real physicalTime,
                                     real* __restrict__ vol){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on Fluid particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w != WALL){
        //gravity
#ifdef _TIMEDEPENDENT_BODYFORCE_    //for sloshing tank
         const real4 body_force = make_real4(bodyForce.x*sin(2.*M_PI*0.496*physicalTime)
                                             ,bodyForce.y, 0., 0.);
#else
         const real4 body_force = make_real4(bodyForce, 0.);
#endif

        //Set ni to provide iterators for particle i
        ni.set(i);

        const real3   ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        const real  rhoi = rho[groupIndex[i]];
        const real    pi = p[groupIndex[i]];
        const real3 veli = vel[groupIndex[i]];
        real3 F1 = real3();

        auto it = ni.begin();

        while(it){
            auto neigh = *it++;

            if(neigh.getGroupIndex() == groupIndex[i]) continue; //skip if same particle

            const real3   rj = make_real3(neigh.getPos());
            const real  rhoj = rho[neigh.getGroupIndex()];
            const real    pj = p[neigh.getGroupIndex()];
            const real massj = mass[neigh.getGroupIndex()];
            const int  typej = neigh.getPos().w;
            real volj;
            //only volume of wall particles into account
            if(typej != WALL)
                volj = massj/rhoj;
            else
                volj = vol[neigh.getGroupIndex()];
            const real3 velj = vel[neigh.getGroupIndex()];
            const real3  rij = box.apply_pbc(ri-rj);
            //low dissipation Riemann problem
            const real3 kernel_grad = rij*kernel.gradient(rij, h, box.boxSize.z);
            const real4 posj = neigh.getPos();
            if(posj.w != WALL){
                real     r2 = dot(rij, rij);
                real   dist = sqrtf(r2);
                real3  _rij = rj - ri;
                real3  e_ij = _rij*1/(dist + 1.0e-15);
                real     ul = dot(e_ij, veli);
                real     ur = dot(e_ij, velj);
                real p_star = (rhoi*pj+rhoj*pi+rhoi*rhoj*c_f*(ul-ur)*
                               fmin(real(3.0)*fmax((ul-ur)/c_f, real(0.0)), real(1.0)))/(rhoi+rhoj);
                real  temp1 = -2.0*p_star*volj/rhoi;
                F1 += temp1*kernel_grad;
            }else{
                //exclude the second term only when Fluid-Wall interaction
                real p_star = (rhoj*pi + rhoi*pj)/(rhoi + rhoj);
                real temp1 = -2.0*p_star*volj/rhoi;
                F1 += temp1*kernel_grad;
            }
            //transport velocity formulation
#ifdef _TRANSPORT_VELOCITY_
            if(posi.w == THIRDBODY){
                const real3 v_tv_i = vel_tv[groupIndex[i]];
                const real3 v_tv_j = vel_tv[neigh.getGroupIndex()];
                const real   massi = mass[groupIndex[i]];
                const real    voli = massi/rhoi;
                const real    coef = real(1.)/massi*(voli*voli + volj*volj);
                // artificial stress tensor Aij (Adami et al. 2013)
                real3 A_ij = real3();
                // x component
                real3 Ax_i = (v_tv_i - veli) * rhoi * veli.x;
                real3 Ax_j = (v_tv_j - velj) * rhoj * velj.x;
                A_ij.x = real(0.5)*dot((Ax_i+Ax_j), kernel_grad);
                // y component
                real3 Ay_i = (v_tv_i - veli) * rhoi * veli.y;
                real3 Ay_j = (v_tv_j - velj) * rhoj * velj.y;
                A_ij.y = real(0.5)*dot((Ay_i+Ay_j), kernel_grad);
                // z component
                real3 Az_i = (v_tv_i - veli) * rhoi * veli.z;
                real3 Az_j = (v_tv_j - velj) * rhoj * velj.z;
                A_ij.z = real(0.5)*dot((Az_i+Az_j), kernel_grad);

                real3 dF_AS = A_ij * coef;
                F1 += dF_AS;
                // background pressure force
                real P_b1 = 5.*0.001*c_f*c_f;
                real     temp3 = real(-2.)*P_b1*volj/rhoi;
                real3  dF_Pb = kernel_grad * temp3;
                F_Pb[groupIndex[i]] += dF_Pb;
            }
#endif
        }
        force[groupIndex[i]] = make_real4(F1, 0);
//        if(posi.w == LIQUID)
        force[groupIndex[i]] += body_force;
    }
}

//Kernel to calculate pressure and Artificial viscosity forces
template<class NeighbourContainer, class Kernel>
__global__ void calcForceArtificial_ker(NeighbourContainer ni,
                                        Kernel kernel,
                                        const real4* __restrict__ sortPos,
                                        const int* __restrict__ groupIndex,
                                        int np, Box box,
                                        real4* __restrict__ force,
                                        real3* __restrict__ vel,
                                        real* __restrict__ rho,
                                        real* __restrict__ mass,
                                        real* __restrict__ p,
                                        real h, real c_f, real3 bodyForce){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i >= np) return;
    //only on Fluid particles
    const real4 posi = cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i);
    if(posi.w == LIQUID){
        //gravity
        const real4 body_force = make_real4(bodyForce, 0.);

        //Set ni to provide iterators for particle i
        ni.set(i);

        const real3   ri = make_real3(cub::ThreadLoad<cub::LOAD_LDG>(sortPos + i));
        const real  rhoi = rho[groupIndex[i]];
        const real    pi = p[groupIndex[i]];
        const real massi = mass[groupIndex[i]];
        const real  voli = massi/rhoi;
        const real3 veli = vel[groupIndex[i]];

        real3 F1 = real3();
        real3 F2 = real3();

        auto it = ni.begin();

        while(it){
            auto neigh = *it++;

            if(neigh.getGroupIndex() == groupIndex[i]) continue; //skip if same particle

            const real3   rj = make_real3(neigh.getPos());
            const real  rhoj = rho[neigh.getGroupIndex()];
            const real    pj = p[neigh.getGroupIndex()];
            const real massj = mass[neigh.getGroupIndex()];
            const real  volj = massj/rhoj;
            const real3 velj = vel[neigh.getGroupIndex()];
            const real3  rij = box.apply_pbc(ri-rj);
            //Artificial viscosity
            //TODO
            //in this AV implemetaion we should use the Cubic kernel which indludes rij in there
            const real3 kernel_grad = /*rij**/rij*kernel.gradient(rij, h, box.boxSize.z);
            const real temp0 = 1./massi*(voli*voli + volj*volj);
            const real   pij = (rhoj*pi + rhoi*pj)/(rhoi + rhoj);
            const real temp1 = -1.*pij*temp0;
            F1 += temp1*kernel_grad;

            const real4   posj = neigh.getPos();
            //TODO
            //if(posj.w == LIQUID){ //free-slip
            const real   alpha = 0.1;
            const real epsilon = 0.001;
            const real   rhoij = (rhoi + rhoj)/2.;
            const real3  velij = veli - velj;
            const real  vij_dr = dot(velij, rij)/(dot(rij, rij)+epsilon*h*h);
            const real    visc = -massj*alpha*c_f*h*vij_dr/rhoij;
            F2 += visc*kernel_grad;
            //}
        }
        force[groupIndex[i]] = make_real4(F1+F2, 0);
        force[groupIndex[i]] += body_force;
    }
}

}//namspace fluidDynamics_ns


//calculate density
template<int densityOpt>
void fluidDynamics::calcDensity(real Dt, real h){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto vel = pd->getVel(access::location::gpu, access::mode::readwrite).raw();
    //If mass is not allocated assume all masses are 1
    real *mass = nullptr;
    if(pd->isMassAllocated())
        mass = pd->getMass(access::location::gpu, access::mode::read).raw();
    auto rho = pd->getRho(access::location::gpu, access::mode::readwrite).raw();
    auto rho0 = pd->getRho0(access::location::gpu, access::mode::readwrite).raw();
    auto pressure = pd->getPressure(access::location::gpu, access::mode::readwrite).raw();
    auto vol = pd->getVol(access::location::gpu, access::mode::readwrite).raw();

    if(densityOpt==densityOption::RiemannDensity){
        fluidDynamics_ns::calcDensityRiemann_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel,
                                                                               sortPos, groupIndex,
                                                                               np, box,
                                                                               vel, rho,
                                                                               mass, pressure,
                                                                               h, Dt, c_f, rho0, vol);
    }else if(densityOpt==densityOption::Continuity){
        fluidDynamics_ns::calcDensityArtificial_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel,
                                                                                  sortPos, groupIndex,
                                                                                  np, box,
                                                                                  vel, rho,
                                                                                  mass, pressure,
                                                                                  h, Dt, c_f, rho0);
    }else {
        throw std::runtime_error("|fluidDynamics| \tdensityOpt is not valid!");
    }
}

//calculate initial number density: sigma0
void fluidDynamics::calcInitNumDensity(real h){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto sigma0 = pd->getSigma0(access::location::gpu, access::mode::readwrite).raw();

    fluidDynamics_ns::calcInitNumDensity_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                           groupIndex, sigma0,
                                                                           np, box, h, rho0_f);
}

//density calculation using summation for free surface cases
void fluidDynamics::densitySumFreeSurface(real h){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto sigma0 = pd->getSigma0(access::location::gpu, access::mode::readwrite).raw();
    auto density = pd->getRho(access::location::gpu, access::mode::readwrite).raw();
    auto rho0 = pd->getRho0(access::location::gpu, access::mode::readwrite).raw();

    fluidDynamics_ns::densitySumFreeSurface_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                              groupIndex, sigma0,
                                                                              density,
                                                                              np, box,
                                                                              h, rho0);
}

//density calculation using summation for gas phase
void fluidDynamics::densitySumLightPhase(real h){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto sigma0 = pd->getSigma0(access::location::gpu, access::mode::readwrite).raw();
    auto density = pd->getRho(access::location::gpu, access::mode::readwrite).raw();
    auto rho0 = pd->getRho0(access::location::gpu, access::mode::readwrite).raw();
    auto pressure = pd->getPressure(access::location::gpu, access::mode::readwrite).raw();
    auto vol = pd->getVol(access::location::gpu, access::mode::readwrite).raw();
    auto mass = pd->getMass(access::location::gpu, access::mode::readwrite).raw();

    fluidDynamics_ns::densitySumLightPhase_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                             groupIndex, sigma0,
                                                                             density,
                                                                             np, box,
                                                                             h, rho0,
                                                                             pressure, c_f,
                                                                             vol, mass);
}

//calculate pressure for wall particles
void fluidDynamics::calcPressureBC(real h){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto d_density = pd->getRho(access::location::gpu, access::mode::readwrite).raw();
    auto pressure = pd->getPressure(access::location::gpu, access::mode::readwrite).raw();
    auto rho0 = pd->getRho0(access::location::gpu, access::mode::readwrite).raw();

    fluidDynamics_ns::calcPressureBC_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                       groupIndex, np, box,
                                                                       d_density, pressure,
                                                                       h, c_f, rho0, bodyForce);
}

//calculate forces
template<int forceOpt>
void fluidDynamics::calcForce(real h, real physicalTime){
    int np = pd->getNumParticles();
    Kernel kernel;
    // get a NeighborContainer
    auto ni = nl->getNeighbourContainer();
    auto sortPos = nl->getPositionIterator();
    auto groupIndex = nl->getGroupIndexIterator();
    auto vel = pd->getVel(access::location::gpu, access::mode::readwrite).raw();
    //If mass is not allocated assume all masses are 1
    real *d_mass = nullptr;
    if(pd->isMassAllocated())
        d_mass = pd->getMass(access::location::gpu, access::mode::read).raw();
    auto rho = pd->getRho(access::location::gpu, access::mode::readwrite).raw();
    auto pressure = pd->getPressure(access::location::gpu, access::mode::readwrite).raw();
    auto force = pd->getForce(access::location::gpu, access::mode::readwrite).raw();
    auto vol = pd->getVol(access::location::gpu, access::mode::readwrite).raw();
#ifdef _TRANSPORT_VELOCITY_
    real    P_b = 5.0f*1.f*c_f*c_f; //the background pressure
    auto vel_tv = pd->getVel_tv(access::location::gpu, access::mode::readwrite).raw();
    auto   F_Pb = pd->getF_Pb(access::location::gpu, access::mode::readwrite).raw();
#else
    real    P_b = 0.0f;
    auto vel_tv = nullptr;
    auto   F_Pb = nullptr;
#endif

    if(forceOpt==forceOption::RiemannForce){
        fluidDynamics_ns::calcForceRiemann_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                             groupIndex, np, box,
                                                                             force, vel,
                                                                             rho, d_mass,
                                                                             pressure,
                                                                             vel_tv,
                                                                             F_Pb, P_b,
                                                                             h, c_f, bodyForce,
                                                                             physicalTime,
                                                                             vol);
    }else if(forceOpt==forceOption::Artificial){
        fluidDynamics_ns::calcForceArtificial_ker<<<np/128+1, 128, 0, stream>>>(ni, kernel, sortPos,
                                                                                groupIndex, np, box,
                                                                                force, vel,
                                                                                rho, d_mass,
                                                                                pressure,
                                                                                h, c_f, bodyForce);
    }else {
        throw std::runtime_error("|fluidDynamics| \tforceOpt is not valid!");
    }
}


}//namespace gpu
