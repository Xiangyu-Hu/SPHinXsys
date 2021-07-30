#ifndef KERNEL_CUH
#define KERNEL_CUH

#include"defines.h"

namespace gpu{
    namespace KernelFunction{
    struct M4CubicSpline{
        inline __device__ __host__ real operator()(real3 r12, real h, real boxZ){
            real r2 = dot(r12, r12);
            real r = sqrtf(r2);

            real q = abs(r)/h;

            if(q >= real(2.0)) return real(0.0);

            real twomq = real(2.0)-q;

            real W = twomq*twomq*twomq;
            if(q <= real(1.0)){
                real onemq = real(1.0) - q;
                real onemq3 = onemq*onemq*onemq;
                W -= real(4.0)*onemq3;
            }

            W *= real(1.0)/(h*h*h*real(4.0)*real(M_PI));
            return W;
        }
        inline __device__ __host__ real3 gradient(real3 r12, real h, real boxZ){

            real r2 = dot(r12, r12);
            real r = sqrtf(r2);

            real invh = real(1.0)/h;

            real q = r*invh;
            if(q >= real(2.0) ) return make_real3(0.0);

            real invh3 = invh*invh*invh;

            real3 gradW = -invh3*invh3*real(3.0)*r12/(real(4.0)*real(M_PI));


            if(q <= real(1.0)){
                gradW *= (real(-4.0)*h + real(3.0)*r);
            }
            else if(q <= real(2.0) ){
                real f = (real(2.0)*h - r);
                gradW *= f*f;

            }
            return gradW;
        }

        static inline  __device__ __host__  real getCutOff(real h){
            return real(2.0)*h;
        }

    };

    struct Wendland_C4{
        inline __device__ __host__ real operator()(real3 r12, real h, real boxZ){
            real  r2 = dot(r12, r12);
            real   r = sqrtf(r2);

            real   q = abs(r)/h;
            real   C = real();
            real val = real();


            if(boxZ == real(0.0))
                C = real(1.0)/real(M_PI)/(h*h)*real(7./4.);//for 2D
            else
                C = real(1.0)/real(M_PI)/(h*h*h)*real(21./16.);//for 3D

            // (1-q/2)^4
            real factor = (real(1.0)-real(0.5)*q)*(real(1.0)-real(0.5)*q)
                    *(real(1.0)-real(0.5)*q)*(real(1.0)-real(0.5)*q);

            if(q < real(2.0))
                val = C*factor*(real(1.0) + real(2.0)*q);
            else
                val = real(0.0);

            return val;
        }
        inline __device__ __host__ real gradient(real3 r12, real h, real boxZ){

            real  r2 = dot(r12, r12);
            real   r = sqrtf(r2);

            real   q = abs(r)/h;
            real   C = real();
            real val = real();

            if(boxZ == real(0.0))
                C = real(1.0)/real(M_PI)/(h*h)*real(7./4.);//for 2D
            else
                C = real(1.0)/real(M_PI)/(h*h*h)*real(21./16.);//for 3D

            // (1-q/2)^3
            real factor = (real(1.0)-real(0.5)*q)*(real(1.0)-real(0.5)*q)
                    *(real(1.0)-real(0.5)*q);

            if(q < real(2.0))
                val = -real(5.0)*C*factor*q/h/r;
            else
                val = real(0.0);

            return val;
        }

        static inline  __device__ __host__  real getCutOff(real h){
            return real(2.0)*h;
        }

    };
  }
}

#endif
