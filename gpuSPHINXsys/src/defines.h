#ifndef DEFINES_H
#define DEFINES_H
#include"cuda_runtime.h"

#ifndef DOUBLE_PRECISION
#define SINGLE_PRECISION
#endif


#define fori(x,y) for(int i=x; i<int(y); i++)
#define forj(x,y) for(int j=x; j<int(y); j++)
namespace gpu{

#if defined SINGLE_PRECISION
using  real  = float;
using  real2 = float2;
using  real3 = float3;
using  real4 = float4;
#define Eps real(1e-7)
#else
using  real  = double;
using  real2 = double2;
using  real3 = double3;
using  real4 = double4;
#define Eps real(1e-16)
#endif

//id for liquid particles
#define LIQUID real(0)
//id for boundary particles
#define WALL real(1)
//id for 3rdBody particles
#define THIRDBODY real(2)

}
#endif
