#include "vector.h"
//=================================================================================================//
namespace SPH
{    
   //=================================================================================================//
     Vec3d OcctToEigen(const gp_Pnt &occt_vector)
   {
    return Vec3d(occt_vector.X(), occt_vector.Y(), occt_vector.Z());
   }
   //=================================================================================================//
     gp_Pnt EigenToOcct(const Vec3d &eigen_vector)
    {
    return gp_Pnt(eigen_vector[0], eigen_vector[1], eigen_vector[2]);
    }
   //=================================================================================================//
   Vec3d OcctVecToEigen(const gp_Vec &occt_vector)
   {
    return Vec3d(occt_vector.X(), occt_vector.Y(), occt_vector.Z());
   }
	//=================================================================================================//

}
