#include "complex_implicit_shape.h"
#include "base_body.h"
#include "sph_system.h"

namespace SPH
{
    ComplexImplicitShape::ComplexImplicitShape(SPHBody *sph_body, ImplicitSurface &implicit_surface) 
    : LevelSetShape(sph_body, implicit_surface, false, false), implicit_surface_(implicit_surface)
    { 
        
    }


    BoundingBox ComplexImplicitShape::findBounds()
    {
		return implicit_surface_.findBounds();
    }

    Vec3d ComplexImplicitShape::findClosestPoint(const Vec3d& input_pnt)
    {
        return implicit_surface_.findClosestPoint(input_pnt);
    }

    bool ComplexImplicitShape::checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED)
    {
        return implicit_surface_.checkContain(input_pnt, BOUNDARY_INCLUDED);
    }

	bool ComplexImplicitShape::checkNotFar(const Vec3d &input_pnt, Real threshold)
    {
        return implicit_surface_.checkNotFar(input_pnt, threshold);
    }

	bool ComplexImplicitShape::checkNearSurface(const Vec3d &input_pnt, Real threshold)
    {
        return implicit_surface_.checkNearSurface(input_pnt, threshold);
    }

	Real ComplexImplicitShape::findSignedDistance(const Vec3d &input_pnt) 
    {
        return implicit_surface_.findSignedDistance(input_pnt);   
    }

	Vec3d ComplexImplicitShape::findNormalDirection(const Vec3d &input_pnt)
    {
        return implicit_surface_.findNormalDirection(input_pnt);
    }

}