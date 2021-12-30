/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/

/**
* @file complex_implicit_surface.h
* @brief Here, we define a data structure to be used as the set of implicit surface
* @author	Huiqiang Yue
*/

#ifndef __COMPLEX_IMPLICIT_SHAPE_H__
#define __COMPLEX_IMPLICIT_SHAPE_H__

#include "implicit_surface.h"
#include "base_geometry.h"
#include "level_set_shape.h"

#include <memory>

namespace SPH
{
    /**
     * Implicit shape defined by level set. This class is just a simple wapper for ImplicitSurface class in order to make 
     * the class compatible with other data structures.
     */
    class ComplexImplicitShape : public LevelSetShape
    {
    public:
        ComplexImplicitShape(SPHBody *sph_body, ImplicitSurface &implicit_surface);
        virtual ~ComplexImplicitShape(){};

        virtual BoundingBox findBounds() override;
        virtual Vec3d findClosestPoint(const Vec3d& input_pnt) override;
 
        virtual bool checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual bool checkNotFar(const Vec3d &input_pnt, Real threshold) override;
		virtual bool checkNearSurface(const Vec3d &input_pnt, Real threshold) override;
		virtual Real findSignedDistance(const Vec3d &input_pnt) override;
		virtual Vec3d findNormalDirection(const Vec3d &input_pnt) override;

    private:

        ImplicitSurface implicit_surface_;

    };
}


#endif // __COMPLEX_IMPLICIT_SHAPE_H__