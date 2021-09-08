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
* @file complex_shape_triangle_mesh.h
* @brief x 
* @details x 
*			x
* @author	Yijin Mao
*/

#ifndef COMPLEX_SHAPE_TRIANGLE_MESH_H
#define COMPLEX_SHAPE_TRIANGLE_MESH_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"
#include "simbody_middle.h"
#include "complex_shape_mesh.h"

#include <iostream>
#include <string>
#include <fstream>

/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH
{
	class TriangleMeshShape;

	class ComplexShapeTriangleMesh : public ComplexShapeMesh
	{
	public:
		ComplexShapeTriangleMesh() : ComplexShapeMesh("ComplexShapeTriangleMesh") {};
		ComplexShapeTriangleMesh(std::string complex_shape_name) : ComplexShapeMesh(complex_shape_name) {};
		virtual ~ComplexShapeTriangleMesh() {};
		virtual BoundingBox findBounds() override;

		void addTriangleMeshShape(TriangleMeshShape* triangle_mesh_shape, ShapeBooleanOps op); 
		void addComplexShapeTriangleMesh(ComplexShapeTriangleMesh* complex_shape_triangle_mesh, ShapeBooleanOps op);
		virtual void addBrick(Vec3d halfsize, int resolution, Vec3d translation, ShapeBooleanOps op);
		virtual void addSphere(Real radius, int resolution, Vec3d translation, ShapeBooleanOps op);
		virtual void addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, ShapeBooleanOps op);
		void addFormSTLFile(std::string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op);

		virtual bool checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED = true);
		virtual bool checkNotFar(const Vec3d& input_pnt, Real threshold);
		virtual bool checkNearSurface(const Vec3d& input_pnt, Real threshold);
		/** Signed distance is negative for point within the complex shape. */
		virtual Real findSignedDistance(const Vec3d& input_pnt);
		/** Normal direction point toward outside of the complex shape. */
		virtual Vec3d findNormalDirection(const Vec3d& input_pnt);
		virtual Vec3d findClosestPoint(const Vec3d& input_pnt);
	protected:
		/** shape container<pointer to geomtry, operation> */
		std::vector<std::pair<TriangleMeshShape*, ShapeBooleanOps>> triangle_mesh_shapes_;
	};
}

#endif //COMPLEX_SHAPE_TRIANGLE_MESH_H
