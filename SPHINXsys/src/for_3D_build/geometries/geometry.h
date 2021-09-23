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
* @file geometry.h
* @brief Here, we define the 3D geometric algortihms. they are based on the polymesh. 
* @details The idea is to define complex geometry by passing stl, obj or other 
*			polymesh files. 
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"
#include "simbody_middle.h"
#include "complex_shape_mesh.h"
#include "complex_shape_triangle_mesh.h"
#include "complex_shape_image_mesh.h"
#include "image_mesh_shape.h"

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

	/**
	 * @brief preclaimed classes.
	 */
	class Kernel;

	class TriangleMeshShape : public Shape
	{
	public:
		//constructor for load stl file from out side
		TriangleMeshShape(std::string file_path_name, Vec3d translation, Real scale_factor);
		// constructor for brick geometry
		TriangleMeshShape(Vec3d halfsize, int resolution, Vec3d translation);
		// constructor for sphere geometry
		TriangleMeshShape(Real radius, int resolution, Vec3d translation);
		//constructor for cylinder geometry
		TriangleMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation);

		SimTK::ContactGeometry::TriangleMesh *getTriangleMesh() { return triangle_mesh_; };
		bool checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED = true);
		Vec3d findClosestPoint(const Vec3d &input_pnt);
		virtual BoundingBox findBounds() override;

	protected:
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh_;

		//generate triangle mesh from polymesh
		SimTK::ContactGeometry::TriangleMesh *generateTriangleMesh(SimTK::PolygonalMesh &ploy_mesh);
	};

	class ComplexShape : public Shape
	{
		Vec3d findClosestPoint(const Vec3d &input_pnt);

	public:
		ComplexShape() : Shape("ComplexShape") { complex_shape_mesh_ = nullptr; };
		ComplexShape(std::string complex_shape_name) : Shape(complex_shape_name) { complex_shape_mesh_ = nullptr; };
		ComplexShape(ComplexShapeMesh*complex_shape_mesh) : Shape("ComplexShape") { complex_shape_mesh_ = complex_shape_mesh; };
		virtual ~ComplexShape() {};
		virtual BoundingBox findBounds() override;

		void addTriangleMeshShape(TriangleMeshShape *triangle_mesh_shape, ShapeBooleanOps op)		
		{
			if(complex_shape_mesh_->getName() == "ComplexShapeTriangleMesh")
			{
				ComplexShapeTriangleMesh *complex_shape_mesh = dynamic_cast<ComplexShapeTriangleMesh*>(complex_shape_mesh_);
				complex_shape_mesh->addTriangleMeshShape(triangle_mesh_shape, op);
			}
		}
		// void addComplexShape(ComplexShape *complex_shape, ShapeBooleanOps op);
		void addBrick(Vec3d halfsize, int resolution, Vec3d translation, ShapeBooleanOps op)		
		{
			if(complex_shape_mesh_->getName() == "ComplexShapeTriangleMesh")
			{
				ComplexShapeTriangleMesh *complex_shape_mesh = dynamic_cast<ComplexShapeTriangleMesh*>(complex_shape_mesh_);
				complex_shape_mesh->addBrick(halfsize, resolution, translation, op);
			}			
		}
		//void addSphere(Real radius, int resolution, Vec3d translation, ShapeBooleanOps op);
		//void addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, ShapeBooleanOps op);
		//void addFormSTLFile(std::string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op);

		virtual bool checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED = true);
		virtual bool checkNotFar(const Vec3d &input_pnt, Real threshold);
		virtual bool checkNearSurface(const Vec3d &input_pnt, Real threshold);
		/** Signed distance is negative for point within the complex shape. */
		virtual Real findSignedDistance(const Vec3d &input_pnt);
		/** Normal direction point toward outside of the complex shape. */
		virtual Vec3d findNormalDirection(const Vec3d &input_pnt);

	protected:
		/** shape container<pointer to geomtry, operation> */

		ComplexShapeMesh* complex_shape_mesh_;
	};
}

#endif //GEOMETRY_3D_H
