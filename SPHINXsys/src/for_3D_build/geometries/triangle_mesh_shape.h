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
* @file triangle_mesh_shape.h
* @brief Here, we define the 3D geometric algortihms. they are based on the polymesh. 
* @details The idea is to define complex geometry by passing stl, obj or other 
*			polymesh files. 
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef TRIANGULAR_MESH_SHAPE_H
#define TRIANGULAR_MESH_SHAPE_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "simbody_middle.h"

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
	class TriangleMeshShape : public Shape
	{
	private:
		UniquePtrKeeper<SimTK::ContactGeometry::TriangleMesh> triangle_mesh_ptr_keeper_;

	public:
		explicit TriangleMeshShape(const std::string &shape_name)
			: Shape(shape_name), triangle_mesh_(nullptr){};

		virtual bool checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec3d findClosestPoint(const Vec3d &input_pnt) override;
		virtual BoundingBox findBounds() override;

		SimTK::ContactGeometry::TriangleMesh *getTriangleMesh() { return triangle_mesh_; };

	protected:
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh_;

		//generate triangle mesh from polymesh
		SimTK::ContactGeometry::TriangleMesh *generateTriangleMesh(SimTK::PolygonalMesh &ploy_mesh);
	};

	class TriangleMeshShapeSTL : public TriangleMeshShape
	{
	public:
		//constructor for load STL file from out side
		explicit TriangleMeshShapeSTL(const std::string &file_path_name, Vec3d translation, Real scale_factor,
									  const std::string &shape_name = "TriangleMeshShapeSTL");
		virtual ~TriangleMeshShapeSTL(){};
	};

	class TriangleMeshShapeBrick : public TriangleMeshShape
	{
	public:
		class ShapeParameters 
		{
		public:
			ShapeParameters() : halfsize_(0), resolution_(0), translation_(0) {};
			Vec3d halfsize_;
			int resolution_;
			Vec3d translation_;
		};
		explicit TriangleMeshShapeBrick(Vec3d halfsize, int resolution, Vec3d translation,
										const std::string &shape_name = "TriangleMeshShapeBrick");
		explicit TriangleMeshShapeBrick(const TriangleMeshShapeBrick::ShapeParameters &shape_parameters,
										const std::string &shape_name = "TriangleMeshShapeBrick");
		virtual ~TriangleMeshShapeBrick(){};
	};

	class TriangleMeshShapeShere : public TriangleMeshShape
	{
	public:
		//constructor for sphere shape
		explicit TriangleMeshShapeShere(Real radius, int resolution, Vec3d translation,
										const std::string &shape_name = "TriangleMeshShapeShere");
		virtual ~TriangleMeshShapeShere(){};
	};

	class TriangleMeshShapeCylinder : public TriangleMeshShape
	{
	public:
		//constructor for cylinder shape
		explicit TriangleMeshShapeCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation,
										   const std::string &shape_name = "TriangleMeshShapeCylinder");
		virtual ~TriangleMeshShapeCylinder(){};
	};
}

#endif //TRIANGULAR_MESH_SHAPE_H
