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
* @version	0.1
*/
#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

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

using namespace std;

namespace SPH {

	/**
	 * @brief preclaimed classes.
	 */
	class Kernel;

	class TriangleMeshShape : public Shape
	{
	public:
		//constructor for load stl file from out side
		TriangleMeshShape(string file_path_name, Vec3d translation, Real scale_factor);
		// constructor for brick geometry
		TriangleMeshShape(Vec3d halfsize, int resolution, Vec3d translation);
		// constructor for sphere geometry
		TriangleMeshShape(Real radius, int resolution, Vec3d translation);
		//constructor for cylinder geometry
		TriangleMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation);

		SimTK::ContactGeometry::TriangleMesh* getTriangleMesh() { return triangle_mesh_; };
		bool checkContain(Vec3d pnt, bool BOUNDARY_INCLUDED = true);
		virtual Vec3d findClosestPoint(Vec3d input_pnt) override;
		virtual void findBounds(Vec3d &lower_bound, Vec3d &upper_bound) override;

	protected:
		SimTK::ContactGeometry::TriangleMesh* triangle_mesh_;

		//generate triangle mesh from polymesh
		SimTK::ContactGeometry::TriangleMesh* generateTriangleMesh(SimTK::PolygonalMesh& ploy_mesh);
	};

	class ComplexShape : public Shape
	{
	public:
		ComplexShape() : Shape("ComplexShape") {};
		ComplexShape(string complex_shape_name) : Shape(complex_shape_name) {};
		virtual ~ComplexShape() {};
		virtual void findBounds(Vec3d& lower_bound, Vec3d& upper_bound) override;

		void addTriangleMeshShape(TriangleMeshShape* triangle_mesh_shape, ShapeBooleanOps op);
		void addComplexShape(ComplexShape* complex_shape, ShapeBooleanOps op);
		void addBrick(Vec3d halfsize, int resolution, Vec3d translation, ShapeBooleanOps op);
		void addSphere(Real radius, int resolution, Vec3d translation, ShapeBooleanOps op);
		void addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, ShapeBooleanOps op);
		void addFormSTLFile(string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op);

		virtual bool checkContain(Vec3d input_pnt, bool BOUNDARY_INCLUDED = true);
		virtual bool checkNotFar(Vec3d input_pnt, Real threshold);
		virtual Vec3d findClosestPoint(Vec3d input_pnt) override;
		virtual Real findSignedDistance(Vec3d input_pnt);
		virtual Vec3d findNormalDirection(Vec3d input_pnt);
		virtual Vecd weightedIntegral(Vecd input_pnt, Kernel * kernel, Real smoothing_length) { return Vecd(1.0); };
		virtual Vecd computeKernelIntegral(Vecd input_pnt, Kernel* kernel);
	protected:
		/** shape container<pointer to geomtry, operation> */
		std::vector<std::pair<TriangleMeshShape*, ShapeBooleanOps>> triangle_mesh_shapes_;
	};
}

