/* ---------------------------------------------------------------------------*
*                       SPHinXsys: geometry                                      *
* ----------------------------------------------------------------------------*
* Geometry is used from construct a region. A region is usually composed of   *
* one or several geometries. Here, we require that the geometries fully       *
* contain or no overlap to each other only.                                   *
* ----------------------------------------------------------------------------*/

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

	class TriangleMeshShape : public Shape
	{
	public:
		//constructor for load stl file from out side
		TriangleMeshShape(string file_path_name, Vec3d translation, Real scale_factor);
		// constructor for brick geometry
		TriangleMeshShape(Vec3d halfsize, int resol, Vec3d translation);
		// constructor for sphere geometry
		TriangleMeshShape(Real radius, int resol, Vec3d translation);
		//constructor for cylinder geometry
		TriangleMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation);

		SimTK::ContactGeometry::TriangleMesh* getTriangleMesh() { return triangle_mesh_; };
		bool checkContain(Vec3d pnt, bool BOUNDARY_INCLUDED = true);
		virtual Vec3d findClosestPoint(Vec3d input_pnt) override;
		virtual void findBounds(Vec3d &lower_bound, Vec3d &upper_bound) override;
		void writePolygonalVertices(int poly_id, string out_folder);

	protected:
		SimTK::ContactGeometry::TriangleMesh* triangle_mesh_;

		//generate triangle mesh from polymesh
		SimTK::ContactGeometry::TriangleMesh* generateTriangleMesh(SimTK::PolygonalMesh& ploy_mesh);
	};

	class ComplexShape : public Shape
	{
	public:
		ComplexShape(string complex_shape_name) : Shape(complex_shape_name) {};
		virtual ~ComplexShape() {};
		bool checkContain(Vec3d pnt, bool BOUNDARY_INCLUDED = true);
		virtual Vec3d findClosestPoint(Vec3d input_pnt) override;
		virtual void findBounds(Vec3d &lower_bound, Vec3d &upper_bound) override;

		void addTriangleMeshShape(TriangleMeshShape* triangle_mesh_shape, ShapeBooleanOps op);
		void addComplexShape(ComplexShape* complex_shape, ShapeBooleanOps op);
		void addBrick(Vec3d halfsize, int resol, Vec3d translation, ShapeBooleanOps op);
		void addSphere(Real radius, int resol, Vec3d translation, ShapeBooleanOps op);
		void addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation, ShapeBooleanOps op);
		void addFormSTLFile(string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op);
		void writePolygonalVertices();

	protected:
		/** shape container<pointer to geomtry, operation> */
		std::vector<std::pair<TriangleMeshShape*, ShapeBooleanOps>> triangle_mesh_shapes_;
		string output_folder_;
	};
}

