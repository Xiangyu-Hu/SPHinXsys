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
#include "array.h"
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

	enum class RegionBooleanOps { add, sub };

	class Geometry  : public Shape
	{
	protected:
		//generate trangle mesh from polymesh
		SimTK::ContactGeometry::TriangleMesh* TriangleMeshFromPolyMesh(SimTK::PolygonalMesh &ploy_mesh);
	public:
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh_;

		//constructor for load stl file from out side
		Geometry(string file_path_name, Vec3d translation, Real scale_factor);
		// constructor for brick geometry
		Geometry(Vec3d halfsize, int resol, Vec3d translation); 
		// constructor for sphere geometry
		Geometry(Real radius, int resol, Vec3d translation); 
		//constructor for cylinder geometry
		Geometry(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation); 

		virtual bool contain(Vec3d pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec3d closestpointonface(Vec3d input_pnt) override;
		virtual void shapebound(Vec3d &lower_bound, Vec3d &upper_bound) override;
		void writePolygonalVertices(int poly_id, string out_folder);
	};

	class Region
	{
	protected:
		//name of the region
		std::string region_name_; 											
		//geometry container
		std::vector<Geometry*> geometries;
		//geometry operation container
		std::vector<RegionBooleanOps> geometryops;							
		//shapes container<pointer to geomerty, operation>
		std::vector<std::pair<Geometry*, RegionBooleanOps>> shapes;
		string output_folder_;
	public:
		Region(string region_name);
		virtual ~Region() {};
		virtual bool contain(Vec3d pnt, bool BOUNDARY_INCLUDED = true);
		virtual void closestpointonface(Vec3d input_pnt, Vec3d& closest_pnt, Real& phi);
		virtual void regionbound(Vec3d &lower_bound, Vec3d &upper_bound);

		void add_geometry(Geometry* geometry, RegionBooleanOps op);
		void add_region(Region *region, RegionBooleanOps op);
		void add_brick(Vec3d halfsize, int resol, Vec3d translation, RegionBooleanOps op);
		void add_sphere(Real radius, int resol, Vec3d translation, RegionBooleanOps op);
		void add_cylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation, RegionBooleanOps op);
		void add_from_STL_file(string file_path_name, Vec3d translation, Real scale_factor, RegionBooleanOps op);

		void done_modeling();
	};
}

