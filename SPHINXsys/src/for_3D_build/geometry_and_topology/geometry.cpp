#include "geometry.h"

using namespace std;

namespace SPH {

	Geometry::Geometry(string filepathname, Vec3d translation) 
		: Shape("Geoemetry")
	{
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlFile(filepathname);
		triangle_mesh_ = TriangleMeshFromPolyMesh(polymesh.transformMesh(translation));
	}

	Geometry::Geometry(Vec3d halfsize, int resol, Vec3d translation) 
		: Shape("Geoemetry")
	{
		SimTK::PolygonalMesh polymesh = PolygonalMesh::createBrickMesh(halfsize, resol);
		triangle_mesh_ = TriangleMeshFromPolyMesh(polymesh.transformMesh(translation));
	}

	Geometry::Geometry(Real radius, int resol, Vec3d translation) 
		: Shape("Geoemetry")
	{
		SimTK::PolygonalMesh polymesh 
			= PolygonalMesh::createSphereMesh(radius, resol);
		triangle_mesh_ = TriangleMeshFromPolyMesh(polymesh.transformMesh(translation));
	}

	Geometry::Geometry(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation) 
		: Shape("Geoemetry")
	{
		SimTK::PolygonalMesh polymesh 
			= PolygonalMesh::createCylinderMesh(axis, radius, halflength, resol);
		triangle_mesh_ = TriangleMeshFromPolyMesh(polymesh.transformMesh(translation));
	}

	SimTK::ContactGeometry::TriangleMesh* Geometry
		::TriangleMeshFromPolyMesh(SimTK::PolygonalMesh &ploy_mesh)
	{
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
		triangle_mesh = new SimTK::ContactGeometry::TriangleMesh(ploy_mesh);
		if (!SimTK::ContactGeometry::TriangleMesh::isInstance(*triangle_mesh)) {
			std::cout << "\n Error the triangle mesh is not valid" << std::endl;
		}
		std::cout << "num of faces:" << triangle_mesh->getNumFaces() << std::endl;

		return triangle_mesh;
	}

	bool Geometry::contain(Vec3d pnt, bool BOUNDARY_INCLUDED /*= true*/)
	{
		bool inside = false;
		SimTK::UnitVec3 normal;
		Vec3d nearest = triangle_mesh_->findNearestPoint(pnt, inside, normal);

		return inside;
	}

	Vec3d Geometry::closestpointonface(Vec3d input_pnt)
	{
		bool inside = false;
		SimTK::UnitVec3 normal;
		Vec3d nearest = triangle_mesh_->findNearestPoint(input_pnt, inside, normal);

		return input_pnt - nearest;
	}

	void Geometry::shapebound(Vec3d &lower_bound, Vec3d &upper_bound)
	{
		size_t number_of_vertices = triangle_mesh_->getNumVertices();
		//initial reference values
		lower_bound = Vec3d(1.0e15);
		upper_bound = Vec3d(-1.0e15);

		for (size_t i = 0; i != number_of_vertices; ++i)
		{
			Vec3d vertex_position = triangle_mesh_->getVertexPosition(i);
			for (size_t j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			}
		}
	}

	Region::Region(string region_name)
	{
		region_name_ = region_name;
	}

	bool Region::contain(Vec3d pnt, bool BOUNDARY_INCLUDED /*= true*/)
	{
		bool exist = false;
		bool inside = false;
		SimTK::UnitVec3 normal;
		for (auto& Rshape : shapes)
		{
			RegionBooleanOps opstring = Rshape.second;
			Geometry* sp = Rshape.first;

			switch (opstring)
			{
			case RegionBooleanOps::add:
			{
				inside = sp->contain(pnt);
				exist = exist || inside;
				break;
			}
			case RegionBooleanOps::sub:
			{
				inside = sp->contain(pnt);
				exist = exist && (!inside);
				break;
			}
			default:
			{
				std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
				break;
			}
			}
		}
		return exist;
	}

	void Region::closestpointonface(Vec3d input_pnt, Vec3d& closest_pnt, Real& phi)
	{
		bool exist = false;
		bool inside = false;
		SimTK::UnitVec3 normal;
	
		//a big positive number
		Real large_number(1.0e15);
		phi = large_number;
		Vec3d p_find(0);

		for (auto& Rshape : shapes)
		{
			Geometry* sp = Rshape.first;

			Vec3d dispalcement = sp->closestpointonface(input_pnt);
			Real distance = dispalcement.norm();
			if (distance < phi)
			{
				phi = distance;
				p_find = dispalcement;
			}

		}

		phi = contain(input_pnt) ? phi : -phi;
		closest_pnt = p_find;
	}

	void Region::add_geometry(Geometry *geometry, RegionBooleanOps op)
	{
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}

	void Region::add_brick(Vec3d halfsize, int resol, Vec3d translation, RegionBooleanOps op)
	{
		Geometry *geometry = new Geometry(halfsize, resol, translation);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}

	void Region::add_sphere(Real radius, int resol, Vec3d translation, RegionBooleanOps op)
	{
		Geometry *geometry = new Geometry(radius, resol, translation);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}

	void Region::add_cylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation, RegionBooleanOps op)
	{
		Geometry *geometry = new Geometry(axis, radius, halflength, resol, translation);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}

	void Region::add_from_STL_file(string file_path_name, Vec3d translation, RegionBooleanOps op)
	{
		Geometry *geometry = new Geometry(file_path_name, translation);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}

	void Region::add_region(Region *region, RegionBooleanOps op)
	{
		switch (op)
		{
		case RegionBooleanOps::add:
		{
			for (auto& Rshape : region->shapes)
			{
				Geometry* sp = Rshape.first;
				RegionBooleanOps opstring = Rshape.second;

				geometries.push_back(sp);
				geometryops.push_back(opstring);
			}
			break;
		}
		case RegionBooleanOps::sub:
		{
			for (auto& Rshape : region->shapes)
			{
				Geometry* sp = Rshape.first;
				RegionBooleanOps opstring = Rshape.second;

				geometries.push_back(sp);
				opstring == RegionBooleanOps::add ?
				geometryops.push_back(RegionBooleanOps::sub) 
					: geometryops.push_back(RegionBooleanOps::add);
			}
			break;
		}
		}
	}

	void Region::done_modeling()
	{
		int faultyShapes = 0;
		bool flag = false;
		for (int i = 0; i < geometries.size(); i++)
		{
			if (flag == false)
			{
				if (geometryops[i] == RegionBooleanOps::add)
				{
					flag = true;
				}
				else
				{
					faultyShapes++;
				}
			}

			if (flag == true)
			{
				shapes.push_back(std::pair<Geometry*, RegionBooleanOps>(geometries[i], geometryops[i]));
			}
		}
		if (faultyShapes > 0) 
		{
			std::cout << "Warning: The first " << faultyShapes << " sub-shapes used to construct reagion '" << region_name_
				<< "' will be neglected as their boolean operations are BodyBooleanOps::sub!" << std::endl;
		}
	}

	void Region::regionbound(Vec3d &lower_bound, Vec3d &upper_bound)
	{
		//initial reference values
		lower_bound = Vec3d(1.0e15);
		upper_bound = Vec3d(-1.0e15);

		for (int i = 0; i < geometries.size(); i++)
		{
			Vec3d shape_lower_bound(0), shape_upper_bound(0);
			geometries[i]->shapebound(shape_lower_bound, shape_upper_bound);
			for (size_t j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], shape_lower_bound[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_upper_bound[j]);
			}
		}
	}
}