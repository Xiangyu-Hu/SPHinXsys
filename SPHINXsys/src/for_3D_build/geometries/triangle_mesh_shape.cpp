#include "triangle_mesh_shape.h"

namespace SPH
{
	//=================================================================================================//
	SimTK::ContactGeometry::TriangleMesh *TriangleMeshShape::generateTriangleMesh(const SimTK::PolygonalMesh &poly_mesh)
	{
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
		triangle_mesh = triangle_mesh_ptr_keeper_.createPtr<SimTK::ContactGeometry::TriangleMesh>(poly_mesh);
		if (!SimTK::ContactGeometry::TriangleMesh::isInstance(*triangle_mesh))
		{
			std::cout << "\n Error the triangle mesh is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		std::cout << "num of faces:" << triangle_mesh->getNumFaces() << std::endl;

		return triangle_mesh;
	}
	//=================================================================================================//
	bool TriangleMeshShape::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
	{

		SimTK::Vec3 in_pnt(pnt[0], pnt[1], pnt[2]);
		SimTK::Vec2 uv_coordinate;
		bool inside = false;
		int face_id;
		SimTK::Vec3 closest_pnt = triangle_mesh_->findNearestPoint(in_pnt, inside, face_id, uv_coordinate);

		StdVec<int> neighbor_face(4);
		neighbor_face[0] = face_id;
		/** go through the neighbor faces. */
		for (int i = 1; i < 4; i++)
		{
			int edge = triangle_mesh_->getFaceEdge(face_id, i - 1);
			int face = triangle_mesh_->getEdgeFace(edge, 0);
			neighbor_face[i] = face != face_id ? face : triangle_mesh_->getEdgeFace(edge, 1);
		}

		SimTK::Vec3 from_face_to_pnt = in_pnt - closest_pnt;
		Real sum_weights = 0.0;
		Real weighted_dot_product = 0.0;
		for (int i = 0; i < 4; i++)
		{
			SimTK::UnitVec3 normal_direction = triangle_mesh_->getFaceNormal(neighbor_face[i]);
			Real dot_product = dot(normal_direction, from_face_to_pnt);
			Real weight = dot_product * dot_product;
			weighted_dot_product += weight * dot_product;
			sum_weights += weight;
		}

		weighted_dot_product /= sum_weights;

		bool weighted_inside = false;
		if (weighted_dot_product < 0.0)
			weighted_inside = true;

		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}

		return weighted_inside;
	}
	//=================================================================================================//
	Vecd TriangleMeshShape::findClosestPoint(const Vecd &probe_point)
	{
		bool inside = false;
		int face_id;
		SimTK::Vec2 norm;
		SimTK::Vec3 closest_pnt = triangle_mesh_->findNearestPoint(SimTK::Vec3(probe_point[0], probe_point[1], probe_point[2]), inside, face_id, norm);
		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		return Vecd(closest_pnt[0], closest_pnt[1], closest_pnt[2]);
	}
	//=================================================================================================//
	BoundingBox TriangleMeshShape::findBounds()
	{
		int number_of_vertices = triangle_mesh_->getNumVertices();
		//initial reference values
		SimTK::Vec3 lower_bound = SimTK::Vec3(Infinity);
		SimTK::Vec3 upper_bound = SimTK::Vec3(-Infinity);

		for (int i = 0; i != number_of_vertices; ++i)
		{
			SimTK::Vec3 vertex_position = triangle_mesh_->getVertexPosition(i);
			for (int j = 0; j != 3; ++j)
			{
				lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			}
		}
		return BoundingBox(Vecd(lower_bound[0],lower_bound[1],lower_bound[2]), Vecd(upper_bound[0],upper_bound[1], upper_bound[2]));
	}
	//=================================================================================================//
	TriangleMeshShapeSTL::TriangleMeshShapeSTL(const std::string &filepathname, Vecd translation, Real scale_factor,
			const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		if (!fs::exists(filepathname))
		{
			std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlFile(filepathname);
		polymesh.scaleMesh(scale_factor);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTK::Vec3(translation[0], translation[1], translation[2])));
	}
	//=================================================================================================//
	TriangleMeshShapeSTL::TriangleMeshShapeSTL(const std::string &filepathname, Mat3d rotation,
			Vecd translation, Real scale_factor, const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		if (!fs::exists(filepathname))
		{
			std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlFile(filepathname);

        polymesh.scaleMesh(scale_factor);
        SimTK::Transform_<Real> transform( SimTK::Rotation_<Real>(SimTK::Mat33(rotation(0,0), rotation(0,1), rotation(0,2), 
						 							   						   rotation(1,0), rotation(1,1), rotation(1,2), 
						 							   						   rotation(2,0), rotation(2,1), rotation(2,2))), 
									SimTK::Vec3(translation[0], translation[1], translation[2]) );
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(transform));
	}
	//=================================================================================================//
	#ifdef __EMSCRIPTEN__	
	TriangleMeshShapeSTL::TriangleMeshShapeSTL(const uint8_t* buffer, Vecd translation, Real scale_factor, 
			const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlBuffer(buffer);
		polymesh.scaleMesh(scale_factor);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(SimTK::Vec3(translation[0], translation[1], translation[2])));
	}
	#endif
	//=================================================================================================//
	TriangleMeshShapeBrick::TriangleMeshShapeBrick(Vecd halfsize, int resolution, Vecd translation,
			const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createBrickMesh(SimTK::Vec3(halfsize[0], halfsize[1], halfsize[2]), resolution);
		triangle_mesh_ = generateTriangleMesh( polymesh.transformMesh(SimTK::Vec3(translation[0], translation[1], translation[2])) );
	}
	//=================================================================================================//
	TriangleMeshShapeBrick::TriangleMeshShapeBrick(const TriangleMeshShapeBrick::ShapeParameters &shape_parameters,
			const std::string &shape_name)
		: TriangleMeshShapeBrick(shape_parameters.halfsize_, shape_parameters.resolution_,
								 shape_parameters.translation_, shape_name) {}
	//=================================================================================================//
	TriangleMeshShapeSphere::TriangleMeshShapeSphere(Real radius, int resolution, Vecd translation,
			const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createSphereMesh(radius, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh( SimTK::Vec3(translation[0], translation[1], translation[2]) ));
	}
	//=================================================================================================//
	TriangleMeshShapeCylinder::TriangleMeshShapeCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution,
			Vecd translation, const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh =
			SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh( SimTK::Vec3(translation[0], translation[1], translation[2]) ));
	}
	//=================================================================================================//
}
