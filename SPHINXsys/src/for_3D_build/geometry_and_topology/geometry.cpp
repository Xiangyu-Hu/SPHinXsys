#include "geometry.h"

using namespace std;

namespace SPH 
{
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(string filepathname, Vec3d translation, Real scale_factor)
		: Shape("TriangleMeshShape")
	{
		if (!fs::exists(filepathname))
		{
			std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlFile(filepathname);
		polymesh.scaleMesh(scale_factor);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(Vec3d halfsize, int resol, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createBrickMesh(halfsize, resol);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(Real radius, int resol, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh 
			= SimTK::PolygonalMesh::createSphereMesh(radius, resol);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh 
			= SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resol);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	SimTK::ContactGeometry::TriangleMesh* TriangleMeshShape
		::generateTriangleMesh(SimTK::PolygonalMesh &ploy_mesh)
	{
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
		triangle_mesh = new SimTK::ContactGeometry::TriangleMesh(ploy_mesh);
		if (!SimTK::ContactGeometry::TriangleMesh::isInstance(*triangle_mesh)) 
		{
			std::cout << "\n Error the triangle mesh is not valid" << std::endl;
		}
		std::cout << "num of faces:" << triangle_mesh->getNumFaces() << std::endl;

		return triangle_mesh;
	}
	//=================================================================================================//
	bool TriangleMeshShape::checkContain(Vec3d pnt, bool BOUNDARY_INCLUDED)
	{

		SimTK::Vec2 uv_coordinate;
		bool inside = false;
		int face_id;
		Vec3d closest_pnt = triangle_mesh_->findNearestPoint(pnt, inside, face_id, uv_coordinate);

		vector<int> neigbor_face(4);
		neigbor_face[0] = face_id;
		/** go throught the neighbor faces. */
		for (int i = 1; i < 4; i++) {
			int edge = triangle_mesh_->getFaceEdge(face_id, i - 1);
			int face = triangle_mesh_->getEdgeFace(edge, 0);
			neigbor_face[i] = face != face_id ? face : triangle_mesh_->getEdgeFace(edge, 1);
		}

		Vec3d from_face_to_pnt = pnt - closest_pnt;
		Real sum_weights = 0.0;
		Real weigthed_dot_product = 0.0;
		for (int i = 0; i < 4; i++) {
			SimTK::UnitVec3 normal_direction = triangle_mesh_->getFaceNormal(neigbor_face[i]);
			Real dot_product = dot(normal_direction, from_face_to_pnt);
			Real weight = dot_product * dot_product;
			weigthed_dot_product += weight * dot_product;
			sum_weights += weight;
		}

		weigthed_dot_product /= sum_weights;

		bool weighted_inside = false;
		if (weigthed_dot_product < 0.0) weighted_inside = true;

		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		return weighted_inside;
	}
	//=================================================================================================//
	Vec3d TriangleMeshShape::findClosestPoint(Vec3d input_pnt)
	{
		bool inside = false;
		int face_id;
		SimTK::Vec2 normal;
		Vec3d closest_pnt;
		closest_pnt = triangle_mesh_->findNearestPoint(input_pnt, inside,face_id, normal);
		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces()) 
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		return closest_pnt;
	}
	//=================================================================================================//
	void TriangleMeshShape::findBounds(Vec3d &lower_bound, Vec3d &upper_bound)
	{
		int number_of_vertices = triangle_mesh_->getNumVertices();
		//initial reference values
		lower_bound = Vec3d(Infinity);
		upper_bound = Vec3d(-Infinity);

		for (int i = 0; i != number_of_vertices; ++i)
		{
			Vec3d vertex_position = triangle_mesh_->getVertexPosition(i);
			for (int j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			}
		}
	}
	//=================================================================================================//
	bool ComplexShape::checkContain(Vec3d pnt, bool BOUNDARY_INCLUDED)
	{
		bool exist = false;
		bool inside = false;
		
		for (auto& Rshape : triangle_mesh_shapes_)
		{
			ShapeBooleanOps opstring = Rshape.second;
			TriangleMeshShape* sp = Rshape.first;

			switch (opstring)
			{
			case ShapeBooleanOps::add:
			{
				inside = sp->checkContain(pnt);
				exist = exist || inside;
				break;
			}
			case ShapeBooleanOps::sub:
			{
				inside = sp->checkContain(pnt);
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
	//=================================================================================================//
	Vec3d ComplexShape::findClosestPoint(Vec3d input_pnt)
	{
		//a big positive number
		Real large_number(Infinity);
		Real dist_min = large_number;
		Vec3d pnt_closest(0);
		Vec3d pnt_found(0);

		for (auto& Rshape : triangle_mesh_shapes_)
		{
			TriangleMeshShape* sp = Rshape.first;
			pnt_found  = sp->findClosestPoint(input_pnt);
			Real dist = (input_pnt - pnt_found).norm();

			if(dist <= dist_min)
			{
				dist_min = dist;
				pnt_closest = pnt_found;
			}
		}

		return pnt_closest;
	}
	//=================================================================================================//
	void ComplexShape::addTriangleMeshShape(TriangleMeshShape* triangle_mesh_shape, ShapeBooleanOps op)
	{
		pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addBrick(Vec3d halfsize, int resol, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(halfsize, resol, translation);
		pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addSphere(Real radius, int resol, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(radius, resol, translation);
		pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resol, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(axis, radius, halflength, resol, translation);
		pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addFormSTLFile(string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(file_path_name, translation, scale_factor);
		pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addComplexShape(ComplexShape* complex_shape, ShapeBooleanOps op)
	{
		switch (op)
		{
		case ShapeBooleanOps::add:
		{
			for (auto& shape_and_op : complex_shape->triangle_mesh_shapes_)
			{
				triangle_mesh_shapes_.push_back(shape_and_op);
			}
			break;
		}
		case ShapeBooleanOps::sub:
		{
			for (auto& shape_and_op : complex_shape->triangle_mesh_shapes_)
			{
				TriangleMeshShape* sp = shape_and_op.first;
				ShapeBooleanOps opstring 
					= shape_and_op.second == ShapeBooleanOps::add ? ShapeBooleanOps::sub : ShapeBooleanOps::add;
				pair<TriangleMeshShape*, ShapeBooleanOps> substract_shape_and_op(sp, opstring);

				triangle_mesh_shapes_.push_back(substract_shape_and_op);
			}
			break;
		}
		}
	}
	//=================================================================================================//
	void ComplexShape::findBounds(Vec3d &lower_bound, Vec3d &upper_bound)
	{
		//initial reference values
		lower_bound = Vec3d(Infinity);
		upper_bound = Vec3d(-Infinity);

		for (int i = 0; i < triangle_mesh_shapes_.size(); i++)
		{
			Vec3d shape_lower_bound(0), shape_upper_bound(0);
			triangle_mesh_shapes_[i].first->findBounds(shape_lower_bound, shape_upper_bound);
			for (int j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], shape_lower_bound[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_upper_bound[j]);
			}
		}
	}
	//=================================================================================================//
}