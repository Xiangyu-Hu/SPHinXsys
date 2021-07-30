#include "geometry.h"

namespace SPH
{
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(std::string filepathname, Vec3d translation, Real scale_factor)
		: Shape("TriangleMeshShape")
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
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(Vec3d halfsize, int resolution, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createBrickMesh(halfsize, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::TriangleMeshShape(Real radius, int resolution, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createSphereMesh(radius, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShape::
		TriangleMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation)
		: Shape("TriangleMeshShape")
	{
		SimTK::PolygonalMesh polymesh =
			SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	SimTK::ContactGeometry::TriangleMesh *TriangleMeshShape::
		generateTriangleMesh(SimTK::PolygonalMesh &ploy_mesh)
	{
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
		triangle_mesh = new SimTK::ContactGeometry::TriangleMesh(ploy_mesh);
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
	bool TriangleMeshShape::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
	{

		SimTK::Vec2 uv_coordinate;
		bool inside = false;
		int face_id;
		Vec3d closest_pnt = triangle_mesh_->findNearestPoint(pnt, inside, face_id, uv_coordinate);

		StdVec<int> neighbor_face(4);
		neighbor_face[0] = face_id;
		/** go throught the neighbor faces. */
		for (int i = 1; i < 4; i++)
		{
			int edge = triangle_mesh_->getFaceEdge(face_id, i - 1);
			int face = triangle_mesh_->getEdgeFace(edge, 0);
			neighbor_face[i] = face != face_id ? face : triangle_mesh_->getEdgeFace(edge, 1);
		}

		Vec3d from_face_to_pnt = pnt - closest_pnt;
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
	Vec3d TriangleMeshShape::findClosestPoint(const Vec3d &input_pnt)
	{
		bool inside = false;
		int face_id;
		SimTK::Vec2 normal;
		Vec3d closest_pnt;
		closest_pnt = triangle_mesh_->findNearestPoint(input_pnt, inside, face_id, normal);
		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		return closest_pnt;
	}
	//=================================================================================================//
	BoundingBox TriangleMeshShape::findBounds()
	{
		int number_of_vertices = triangle_mesh_->getNumVertices();
		//initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);

		for (int i = 0; i != number_of_vertices; ++i)
		{
			Vec3d vertex_position = triangle_mesh_->getVertexPosition(i);
			for (int j = 0; j != 3; ++j)
			{
				lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
	bool ComplexShape::checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED)
	{
		bool exist = false;
		bool inside = false;

		for (auto &each_shape : triangle_mesh_shapes_)
		{
			TriangleMeshShape *sp = each_shape.first;
			ShapeBooleanOps operation_string = each_shape.second;

			switch (operation_string)
			{
			case ShapeBooleanOps::add:
			{
				inside = sp->checkContain(input_pnt);
				exist = exist || inside;
				break;
			}
			case ShapeBooleanOps::sub:
			{
				inside = sp->checkContain(input_pnt);
				exist = exist && (!inside);
				break;
			}
			default:
			{
				std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				throw;
			}
			}
		}
		return exist;
	}
	//=================================================================================================//
	Vec3d ComplexShape::findClosestPoint(const Vec3d &input_pnt)
	{
		//a big positive number
		Real large_number(Infinity);
		Real dist_min = large_number;
		Vec3d pnt_closest(0);
		Vec3d pnt_found(0);

		for (auto &each_shape : triangle_mesh_shapes_)
		{
			TriangleMeshShape *sp = each_shape.first;
			pnt_found = sp->findClosestPoint(input_pnt);
			Real dist = (input_pnt - pnt_found).norm();

			if (dist <= dist_min)
			{
				dist_min = dist;
				pnt_closest = pnt_found;
			}
		}

		return pnt_closest;
	}
	//=================================================================================================//
	Real ComplexShape::findSignedDistance(const Vec3d &input_pnt)
	{
		Real distance_to_surface = (input_pnt - findClosestPoint(input_pnt)).norm();
		return checkContain(input_pnt) ? -distance_to_surface : distance_to_surface;
	}
	//=================================================================================================//
	Vec3d ComplexShape::findNormalDirection(const Vec3d &input_pnt)
	{
		bool is_contain = checkContain(input_pnt);
		Vecd displacement_to_surface = findClosestPoint(input_pnt) - input_pnt;
		while (displacement_to_surface.norm() < Eps)
		{
			Vecd jittered = input_pnt; //jittering
			for (int l = 0; l != input_pnt.size(); ++l)
				jittered[l] = input_pnt[l] + (((Real)rand() / (RAND_MAX)) - 0.5) * 100.0 * Eps;
			if (checkContain(jittered) == is_contain)
				displacement_to_surface = findClosestPoint(jittered) - jittered;
		}
		Vecd direction_to_surface = displacement_to_surface.normalize();
		return is_contain ? direction_to_surface : -1.0 * direction_to_surface;
	}
	//=================================================================================================//
	void ComplexShape::addTriangleMeshShape(TriangleMeshShape *triangle_mesh_shape, ShapeBooleanOps op)
	{
		std::pair<TriangleMeshShape *, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addBrick(Vec3d halfsize, int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape *triangle_mesh_shape = new TriangleMeshShape(halfsize, resolution, translation);
		std::pair<TriangleMeshShape *, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addSphere(Real radius, int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape *triangle_mesh_shape = new TriangleMeshShape(radius, resolution, translation);
		std::pair<TriangleMeshShape *, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength,
								   int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		/** Here SimTK::UnitVec3 give the direction of the cylinder, viz. SimTK::UnitVec3(0,0,1) create a cylinder in z-axis.*/
		TriangleMeshShape *triangle_mesh_shape =
			new TriangleMeshShape(axis, radius, halflength, resolution, translation);
		std::pair<TriangleMeshShape *, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addFormSTLFile(std::string file_path_name, Vec3d translation,
									  Real scale_factor, ShapeBooleanOps op)
	{
		TriangleMeshShape *triangle_mesh_shape =
			new TriangleMeshShape(file_path_name, translation, scale_factor);
		std::pair<TriangleMeshShape *, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShape::addComplexShape(ComplexShape *complex_shape, ShapeBooleanOps op)
	{
		switch (op)
		{
		case ShapeBooleanOps::add:
		{
			for (auto &shape_and_op : complex_shape->triangle_mesh_shapes_)
			{
				triangle_mesh_shapes_.push_back(shape_and_op);
			}
			break;
		}
		case ShapeBooleanOps::sub:
		{
			for (auto &shape_and_op : complex_shape->triangle_mesh_shapes_)
			{
				TriangleMeshShape *sp = shape_and_op.first;
				ShapeBooleanOps operation_string =
					shape_and_op.second == ShapeBooleanOps::add ? ShapeBooleanOps::sub : ShapeBooleanOps::add;
				std::pair<TriangleMeshShape *, ShapeBooleanOps> substract_shape_and_op(sp, operation_string);
				triangle_mesh_shapes_.push_back(substract_shape_and_op);
			}
			break;
		}
		case ShapeBooleanOps::sym_diff:
		{
			std::cout << "\n FAILURE: the boolean operation: ShapeBooleanOps::sym_diff is not applicable for 3D geometry!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		case ShapeBooleanOps::intersect:
		{
			std::cout << "\n FAILURE: the boolean operation: ShapeBooleanOps::intersect is not applicable for 3D geometry!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		default:
		{
			std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		}
	}
	//=================================================================================================//
	bool ComplexShape::checkNotFar(const Vec3d &input_pnt, Real threshold)
	{
		return checkContain(input_pnt) || checkNearSurface(input_pnt, threshold) ? true : false;
	}
	//=================================================================================================//
	bool ComplexShape::checkNearSurface(const Vec3d &input_pnt, Real threshold)
	{
		return getMaxAbsoluteElement(input_pnt - findClosestPoint(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	BoundingBox ComplexShape::findBounds()
	{
		//initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);

		for (size_t i = 0; i < triangle_mesh_shapes_.size(); i++)
		{
			BoundingBox shape_bounds = triangle_mesh_shapes_[i].first->findBounds();
			for (int j = 0; j != 3; ++j)
			{
				lower_bound[j] = SMIN(lower_bound[j], shape_bounds.first[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_bounds.second[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
}