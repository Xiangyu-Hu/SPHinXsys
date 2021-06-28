/**
 * @file 	base_geometry.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "base_geometry.h"

namespace SPH {
	//=================================================================================================//
	Branch::Branch(Vecd init_point, Vecd auxillary_point, Tree* tree) :
		Edge<size_t, IndexVector>(),
		is_end_(false)
	{
		in_edge_ = 0;
		id_ = 0;
		inner_points_.push_back(tree->points_.size() - 1);
		Vecd displacement = auxillary_point - init_point;
		end_direction_ = displacement / (displacement.norm() + TinyReal);
	}
	//=================================================================================================//
	Branch::Branch(size_t parent_id, Tree* tree) : 
		Edge<size_t, IndexVector>(parent_id),
		is_end_(false)
	{
		id_ = tree->branches_.size();
		tree->branches_[parent_id]->out_edge_.push_back(id_);
	}
	//=================================================================================================//
	size_t Structure::EdgeLocation(size_t point_index)
	{
		return point_index < points_.size() ? edge_locations_[point_index] : MaxSize_t;
	}
	//=================================================================================================//
	Tree::Tree(Vecd init_point, Vecd auxillary_point) : Structure(), last_branch_id_(0)
	{
		points_.push_back(init_point);
		edge_locations_.push_back(0);
		face_locations_.push_back(0);

		Branch* first_branch = new Branch(init_point, auxillary_point, this);
		branches_.push_back(first_branch);
	}
	//=================================================================================================//
	Tree::~Tree() 
	{
		for (size_t i = 0; i != branches_.size(); i++) delete branches_[i];
	}
	//=================================================================================================//
	void Tree::addANewBranch(Branch* branch)
	{
		branches_.push_back(branch);
		last_branch_id_ = branch->id_;
	}
	//=================================================================================================//
	void Tree::addANewBranchInnerPoints(Branch* branch, Vecd new_point, Vecd end_direction)
	{
		points_.push_back(new_point);
		edge_locations_.push_back(branch->id_);
		branch->inner_points_.push_back(points_.size() - 1);
		branch->end_direction_ = end_direction;
	}
	//=================================================================================================//
}
