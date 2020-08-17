/**
 * @file 	branch.cpp
 * @brief 	Here, Functions belong to BaseBranch are given.
 * @author	chi ZHang and Xiangyu Hu
 * @version	0.2.0
 */
 #include "kd_tree.hpp"
 #include "network.h"
 //=================================================================================================//
namespace SPH
{
    //=================================================================================================//
    Node::Node(Point init_node)
	{
		nodes_.push_back(init_node);
		last_node_ = 0;
        nodes_idx_.push_back(last_node_);
		kd_tree_ = new KDTree(nodes_);
	}
    //=================================================================================================//
    std::vector<int> Node::addNode(std::vector<Point>& pts)
    {
        std::vector<int> node_id;
        for (const Point& pnt : pts)
        {
            nodes_.push_back(pnt);
            last_node_ ++;
            nodes_idx_.push_back(last_node_);
            node_id.push_back(last_node_);
        }

        kd_tree_ = new KDTree(nodes_);
        return node_id;
    }
    //=================================================================================================//
    Real Node::getDistanceFromPoint(Point pt)
    {
        Point nearestpnt = kd_tree_->nearestPoint(pt);
        return (nearestpnt - pt).norm();
    }
    //=================================================================================================//
    Real Node::getDistanceFromNode(int node_idx)
    {
        Point pnt = nodes_[node_idx];
        Point nearestpnt = kd_tree_->nearestPoint(pnt);
        return (nearestpnt - pnt).norm();
    }
    //=================================================================================================//
    Real Node::checkCollitionOnKDTree(Point pnt)
    {
        Point nearestpnt = collision_kd_tree_->nearestPoint(pnt);
        return (nearestpnt - pnt).norm();
    }
    //=================================================================================================//
    void Node::updateCollisionTree(std::vector<int>& exclude_nodes)
    {   
        std::vector<Point> collision_check;
        std::vector<int> idx_to_consider;

        idx_to_consider = differenceVector(nodes_idx_, exclude_nodes);

        for(const int& i : idx_to_consider)
        {
            collision_check.push_back(nodes_[i]);
        }
        /** Update collision tree. */
        collision_kd_tree_ = new KDTree(collision_check);
    } 
    //=================================================================================================//
    pointIndex Node::getCollisionNode(Point pnt)
    {
        std::pair<Point, size_t> point_idx = collision_kd_tree_->nearestPointIndex(pnt);
        return point_idx;
    }
    //=================================================================================================//
	Branch::Branch(SPHBody* body,int parent_id, int init_node_idx, Vecd init_dir, Real l, Real angle, Real w, 
                    Node* node, std::vector<int> &family_nodes_idx,  int n_segments)
	{
        body_= body;
        levelset_mesh_ = body_->levelset_mesh_;
        parent_branch_ = parent_id;
        node_idxs_.push_back(init_node_idx);
        branch_nodes_.push_back(node->nodes_[init_node_idx]);
        segment_length_ = l / (Real(n_segments));
        node->updateCollisionTree(family_nodes_idx);
        
        Vecd init_node_pos = branch_nodes_[0];
        Vecd surface_norm = levelset_mesh_->probeNormalDirection(init_node_pos);
        surface_norm /= surface_norm.norm() + TinyReal;
        Vecd in_plane = -getCrossProduct(init_dir, surface_norm);

        Vecd dir = cos(angle) * init_dir + sin(angle) * in_plane;
        dir /= dir.norm() + TinyReal;
        Vecd grad = node->getGradient(branch_nodes_[0], 0.01);
        segment_dir_ = (w * grad + dir) / ((w * grad + dir).norm() + TinyReal);

        for(int i = 1; i < n_segments; i++)
        {
            Point new_node = creatNewBranchNode(branch_nodes_[i-1], segment_dir_);
            Real dist_to_collision_node =  node->checkCollitionOnKDTree(new_node);
            if(dist_to_collision_node < 2.0 * segment_length_)
            {
                growing_ = false;
                std::cout<< "Branch Collision Detected, Break! " << std::endl;
                break;
            }
            branch_nodes_.push_back(new_node);

            surface_norm = levelset_mesh_->probeNormalDirection(new_node);
            surface_norm /= surface_norm.norm() + TinyReal;
            /** Project grad to surface. */
            grad = node->getGradient(new_node, 0.01);
            grad -= dot(grad, surface_norm) * surface_norm;
            dir = (w * grad + segment_dir_) / ((w * grad + segment_dir_).norm() + TinyReal);
            segment_dir_= dir;
        }
        std::vector<Point> new_created_branch_nodes =  subvector(branch_nodes_, 1, branch_nodes_.size()-1);
        std::vector<int> new_node_idxs = node->addNode(new_created_branch_nodes);
        node_idxs_.insert(node_idxs_.end(), new_node_idxs.begin(), new_node_idxs.end());
    }
    //=================================================================================================//
    Point Branch::creatNewBranchNode(Point init_node, Vecd dir)
    {
        Vecd pnt_to_project = init_node + dir * segment_length_;

		Real phi = levelset_mesh_->probeLevelSet(pnt_to_project);
		Vecd unit_normal = levelset_mesh_->probeNormalDirection(pnt_to_project);
		unit_normal /= unit_normal.norm() + TinyReal;
		
		return pnt_to_project - phi * unit_normal;
    }
    //=================================================================================================//
	NetworkTree::NetworkTree(SPHBody* body, Point starting_pnt, Point second_pnt, Real init_l, Real seg_l)
	{
        std::vector<int> family_node_idxs;
        body_= body;
        levelset_mesh_ = body_->levelset_mesh_;

        node_ = new Node(starting_pnt);
        family_node_idxs.push_back(0);
        Real length_ = init_l;			
		Real length_std_ = sqrt(0.2) * length_;	
		Real length_min_ = 0.1 * length_; 
        
        Vecd init_dir = (second_pnt - starting_pnt) / ( (second_pnt - starting_pnt).norm() + TinyReal);
        int n_seg = int(init_l / seg_l);
        int last_branch = 0;
        branches_.push_back(new Branch(body_, last_branch, 0, init_dir, init_l, 0.0, 0.0, node_, family_node_idxs, n_seg));
        std::cout<< "Branch No :" << last_branch << std::endl;

        std::vector<int> branches_to_grow;
        std::vector<int> new_branches_to_grow;
        Real angle_to_use;
        Real l_to_use;
        branches_to_grow.push_back(last_branch);

        if(fascicles_)
        {
            branches_to_grow.clear();
            family_node_idxs.clear();
            family_node_idxs.insert(family_node_idxs.end(),branches_[0]->node_idxs_.begin(), branches_[0]->node_idxs_.end());

            for(int i = 0; i < 2; i++)
            {
                angle_to_use = fascicle_angles_[i];
                l_to_use = fascicle_length_[i];
                n_seg = int(l_to_use / seg_l);
                branches_.push_back(new Branch(body_, 0, branches_[0]->node_idxs_.back(), branches_[0]->segment_dir_, l_to_use, 
                                        angle_to_use, 0.0, node_, family_node_idxs, n_seg));
                last_branch ++;
                family_node_idxs.insert(family_node_idxs.end(),branches_[last_branch]->node_idxs_.begin(), branches_[last_branch]->node_idxs_.end());
                branches_to_grow.push_back(last_branch);
                branches_[0]->child_branchs_.push_back(last_branch);
            }       
        }

        for(int i = 0; i < n_it_; i++)
        {
            new_branches_to_grow.clear();
            shuffle(branches_to_grow);
            for(int j = 0; j < branches_to_grow.size(); j++)
            {
                int grow_id = branches_to_grow[j];
                angle_to_use = -2.0 * branch_angle_ * (((double)rand() / (RAND_MAX)) - 0.5);
                for(int k = 0; k < 2; k++)
                {
                    family_node_idxs.clear();
                    family_node_idxs.insert(family_node_idxs.end(),branches_[grow_id]->node_idxs_.begin(), branches_[grow_id]->node_idxs_.end());
                    if(k > 0)
                        family_node_idxs.insert(family_node_idxs.end(), branches_[last_branch]->node_idxs_.begin(), branches_[last_branch]->node_idxs_.end());
                    /** Add new branch. */
                    Real l = length_ + rand_norm(0.0, length_std_);
                    l_to_use = l < length_min_ ? length_min_ : l;
                    n_seg = int(l_to_use / seg_l);
                    branches_.push_back(new Branch(body_,grow_id, branches_[grow_id]->node_idxs_.back(), branches_[grow_id]->segment_dir_, l_to_use, 
                                        angle_to_use, branch_w_, node_, family_node_idxs, n_seg));
                    last_branch ++;
                    branches_[grow_id]->child_branchs_.push_back(last_branch);
                    std::cout<< "Branch No :" << last_branch << std::endl;
                    if(branches_[last_branch]->growing_)
                    {
                        new_branches_to_grow.push_back(last_branch);
                    }

                    angle_to_use *= -1.0;
                }
            }
            branches_to_grow.clear();
            branches_to_grow = new_branches_to_grow;
        }

        generateParticlesForSPHBody();
    }
    //=================================================================================================//
    void NetworkTree::generateParticlesForSPHBody()
    {
        std::cout << "Now converting network node to particles..." << "\n" << std::endl;
        for(int i = 0; i< node_->nodes_.size(); i++ )
        {
            body_->body_input_points_volumes_.push_back(make_pair(node_->nodes_[i], 0.0));
        }
        std::cout << body_->body_input_points_volumes_.size() << "particles are added to body" << "\n" << std::endl;
    }
    //=================================================================================================//
    void NetworkTree::buildInnerConfigurationgForSPHBody()
    {
        // parallel_for(blocked_range<size_t>(1, branches_.size()),
		// 	[&](const blocked_range<size_t>& r) 
        //     {
		// 		for (size_t num = r.begin(); num != r.end(); ++num) 
        //         {
        //             for(int i = 1; i < branches_[num]->node_idxs_.size(); i++)
        //             {   std::vector<int> neighbor_nodes_;
        //                 if(i == 1)
        //                 {
        //                     int parent_id = branches_[num]->parent_branch_;
        //                     neighbor_nodes_.push_back(branches_[parent_id]->node_idxs_.end());
        //                     neighbor_nodes_.push_back(branches_[parent_id]->node_idxs_.end()[-2]);
        //                     for(int j = 1; j != 3; j++)
        //                     {
        //                         neighbor_nodes_.push_back(branches_[num]->node_idxs_[i + j]);
        //                     }
        //                 }

        //                 if(i == 2)
        //                 {
        //                     int parent_id = branches_[num]->parent_branch_;
        //                     neighbor_nodes_.push_back(branches_[parent_id]->node_idxs_.end());
        //                     neighbor_nodes_.push_back(branches_[parent_id]->node_idxs_.end()[-2]);
        //                     for(int j = 1; j != 3; j++)
        //                     {
        //                         neighbor_nodes_.push_back(branches_[num]->node_idxs_[i + j]);
        //                     }
        //                 }

        //                 int particle_id = branches_[num]->node_idxs_[i];
        //                 Neighborhood& neighborhood = body_->inner_configuration[particle_id];
		// 			    NeighborList& neighbor_list = std::get<0>(neighborhood);
		// 			    size_t previous_count_of_neigbors = std::get<2>(neighborhood);

        //             }
        //             int particle_id = branches_[num]->node_idxs_[];

		// 			Neighborhood& neighborhood = inner_configuration[num];
		// 			NeighborList& neighbor_list = std::get<0>(neighborhood);
		// 			size_t previous_count_of_neigbors = std::get<2>(neighborhood);

		// 			for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
		// 				for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
		// 				{
		// 					ConcurrentListDataVector& target_particles = cell_linked_lists_[l][m].particle_data_lists_;
		// 					for (size_t n = 0; n != target_particles.size(); ++n)
		// 					{
		// 						//displacement pointing from neighboring particle to origin particle
		// 						Vecd displacement = base_particle_data[num].pos_n_ - target_particles[n].second;
		// 						if (displacement.norm() <= cutoff_radius_ && num != target_particles[n].first)
		// 						{
		// 							std::get<1>(neighborhood) >= neighbor_list.size() ?
		// 								neighbor_list.emplace_back(new NeighborRelation(base_particle_data, *kernel_,
		// 									displacement, num, target_particles[n].first))
		// 								: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data,
		// 									*kernel_, displacement, num, target_particles[n].first);
		// 							std::get<1>(neighborhood)++;
		// 						}
		// 					}
		// 				}
		// 			std::get<2>(neighborhood) = std::get<1>(neighborhood);
		// 			std::get<1>(neighborhood) = 0;
		// 	}
		// }, ap);
    }
    //=================================================================================================//
}
//=================================================================================================//