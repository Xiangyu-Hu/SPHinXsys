/**
 * @file 	branch.h
 * @brief 	This is Branch class.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.2
 */
#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <algorithm>
#include <vector>
#include <iterator>

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <fstream>
#include <random>
using namespace std;

namespace SPH 
{
	class KDTree;
	class BaseLevelSet;
	
	using pointIndex = typename std::pair<Point, size_t>;

	class IntComparer 
	{
	public:
		IntComparer(){}
		~IntComparer(){}
		bool operator() (int i,int j) { return (i<j);}
	};
    /**
     * @class Node
     * @brief A class containing the nodes of the branches plus some fuctions to compute distance related quantities
     */
    class Node
	{
	protected:
			
	public:
		/** Default constructor. */
		Node(Point init_node);
        ~Node(){};
        /** Attributes */
		std::vector<Point> nodes_;	/**< list of points containing the coordinates of the nodes.*/
		std::vector<int> nodes_idx_; /**< lists of the node index. */
        int last_node_;         	/**< last added node. */
        IndexVector end_nodes_;     /**< a list containing the indices of all end nodes (nodes that are not connected). */
		KDTree* kd_tree_;			/**< a k-d tree to compute the distance from any point to the closest node in the tree, except from the brother and mother branches. 
										It is used to check collision between branches. */
		KDTree* collision_kd_tree_;
		IntComparer comparer_;			/**< Comparer for sorting. */

		/**
		 *@brief This function stores a list of nodes of a branch and returns the node indices. It also updates the tree to compute distances.
		 *@param[in] pts a list of arrays containing the coordinates of the nodes of one branch.
		 *@param[out] the indices of the added nodes.
		 */
		std::vector<int> addNode(std::vector<Point>& pts);
		/**
		 *@brief Ths function returns the distance from any point to the closest node in the tree.
		 *param[in] pt the point to calculate the distance from.
		 @param[out]   the distance between point and the closest node in the tree.
		 */
		Real getDistanceFromPoint(Point pt);
		/**
		 *@brief This function returns the distance from any node to the closest node in the tree.
		 *@param[in] node_idx the index of the node to calculate the distance from.
		 *@param[out]  the distance between specified node and the closest node in the tree.
		 */
		Real getDistanceFromNode(int node_idx);
		/**
		 *@brief This function updates the collision_tree excluding a list of nodes from all the nodes in the tree. 
		 * 		 If all the existing nodes are excluded, one distant node is added.
		 *@param[in] branch_nodes contains the nodes to exclude from the tree. Usually it should be the mother and the brother branch nodes.
		 */
		void updateCollisionTree(std::vector<int>& exclude_nodes);
		/**
		 *@brief This function returns the distance from any node to the closest node in the tree.
		 *@param[in] node_idx the index of the node to calculate the distance from.
		 *@param[out]  the distance between specified node and the closest node in the tree.
		 */
		Real checkCollitionOnKDTree(Point pnt);
		/**
		 *@brief  This function returns the distance between one point and the closest node in the tree and the index of 
		 *			the closest node using the collision_tree.
		 *@param[in] the point to calculate the distance from.
		 *@param[out] (distance to the closest node, index of the closest node)
		 */
		pointIndex getCollisionNode(Point pnt);
		/**
		 *@brief This function returns the gradient of the distance from the existing points of the tree from any point. 
		 *		 It uses a central finite difference approximation.
		 *@param[in] pt the point to calculate the gradient of the distance from.
		 *@param[out] gradient of the distance.
		 */
		Vecd getGradient(Point pt, Real delta);
    };
    /**
	 * @class BaseBranch
	 * @brief Class that contains a branch of the tree.  
	 */
	class Branch
	{
	protected:
		SPHBody* body_;
		BaseLevelSet* levelset_mesh_;
	public:
		Branch(SPHBody* body, 		/**< SPH Body */
				int parent_id,		/**< Parent branch ID. */
				int init_node_idx, 	/**< initial node to grow the branch. An index that refers to a node in the nodes.nodes. */
				Vecd init_dir, 		/**< initial direction to grow the branch. It refers to the direction of the last segment of the mother brach. */
				Real l, 			/**< length of the branch */
				Real angle, 		/**< angle (rad) with respect to the init_dir in the plane of surface */
				Real w, 			/**< repulsitivity parameter. Controls how much the branches repel each other */
				Node* node, 		/**< the object of the class nodes that contains all the nodes of the existing branches. */
				std::vector<int> &family_nodes_idx, /**< the nodes of the brother and mother branches, to be excluded from the collision detection. */
				int n_segments		/**< number of segments to divide the branch */
		);
		virtual ~Branch() {};

		std::vector<int> child_branchs_; 	/**< Contains the indexes of the child branches. It is not assigned when created*/
		int parent_branch_;					/**< Index of parent branch. */
		Vecd segment_dir_;					/**< Direction of the last segment of the branch.*/
		std::vector<int> node_idxs_;		/**< Contains the node indices of the branch. */
		std::vector<Point> branch_nodes_;	/**< Contains the node indices of the branch. */
		bool growing_ = true;				/**< False if the branch collide or is out of the surface. True otherwise.*/
		Real segment_length_;				/**< approximated length of segments. */
		/**
		 *@brief Functions that creats a new node in the mesh surface and it to the queue is it lies in the surface.
		 *@param[in] init_node vector that contains the coordinates of the last node added in the branch.
		 * 			 vector that contains the coordinates of the last node added in the branch.
		 *@param[in] dir a vector that contains the direction from the init_node to the node to project.
		 *@param[out] bool true if the new node is in the triangle.
		 */
		Point creatNewBranchNode(Point init_node, Vecd dir);
    };
	/**
	 * @class NetworkTree
	 * @brief Class that contains a NetworkTree
	 */
	class NetworkTree
	{
	protected:
		SPHBody* body_;
		BaseLevelSet* levelset_mesh_;;
		
		int n_it_ =15; 				/**< Number of iterations (generations of branches. */
		Real branch_angle_ = 0.15;	/**< angle with respect to the direction of the previous branch and the new branch. */
		Real branch_w_ = 0.2; 		/**< epulsivity parameter. */
		Real length_;				/**< edian length of the branches. */
		Real length_std_;			/**< standard deviation of the length. Set to zero to avoid random lengths */
		Real length_min_; 			/**< inimum length of the branches. To avoid randomly generated negative lengths*/
		bool fascicles_	= true;			/**< Create fascicles? */
		std::vector<Real> fascicle_angles_ = {-1.5, 0.2}; 	/**< angles with respect to the initial branches of the fascicles.*/
		std::vector<Real> fascicle_length_ = {0.5, 0.5}; 	/**< length  of the fascicles. Include one per fascicle to include.*/
	public:
		NetworkTree(SPHBody* body, 		/**< Pointer to SPH Body. */
					Point starting_pnt, /**< Starting point of net work. */
					Point second_pnt, 	/**< this point is only used to calculate the initial direction of the tree and is not included in the tree. */
					Real init_l, 		/**< Initial length of the first branch. */
					Real seg_l			/**< Length of the segments (approximately, because the lenght of the branch is random. */
					);
		virtual ~NetworkTree() {};
		Node* node_;						/**< Node container that contains all the nodes of the tree. */
		std::vector<Branch*> branches_;		/**< Branch container that contains all the nodes of the tree.*/
		int last_branch;					/**< index of last Branch. */
		std::vector<int> branch_to_grow_;	/**< index of Branch to grow. */
		/**
		 *@brief This function creats particles for SPH body based on the nodes growed in the network.
		 */
		void generateParticlesForSPHBody();
		/**
		 *@brief This function build inner configuration for SPH body in 1D format.
		 *		 The neighboring particles consists of 4 particles, all from same branch, two from parent branch and two from local branch,
		 *		 or two from local branch and another two from child branchs.
		 *		 Note, particles located at end of branch only have two neighboring particles.
		 */
		void buildInnerConfigurationgForSPHBody();
	};
}