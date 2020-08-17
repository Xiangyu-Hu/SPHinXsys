/**
* @file 	kd-tree.hpp
* @brief 	This is an adaptation of the KD-tree implementation
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_data_package.h"
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include <cmath>
#include <iterator>
#include <limits>

namespace SPH 
{
    using indexArr = std::vector<size_t>;
    using pointIndex = typename std::pair<Point, size_t>;
    using pointVec = std::vector<Point>;
    using pointIndexArr = typename std::vector<pointIndex>;

    /**
     *@class KDNode
     *@brief Node for KDTree
     */
    class KDNode 
    {
    public:
        /** std::shared_ptr is a smart pointer that retains shared ownership of an object through a pointer. 
         *  Several shared_ptr objects may own the same object. 
         */
        using KDNodePtr = std::shared_ptr<KDNode>;
        size_t index_;
        Point x_;
        KDNodePtr left_;
        KDNodePtr right_;
        bool is_empty_;
        /** Constructor. */
        KDNode(){is_empty_ = true;};
        KDNode(const Point &, const size_t &, const KDNodePtr &,const KDNodePtr &);
        KDNode(const pointIndex &, const KDNodePtr &, const KDNodePtr &);
        ~KDNode(){};

        Real coord(const size_t &);

        // conversions
        bool  is_empty();
        Point  position();
        size_t index();
        pointIndex posIndex();
    };
    /** Pointer to a empty KDNode. */
    using KDNodePtr = std::shared_ptr<KDNode>;
    KDNodePtr NewKDNodePtr();
    /** Get distance from two point or KD node. */
    inline double dist(const Point &, const Point &);
    inline double dist(const KDNodePtr &, const KDNodePtr &);
    /**
     *@class IndexComparer
     *@brief Compare for sorting.
     */
    class IndexComparer 
    {
    public:
        size_t idx_;
        explicit IndexComparer(size_t idx);
        inline bool compareIndex(const std::pair<Point, size_t> &, const std::pair<Point, size_t> & );
    };
    /**
     *@function sort_on_idx
     *@brief sorting on index.
     */
    inline void sortOnIndex(const pointIndexArr::iterator &, const pointIndexArr::iterator &, size_t idx);
    /**
    * @class KDTree
    * k-d tree implementation, based on the C version at rosettacode.org.
    */
    class KDTree 
    {
        KDNodePtr root;
        KDNodePtr leaf;

        KDNodePtr makeTree(const pointIndexArr::iterator &begin,const pointIndexArr::iterator &end,
                            const size_t &length, const size_t &level);

    public:
        KDTree() = default;
        explicit KDTree(pointVec point_array);

    private:
        KDNodePtr nearest_(const KDNodePtr &branch,const Point &pt, const size_t &level,const KDNodePtr &best, const Real &best_dist);
        KDNodePtr nearest_(const Point &pt);

    public:
        Point nearestPoint(const Point &pt);
        size_t nearestIndex(const Point &pt);
        pointIndex nearestPointIndex(const Point &pt);

    private:
        pointIndexArr neighborhood_(const KDNodePtr &branch, const Point &pt, const Real &rad,const size_t &level);
    public:
        pointIndexArr neighborhood(const Point &pt, const Real &rad);
        pointVec neighborhoodPoints(const Point &pt, const Real &rad);
        indexArr neighborhoodIndices(const Point &pt,const Real &rad);
    };
}