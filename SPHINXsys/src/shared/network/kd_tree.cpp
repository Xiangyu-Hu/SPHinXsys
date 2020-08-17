/*
 *@file KDTree.cpp
 *@brief This is an adaptation of the KD-tree implementation in rosetta code
 *        https://rosettacode.org/wiki/K-d_tree
 *@author Chi ZHANG 
 *@version 0.1
 */
#include "kd_tree.hpp"
//=============================================================================================//
namespace SPH 
{
    //=============================================================================================//
    KDNode::KDNode(const Point &pt, const size_t &idx, const KDNodePtr &left, const KDNodePtr &right) 
    {
        is_empty_ = false;
        x_ = pt;
        index_ = idx;
        left_ = left;
        right_ = right;
    }
    //=============================================================================================//
    KDNode::KDNode(const pointIndex &pi, const KDNodePtr &left, const KDNodePtr &right) 
    {
        is_empty_ = false;
        x_ = pi.first;
        index_ = pi.second;
        left_ = left;
        right_ = right;
    }
    //=============================================================================================//
    Real KDNode::coord(const size_t &idx) 
    {
        return x_.get(idx); 
    }
    //=============================================================================================//
    bool KDNode::is_empty() 
    {
        return is_empty_;
    }
    //=============================================================================================//
    Point KDNode::position() 
    { 
        return x_; 
    }
    //=============================================================================================//
    size_t KDNode::index() 
    { 
        return index_; 
    }
    //=============================================================================================//
    pointIndex KDNode::posIndex() 
    { 
        return pointIndex(x_, index_); 
    }
    //=============================================================================================//
    KDNodePtr NewKDNodePtr() 
    {
        KDNodePtr mynode = std::make_shared<KDNode>();
        return mynode;
    }
    //=============================================================================================//
    inline Real dist2(const Point &a, const Point &b) 
    {
        Real distc = 0;
        for (size_t i = 0; i < a.size(); i++) 
        {
            Real di = a[i] - b[i];
            distc += di * di;
        }
        return distc;
    }
    //=============================================================================================//
    inline Real dist2(const KDNodePtr &a, const KDNodePtr &b) 
    {
        return dist2(a->x_, b->x_);
    }
    //=============================================================================================//
    IndexComparer::IndexComparer(size_t idx) : idx_{idx} {};
    //=============================================================================================//
    inline bool IndexComparer::compareIndex(const pointIndex &a,const pointIndex &b ) 
    {
        return (a.first[idx_] < b.first[idx_]); 
    }
    //=============================================================================================//
    inline void sortOnIndex(const pointIndexArr::iterator &begin, const pointIndexArr::iterator &end,size_t idx) 
    {
        IndexComparer comp(idx);
        comp.idx_ = idx;

        using std::placeholders::_1;
        using std::placeholders::_2;

        std::sort(begin, end, std::bind(&IndexComparer::compareIndex, comp, _1, _2));
    }

    using pointVec = std::vector<Point>;
    //=============================================================================================//
    KDNodePtr KDTree::makeTree(const pointIndexArr::iterator &begin, const pointIndexArr::iterator &end,    
        const size_t &length,const size_t &level) 
    {
        if (begin == end) 
        {
            return NewKDNodePtr();
        }

        size_t dim = begin->first.size();

        if (length > 1) 
        {
            sortOnIndex(begin, end, level);
        }

        auto middle = begin + (length / 2);

        auto l_begin = begin;
        auto l_end = middle;
        auto r_begin = middle + 1;
        auto r_end = end;

        size_t l_len = length / 2;
        size_t r_len = length - l_len - 1;

        KDNodePtr left;
        if (l_len > 0 && dim > 0) {
            left = makeTree(l_begin, l_end, l_len, (level + 1) % dim);
        } else {
            left = leaf;
        }
        KDNodePtr right;
        if (r_len > 0 && dim > 0) {
            right = makeTree(r_begin, r_end, r_len, (level + 1) % dim);
        } else {
            right = leaf;
        }

        KDNodePtr root_ptr =  std::make_shared<KDNode>(*middle, left, right);
        return root_ptr;
    }
    //=============================================================================================//
    KDTree::KDTree(pointVec point_array) 
    {
        leaf = std::make_shared<KDNode>();
        // iterators
        pointIndexArr arr;
        for (size_t i = 0; i < point_array.size(); i++) 
        {
            arr.push_back(pointIndex(point_array.at(i), i));
        }

        auto begin = arr.begin();
        auto end = arr.end();

        size_t length = arr.size();
        size_t level = 0;
        root = KDTree::makeTree(begin, end, length, level);
    }
    //=============================================================================================//
    KDNodePtr KDTree::nearest_(const KDNodePtr &branch, const Point &pt, const size_t &level,
        const KDNodePtr &best, const Real &best_dist) 
    {
        Real d, dx, dx2;

        if ((*branch).is_empty()) 
        {
            return NewKDNodePtr();
        }

        Point branch_pt = (*branch).position();
        size_t dim = branch_pt.size();
        d = dist2(branch_pt, pt);
        dx = branch_pt.get(level) - pt.get(level);
        dx2 = dx * dx;

        KDNodePtr best_l = best;

        Real best_dist_l = best_dist;

        if (d < best_dist) {
            best_dist_l = d;
            best_l = branch;
        }

        size_t next_lv = (level + 1) % dim;
        KDNodePtr section;
        KDNodePtr other;

        // select which branch makes sense to check
        if (dx > 0) {
            section = branch->left_;
            other = branch->right_;
        } else {
            section = branch->right_;
            other = branch->left_;
        }

        // keep nearest neighbor from further down the tree
        KDNodePtr further = nearest_(section, pt, next_lv, best_l, best_dist_l);
        if (!(*further).is_empty())
        {
            Real dl = dist2(further->x_, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
            }
        }
        // only check the other branch if it makes sense to do so
        if (dx2 < best_dist_l) 
        {
            further = nearest_(other, pt, next_lv, best_l, best_dist_l);
            if (!(*further).is_empty())
            {
                Real dl = dist2(further->x_, pt);
                if (dl < best_dist_l) {
                    best_dist_l = dl;
                    best_l = further;
                }
            }
        }

        return best_l;
    };
    //=============================================================================================//
    KDNodePtr KDTree::nearest_(const Point &pt)
    {
        size_t level = 0;
        // KDNodePtr best = branch;
        Real branch_dist = dist2((*root).position(), pt);

        return nearest_(root,          // beginning of tree
                        pt,            // point we are querying
                        level,         // start from level 0
                        root,          // best is the root
                        branch_dist);  // best_dist = branch_dist
    };
    //=============================================================================================//
    Point KDTree::nearestPoint(const Point &pt) 
    {
        KDNodePtr nearest_kd = nearest_(pt);
        return ((*nearest_kd).position());
    }
    //=============================================================================================//
    size_t KDTree::nearestIndex(const Point &pt) 
    {
        KDNodePtr nearest_kd = nearest_(pt);
        return ((*nearest_kd).index());
    }
    //=============================================================================================//
    pointIndex KDTree::nearestPointIndex(const Point &pt) 
    {
        KDNodePtr Nearest = nearest_(pt);
        return pointIndex((*Nearest).position(), (*Nearest).index());
    }
    //=============================================================================================//
    pointIndexArr KDTree::neighborhood_(const KDNodePtr &branch, const Point &pt, 
        const Real &rad, const size_t &level) 
    {
        Real d, dx, dx2;

        if ((*branch).is_empty()) 
        {
            return pointIndexArr();
        }

        size_t dim = pt.size();

        Real r2 = rad * rad;

        d = dist2((*branch).position(), pt);
        dx = (*branch).position().get(level) - pt.get(level);
        dx2 = dx * dx;

        pointIndexArr nbh, nbh_s, nbh_o;
        if (d <= r2) {
            nbh.push_back((*branch).posIndex());
        }

        KDNodePtr section;
        KDNodePtr other;
        if (dx > 0) {
            section = branch->left_;
            other = branch->right_;
        } else {
            section = branch->right_;
            other = branch->left_;
        }

        nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
        if (dx2 < r2) {
            nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
            nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
        }

        return nbh;
    };
    //=============================================================================================//
    pointIndexArr KDTree::neighborhood(const Point &pt, const Real &rad) 
    {
        size_t level = 0;
        return neighborhood_(root, pt, rad, level);
    }
    //=============================================================================================//
    pointVec KDTree::neighborhoodPoints(const Point &pt,const Real &rad) 
    {
        size_t level = 0;
        pointIndexArr nbh = neighborhood_(root, pt, rad, level);
        pointVec nbhp;
        nbhp.resize(nbh.size());
        std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                    [](pointIndex x) { return x.first; });
        return nbhp;
    }
    //=============================================================================================//
    indexArr KDTree::neighborhoodIndices(const Point &pt, const Real &rad) 
    {
        size_t level = 0;
        pointIndexArr nbh = neighborhood_(root, pt, rad, level);
        indexArr nbhi;
        nbhi.resize(nbh.size());
        std::transform(nbh.begin(), nbh.end(), nbhi.begin(),
                    [](pointIndex x) { return x.second; });
        return nbhi;
    }
    //=============================================================================================//
}
//=============================================================================================//