/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	body_part_by_cell_tracing.h
 * @brief 	This is the base classes of body parts by cell with tracing method.
 * @details	By default, the tracing method is empty, but with real tracing method, 
            it can provide the position of each cell under the given imaginary motion 
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef BODY_PART_BY_CELL_TRACING_H
#define BODY_PART_BY_CELL_TRACING_H

#include "base_body.h"
#include "base_body_part.h"
#include <string>

namespace SPH
{
/**
 * @class BodyPart
 * @brief An auxillary class for SPHBody to indicate a part of the body.
 */
using namespace std::placeholders;

class BodyPart;

class BaseTracingMethod
{
 public:
     BaseTracingMethod() {};
     virtual ~BaseTracingMethod(){};

     virtual Vecd tracingPosition (Vecd previous_position, Real current_time = 0.0)
     {
         Vecd current_position = previous_position;
         return current_position;
     }

     virtual Vecd updateNormalForVector(Vecd previous_vector)
     {
         return previous_vector;
     }
};

/**
 * @class BodyPartByCell
 * @brief A body part with a collection of cell lists.
 */
//template <class TracingMethodType = BaseTracingMethod>
//class BodyPartByCellWithTracing : public BodyPart, public TracingMethodType
//{
//  public:
//    ConcurrentCellLists body_part_cells_; /**< Collection of cells to indicate the body part. */
//    ConcurrentCellLists &LoopRange() { return body_part_cells_; };
//    size_t SizeOfLoopRange();
//
//    BodyPartByCellWithTracing(RealBody& real_body, const std::string& body_part_name)
//        : BodyPart(real_body, body_part_name), cell_linked_list_(real_body.getCellLinkedList())
//    {
//        tracing_cell_method_ = std::bind(&TracingMethodType::tracingPosition, this, _1, _2);
//    };
//    virtual ~BodyPartByCellWithTracing(){};
//
//  protected:
//    BaseCellLinkedList &cell_linked_list_;
//    typedef std::function<bool(Vecd, Real)> TaggingCellMethod;
//    typedef std::function<Vecd(Vecd, Real)> TracingCellMethod;
//    TracingCellMethod tracing_cell_method_;
//    void tagCells(TaggingCellMethod &tagging_cell_method, TracingCellMethod &tracing_cell_method);
//};


/**
 * @class NearShapeSurface
 * @brief A body part with the cell lists near the surface of a prescribed shape.
 */
//template <class TracingMethodType= BaseTracingMethod>
class NearShapeSurfaceTracing : public BodyPartByCell
{
private:
    UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;
public:
    LevelSetShape& level_set_shape_;
    BaseTracingMethod& tracing_cell_method_base_;

    NearShapeSurfaceTracing(RealBody& real_body, SharedPtr<Shape> shape_ptr, BaseTracingMethod& tracing_cell_method_base);
    //explicit NearShapeSurfaceTracing(RealBody &real_body, BaseTracingMethod& tracing_cell_method_base);
    //NearShapeSurfaceTracing(RealBody &real_body, const std::string &shape_name, BaseTracingMethod& tracing_cell_method_base);
    virtual ~NearShapeSurfaceTracing(){};
    void updateCellList();
private:
    bool checkNearSurface(Vecd cell_position, Real threshold);
    typedef std::function<Vecd(Vecd, Real)> TracingCellMethod;
    //TracingCellMethod tracing_cell_method_;
    TaggingCellMethod tagging_cell_method_;
    
};

} // namespace SPH
#endif // BODY_PART_BY_CELL_TRACING_H
