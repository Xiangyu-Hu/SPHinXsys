#include "body_part_by_cell_tracing.h"

#include "base_particles.hpp"
namespace SPH
{
    //=============================================================================================//
    //template <class TracingMethodType>
    //size_t BodyPartByCellWithTracing<TracingMethodType>::SizeOfLoopRange()
    //{
    //	size_t size_of_loop_range = 0;
    //    for (size_t i = 0; i != body_part_cells_.size(); ++i)
    //    {
    //        size_of_loop_range += body_part_cells_[i]->size();
    //    }
    //    return size_of_loop_range;
    //}
    ////=================================================================================================//
    //template <class TracingMethodType>
    //void BodyPartByCellWithTracing<TracingMethodType>::tagCells(TaggingCellMethod &tagging_cell_method, TracingCellMethod &tracing_cell_method)
    //{
    //    cell_linked_list_.tagBodyPartByCell(body_part_cells_, tagging_cell_method);
    //}
    //=================================================================================================//
    //template <class TracingMethodType>
    NearShapeSurfaceTracing::
        NearShapeSurfaceTracing(RealBody &real_body, SharedPtr<Shape> shape_ptr, BaseTracingMethod& tracing_cell_method_base)
        : BodyPartByCell(real_body, shape_ptr->getName()), //tracing_cell_method_(std::bind(&TracingMethodType::tracingPosition, this, _1, _2)),
        tagging_cell_method_(std::bind(&NearShapeSurfaceTracing::checkNearSurface, this, _1, _2)), tracing_cell_method_base_(tracing_cell_method_base),
          level_set_shape_(level_set_shape_keeper_.createRef<LevelSetShape>(real_body, *shape_ptr.get(), true))
    {
        //TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceTracing<TracingMethodType>::checkNearSurface, this, _1, _2);
        tagCells(tagging_cell_method_);
        
        
    }
    //=================================================================================================//
    //template <class TracingMethodType>
    //NearShapeSurfaceTracing::NearShapeSurfaceTracing(RealBody &real_body, BaseTracingMethod& tracing_cell_method_base)
    //    : BodyPartByCell(real_body, "NearShapeSurface"),//tracing_cell_method_(std::bind(&TracingMethodType::tracingPosition, this, _1, _2)),
    //    tagging_cell_method_(std::bind(&NearShapeSurfaceTracing::checkNearSurface, this, _1, _2)), tracing_cell_method_base_(tracing_cell_method_base),
    //      level_set_shape_(DynamicCast<LevelSetShape>(this, *real_body.body_shape_))
    //{
    //    //TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceTracing<TracingMethodType>::checkNearSurface, this, _1, _2);
    //    tagCells(tagging_cell_method_);
    //}
    //=================================================================================================//
    //template <class TracingMethodType>
    //NearShapeSurfaceTracing::NearShapeSurfaceTracing(RealBody &real_body, const std::string &shape_name, BaseTracingMethod& tracing_cell_method_base)
    //    : BodyPartByCell(real_body, shape_name), //tracing_cell_method_(std::bind(&TracingMethodType::tracingPosition, this, _1, _2)),
    //    tagging_cell_method_(std::bind(&NearShapeSurfaceTracing::checkNearSurface, this, _1, _2)),tracing_cell_method_base_(tracing_cell_method_base),
    //      level_set_shape_(DynamicCast<LevelSetShape>(this, *DynamicCast<ComplexShape>(this, real_body.body_shape_)->getShapeByName(shape_name)))
    //{
    //    //TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceTracing<TracingMethodType>::checkNearSurface, this, _1, _2);
    //    tagCells(tagging_cell_method_);
    //}
    //=================================================================================================//
    //template <class TracingMethodType>
    bool NearShapeSurfaceTracing::checkNearSurface(Vecd cell_position, Real threshold)
    {
        return level_set_shape_.checkNearSurface(tracing_cell_method_base_.tracingPosition(cell_position, 0.0), threshold);
    }
    //=================================================================================================//
    //template <class TracingMethodType>
    void NearShapeSurfaceTracing::updateCellList()
    {
        this->body_part_cells_.clear();
        tagCells(tagging_cell_method_);
    }
    //=================================================================================================//
} // namespace SPH
