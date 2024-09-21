#include "relax_dynamics_surface.h"
#include "level_set_shape.h"
#include "surface_shape.h"

#include <opencascade/GeomAPI_ProjectPointOnSurf.hxx>
#include <opencascade/gp_Pnt.hxx>

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
ShapeSurfaceConstraint::ShapeSurfaceConstraint(RealBody &real_body_, Shape &shape)
    : LocalDynamics(real_body_),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      shape_(&shape) {}
//=================================================================================================//
void ShapeSurfaceConstraint::update(size_t index_i, Real dt)
{
    pos_[index_i] = shape_->findClosestPoint(pos_[index_i]);
}
//=================================================================================================//
RelaxationStepInnerFirstHalf::
    RelaxationStepInnerFirstHalf(BaseInnerRelation &inner_relation)
    : BaseDynamics<void>(), real_body_(inner_relation.real_body_),
      inner_relation_(inner_relation), relaxation_acceleration_inner_(inner_relation) {}
//=================================================================================================//
void RelaxationStepInnerFirstHalf::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_.exec();
}

//=================================================================================================//
RelaxationStepInnerSecondHalf::
    RelaxationStepInnerSecondHalf(BaseInnerRelation &inner_relation, Shape &shape)
    : BaseDynamics<void>(), real_body_(inner_relation.real_body_),
      get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
      surface_bounding_(*real_body_, shape) {}
//=================================================================================================//
void RelaxationStepInnerSecondHalf::exec(Real dt)
{
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_bounding_.exec();
}

//=================================================================================================//
SurfaceNormalDirection::SurfaceNormalDirection(SPHBody &sph_body, Shape &shape)
    : LocalDynamics(sph_body),
      surface_shape_(DynamicCast<SurfaceShape>(this, &shape)),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      n_(particles_->registerStateVariable<Vecd>("NormalDirection")) {}

//=================================================================================================//
void SurfaceNormalDirection::update(size_t index_i, Real dt)
{
    Extrema_ExtAlgo Algo = Extrema_ExtAlgo_Tree;
    gp_Vec tangent_u, tangent_v;
    gp_Vec norm;
    Vecd normal_direction_;

    gp_Pnt point1 = EigenToOcct(pos_[index_i]);
    GeomAPI_ProjectPointOnSurf pointonsurf(point1, surface_shape_->surface_, Algo);

    Standard_Real u;
    Standard_Real v;
    gp_Pnt Point;

    pointonsurf.Parameters(1, u, v);

    surface_shape_->surface_->D1(u, v, Point, tangent_u, tangent_v);
    norm = tangent_u.Crossed(tangent_v);
    normal_direction_ = OcctVecToEigen(norm);
    n_[index_i] = normal_direction_.normalized();
}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
