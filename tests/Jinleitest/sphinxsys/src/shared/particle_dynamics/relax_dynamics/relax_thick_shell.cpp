#include "relax_thick_shell.h"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
ShellMidSurfaceBounding::ShellMidSurfaceBounding(NearShapeSurface &body_part)
    : BaseLocalDynamics<BodyPartByCell>(body_part),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      constrained_distance_(0.5 * sph_body_.getSPHAdaptation().MinimumSpacing()),
      particle_spacing_ref_(sph_body_.getSPHAdaptation().MinimumSpacing()),
      level_set_shape_(DynamicCast<LevelSetShape>(this, &sph_body_.getInitialShape())) {}
//=================================================================================================//
void ShellMidSurfaceBounding::update(size_t index_i, Real dt)
{
    Vecd none_normalized_normal = level_set_shape_->findLevelSetGradient(pos_[index_i]);
    Vecd normal = level_set_shape_->findNormalDirection(pos_[index_i]);
    Real factor = 0.2 * none_normalized_normal.norm();
    pos_[index_i] -= factor * constrained_distance_ * normal;
}
//=================================================================================================//
ShellNormalDirectionPrediction::
    ShellNormalDirectionPrediction(BaseInnerRelation &inner_relation,
                                   Real thickness, Real consistency_criterion)
    : BaseDynamics<void>(),
      convergence_criterion_(cos(0.01 * Pi)),
      consistency_criterion_(consistency_criterion),
      normal_prediction_(inner_relation.getSPHBody(), thickness),
      normal_prediction_convergence_check_(inner_relation.getSPHBody(), convergence_criterion_),
      consistency_correction_(inner_relation, consistency_criterion_),
      consistency_updated_check_(inner_relation.getSPHBody()),
      smoothing_normal_(inner_relation) {}
//=================================================================================================//
void ShellNormalDirectionPrediction::exec(Real dt)
{
    predictNormalDirection();
    correctNormalDirection();
    predictNormalDirection();
    smoothing_normal_.exec();
}
//=================================================================================================//
void ShellNormalDirectionPrediction::predictNormalDirection()
{
    bool prediction_convergence = false;
    size_t ite_predict = 0;
    while (!prediction_convergence)
    {
        normal_prediction_.exec();
        prediction_convergence = normal_prediction_convergence_check_.exec();
        if (ite_predict > 100)
        {
            std::cout << "\n Error: class ShellNormalDirectionPrediction normal prediction not converged after 100 iterations." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        ite_predict++;
    }
    std::cout << "\n Information: normal direction prediction converged after '" << ite_predict << "' steps." << std::endl;
}
//=================================================================================================//
void ShellNormalDirectionPrediction::correctNormalDirection()
{
    bool consistency_updated = false;
    size_t ite_updated = 0;
    while (!consistency_updated)
    {
        consistency_correction_.exec();
        consistency_updated = consistency_updated_check_.exec();
        if (ite_updated > 100)
        {
            std::cout << "\n Error: class ShellNormalDirectionPrediction normal consistency not updated  after 100 iterations." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        ite_updated++;
    }
    std::cout << "\n Information: normal consistency updated after '" << ite_updated << "' steps." << std::endl;
}
//=================================================================================================//
ShellNormalDirectionPrediction::NormalPrediction::NormalPrediction(SPHBody &sph_body, Real thickness)
    : LocalDynamics(sph_body), thickness_(thickness),
      level_set_shape_(DynamicCast<LevelSetShape>(this, &sph_body.getInitialShape())),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      n_temp_(particles_->registerStateVariable<Vecd>(
          "PreviousNormalDirection", [&](size_t i) -> Vecd
          { return n_[i]; })) {}
//=================================================================================================//
void ShellNormalDirectionPrediction::NormalPrediction::update(size_t index_i, Real dt)
{
    n_temp_[index_i] = n_[index_i];
    n_[index_i] = level_set_shape_->findNormalDirection(pos_[index_i] + 0.3 * thickness_ * n_temp_[index_i]);
}
//=================================================================================================//
ShellNormalDirectionPrediction::PredictionConvergenceCheck::
    PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion)
    : LocalDynamicsReduce<ReduceAND>(sph_body),
      convergence_criterion_(convergence_criterion),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      n_temp_(particles_->getVariableDataByName<Vecd>("PreviousNormalDirection")) {}
//=================================================================================================//
bool ShellNormalDirectionPrediction::PredictionConvergenceCheck::reduce(size_t index_i, Real dt)
{
    return n_[index_i].dot(n_temp_[index_i]) > convergence_criterion_;
}
//=================================================================================================//
ShellNormalDirectionPrediction::ConsistencyCorrection::
    ConsistencyCorrection(BaseInnerRelation &inner_relation, Real consistency_criterion)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      consistency_criterion_(consistency_criterion),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      updated_indicator_(particles_->registerStateVariable<int>(
          "UpdatedIndicator", [&](size_t i) -> int
          { return 0; }))
{
    updated_indicator_[particles_->TotalRealParticles() / 3] = 1;
}
//=================================================================================================//
void ShellNormalDirectionPrediction::ConsistencyCorrection::interaction(size_t index_i, Real dt)
{
    mutex_modify_neighbor_.lock();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (updated_indicator_[index_i] == 1)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (updated_indicator_[index_j] == 0)
            {
                updated_indicator_[index_j] = 1;
                if (n_[index_i].dot(n_[index_j]) < consistency_criterion_)
                {
                    if (n_[index_i].dot(-n_[index_j]) < consistency_criterion_)
                    {
                        n_[index_j] = n_[index_i];
                        updated_indicator_[index_j] = 2;
                    }
                    else
                    {
                        n_[index_j] = -n_[index_j];
                        updated_indicator_[index_j] = 1;
                    }
                }
            }
        }
    }
    mutex_modify_neighbor_.unlock();
}
//=================================================================================================//
ShellNormalDirectionPrediction::ConsistencyUpdatedCheck::ConsistencyUpdatedCheck(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceAND>(sph_body),
      updated_indicator_(particles_->getVariableDataByName<int>("UpdatedIndicator")) {}
//=================================================================================================//
bool ShellNormalDirectionPrediction::ConsistencyUpdatedCheck::reduce(size_t index_i, Real dt)
{
    return updated_indicator_[index_i] != 0;
}
//=================================================================================================//
ShellNormalDirectionPrediction::SmoothingNormal::
    SmoothingNormal(BaseInnerRelation &inner_relation)
    : ParticleSmoothing<Vecd>(inner_relation, "NormalDirection") {};
//=================================================================================================//
void ShellNormalDirectionPrediction::SmoothingNormal::update(size_t index_i, Real dt)
{
    ParticleSmoothing<Vecd>::update(index_i, dt);
    smoothed_[index_i] /= temp_[index_i].norm() + TinyReal;
}
//=================================================================================================//
ShellRelaxationStep::ShellRelaxationStep(BaseInnerRelation &inner_relation)
    : BaseDynamics<void>(),
      real_body_(DynamicCast<RealBody>(this, inner_relation.getSPHBody())),
      inner_relation_(inner_relation), near_shape_surface_(real_body_),
      relaxation_residue_(inner_relation),
      relaxation_scaling_(real_body_), position_relaxation_(real_body_),
      mid_surface_bounding_(near_shape_surface_) {}
//=================================================================================================//
void ShellRelaxationStep::exec(Real ite_p)
{
    real_body_.updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_residue_.exec();
    Real scaling = relaxation_scaling_.exec();
    position_relaxation_.exec(scaling);
    mid_surface_bounding_.exec();
}
//=================================================================================================//
} // namespace relax_dynamics
} // namespace SPH
