#ifndef BASE_BODY_PART_HPP
#define BASE_BODY_PART_HPP

#include "base_body_part.h"

#include "base_particles.h"

namespace SPH
{
//=================================================================================================//
template <typename T, typename = void>
struct has_setupBaseParticles : std::false_type
{
};
//=================================================================================================//
template <typename T>
struct has_setupBaseParticles<T, std::void_t<decltype(&T::setupBaseParticles)>> : std::true_type
{
};
//=================================================================================================//
template <typename TargetCriterion>
class BodyPart::TargetParticleMask : public TargetCriterion
{
  public:
    template <class ExecutionPolicy, typename EnclosureType, typename... Args>
    TargetParticleMask(ExecutionPolicy &ex_policy, EnclosureType &encloser, Args &&...args)
        : TargetCriterion(std::forward<Args>(args)...), part_id_(encloser.part_id_),
          body_part_id_(encloser.dv_body_part_id_->DelegatedData(ex_policy)) {}
    ~TargetParticleMask() {}

    template <typename... Args>
    bool operator()(UnsignedInt target_index, Args &&...args)
    {
        return (body_part_id_[target_index] == part_id_) &&
               TargetCriterion::operator()(target_index, std::forward<Args>(args)...);
    }

  protected:
    int part_id_;
    int *body_part_id_;
};
//=================================================================================================//
template <typename TagCriteria>
BodyPartByParticle::BodyPartByParticle(SPHBody &sph_body, TagCriteria criteria)
    : BodyPartByParticle(sph_body)
{
    if constexpr (has_setupBaseParticles<TagCriteria>::value)
    {
        criteria.setupBaseParticles(base_particles_);
    }
    TaggingParticleMethod tagging_method = criteria;
    tagParticles(tagging_method);
}
//=================================================================================================//
template <typename DataType>
VariableRangeTagCriteria<DataType>::VariableRangeTagCriteria(
    const std::string &variable_name, DataType lower_bound, DataType upper_bound)
    : variable_name_(variable_name), lower_bound_(lower_bound),
      upper_bound_(upper_bound), variable_(nullptr)
{
    if (lower_bound_ > upper_bound_)
    {
        throw std::invalid_argument("Lower bound must be less than or equal to upper bound.");
    }
}
//=================================================================================================//
template <typename DataType>
void VariableRangeTagCriteria<DataType>::setupBaseParticles(BaseParticles &base_particles)
{
    variable_ = base_particles.template getVariableDataByName<DataType>(variable_name_);
}
//=================================================================================================//
template <typename DataType>
bool VariableRangeTagCriteria<DataType>::operator()(size_t index_i) const
{
    return (lower_bound_ <= variable_[index_i]) &&
           (variable_[index_i] <= upper_bound_);
}
//=================================================================================================//
template <typename... Args>
BodyPartByRealVar::BodyPartByRealVar(SPHBody &sph_body, Args &&...args)
    : BodyPartByParticle(sph_body, VariableRangeTagCriteria<Real>(std::forward<Args>(args)...)) {}
//=================================================================================================//
template <typename... Args>
BodyPartByIntVar::BodyPartByIntVar(SPHBody &sph_body, Args &&...args)
    : BodyPartByParticle(sph_body, VariableRangeTagCriteria<int>(std::forward<Args>(args)...)) {}
//=================================================================================================//
} // namespace SPH
#endif // BASE_BODY_PART_HPP