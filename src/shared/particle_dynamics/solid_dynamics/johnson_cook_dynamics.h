#ifndef JOHNSON_COOK_DYNAMICS_H
#define JOHNSON_COOK_DYNAMICS_H

#include "elastic_dynamics.h"
#include "johnson_cook_solid.h"

namespace SPH
{
namespace solid_dynamics
{
/** Uniform total-Lagrangian integration; failed material transfers no traction. */
class JohnsonCookIntegration1stHalf : public Integration1stHalf
{
  public:
    explicit JohnsonCookIntegration1stHalf(BaseInnerRelation &inner_relation);
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    JohnsonCookSolid &johnson_cook_solid_;
    int *material_active_;
};

/** Keeps the failed deformation history fixed and retains ballistic mass. */
class JohnsonCookIntegration2ndHalf : public Integration2ndHalf
{
  public:
    explicit JohnsonCookIntegration2ndHalf(BaseInnerRelation &inner_relation);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *material_active_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // JOHNSON_COOK_DYNAMICS_H
