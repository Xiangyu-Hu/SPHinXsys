#include "complex_body_relation.h"
#include "base_particle_dynamics.h"
#include "contact_body_relation.h"
#include "inner_body_relation.h"

namespace SPH
{
//=================================================================================================//
ComplexRelation::
    ComplexRelation(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
    : SPHRelation(inner_relation.getSPHBody()),
      inner_relation_(inner_relation)
{
    contact_relations_.push_back(&contact_relation);
}
//=================================================================================================//
ComplexRelation::
    ComplexRelation(BaseInnerRelation &inner_relation, StdVec<BaseContactRelation *> contact_relations)
    : SPHRelation(inner_relation.getSPHBody()),
      inner_relation_(inner_relation)
{
    for (size_t k = 0; k != contact_relations.size(); ++k)
    {
        if (&inner_relation.getSPHBody() != &contact_relations[k]->getSPHBody())
        {
            std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        contact_relations_.push_back(contact_relations[k]);
    }
}
//=================================================================================================//
void ComplexRelation::updateConfiguration()
{
    inner_relation_.updateConfiguration();
    for (size_t k = 0; k != contact_relations_.size(); ++k)
        contact_relations_[k]->updateConfiguration();
}
//=================================================================================================//
} // namespace SPH
