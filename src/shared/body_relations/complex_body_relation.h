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
 * @file 	complex_body_relation.h
 * @brief 	The topological relations within one body and to other bodies.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef COMPLEX_BODY_RELATION_H
#define COMPLEX_BODY_RELATION_H

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class ComplexRelation
 * @brief The relation combined an inner and a contact body relation.
 * The interaction is in a inner-boundary-condition fashion. Here inner interaction is
 * different from contact interaction.
 */
class ComplexRelation : public SPHRelation
{
  private:
    UniquePtrKeeper<BaseInnerRelation> base_inner_relation_ptr_keeper_;
    UniquePtrKeeper<BaseContactRelation> base_contact_relation_ptr_keeper_;

  protected:
    BaseInnerRelation &inner_relation_;
    BaseContactRelation &contact_relation_;

  public:
    BaseInnerRelation &getInnerRelation() { return inner_relation_; };
    BaseContactRelation &getContactRelation() { return contact_relation_; };
    RealBodyVector contact_bodies_;
    ParticleConfiguration &inner_configuration_;
    StdVec<ParticleConfiguration> &contact_configuration_;

    ComplexRelation(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation);
    ComplexRelation(RealBody &real_body, RealBodyVector contact_bodies);
    ComplexRelation(BaseInnerRelation &inner_relation, RealBodyVector contact_bodies);
    ComplexRelation(RealBody &real_body, BodyPartVector contact_body_parts);
    virtual ~ComplexRelation(){};

    virtual void resizeConfiguration() override;
    virtual void updateConfiguration() override;
};
} // namespace SPH
#endif // COMPLEX_BODY_RELATION_H