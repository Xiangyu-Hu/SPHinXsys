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
 * @brief The relation combined an inner and one or several contact body relation.
 * Note that this relation are not used for construct local dynamics,
 * which are only done by using inner and contact relations.
 * This is temporary class used for updating several configuration together,
 * before the automation on updating is done.
 */
class ComplexRelation : public SPHRelation
{
  protected:
    BaseInnerRelation &inner_relation_;
    StdVec<BaseContactRelation *> contact_relations_;

  public:
    ComplexRelation(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation);
    ComplexRelation(BaseInnerRelation &inner_relation, StdVec<BaseContactRelation *> contact_relations);
    virtual ~ComplexRelation(){};

    virtual void updateConfiguration() override;
};
} // namespace SPH
#endif // COMPLEX_BODY_RELATION_H