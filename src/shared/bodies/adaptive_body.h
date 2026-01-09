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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    adaptive_body.h
 * @brief 	This is the class for templated bodies.
 * @author	Xiangyu Hu
 */

#ifndef ADAPTIVE_BODY_H
#define ADAPTIVE_BODY_H

#include "base_body.h"

namespace SPH
{
template <typename...>
class AdaptiveBody;

template <class AdaptationType, class BaseBodyType>
class AdaptiveBody<AdaptationType, BaseBodyType> : public BaseBodyType
{
    AdaptationType adaptation_;

  public:
    typedef AdaptationType Adaptation;
    using BaseAdaptation = typename AdaptationType::BaseAdaptation;
    using SpacingAdaptation = typename AdaptationType::SpacingAdaptation;

    template <typename... Args>
    AdaptiveBody(SPHSystem &sph_system, AdaptationType adaptation, Args &&...args)
        : BaseBodyType(sph_system, std::forward<Args>(args)...), adaptation_(adaptation)
    {
        this->sph_adaptation_ = &adaptation_;
    };

    virtual ~AdaptiveBody() {};
    AdaptationType &getAdaptation() { return adaptation_; };

    virtual void createCellLinkedListPtr() override
    {
        this->cell_linked_list_ptr_ = adaptation_.createFinestCellLinkedList(
            this->getSPHSystemBounds(), *this->base_particles_);
    };
};
} // namespace SPH
#endif // ADAPTIVE_BODY_H