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
 * @file surface_indication.h
 * @brief Here, we define the algorithm classes for indicating material surfaces.
 * @details Note that, SPHinXsys does not require surface indication
 * for simulating general free surface flow problems.
 * However, some other applications may use this function,
 * such as transport velocity formulation,
 * for masking some function which is only applicable for the bulk of the fluid body.
 * Currently, indicator used 0 for bulk, 1 for free surface indicator,
 * other to be defined.
 * @author	Xiangyu Hu
 */

#ifndef SURFACE_INDICATION_H
#define SURFACE_INDICATION_H

#include "base_general_dynamics.h"

namespace SPH
{
template <typename... InteractionTypes>
class FreeSurfaceIndication;

template <class DataDelegationType>
class FreeSurfaceIndication<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit FreeSurfaceIndication(BaseRelationType &base_relation);
    virtual ~FreeSurfaceIndication(){};

  protected:
    int *indicator_;
    Real *pos_div_, *Vol_;
    Real threshold_by_dimensions_;
};

template <>
class FreeSurfaceIndication<Inner<>>
    : public FreeSurfaceIndication<DataDelegateInner>
{
  public:
    explicit FreeSurfaceIndication(BaseInnerRelation &inner_relation);
    virtual ~FreeSurfaceIndication(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    bool isVeryNearFreeSurface(size_t index_i);
};

template <>
class FreeSurfaceIndication<Inner<SpatialTemporal>>
    : public FreeSurfaceIndication<Inner<>>
{
  public:
    explicit FreeSurfaceIndication(BaseInnerRelation &inner_relation);
    virtual ~FreeSurfaceIndication(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int *previous_surface_indicator_;
    bool isNearPreviousFreeSurface(size_t index_i);
};
using SpatialTemporalFreeSurfaceIndicationInner = FreeSurfaceIndication<Inner<SpatialTemporal>>;

template <>
class FreeSurfaceIndication<Contact<>>
    : public FreeSurfaceIndication<DataDelegateContact>
{
  public:
    explicit FreeSurfaceIndication(BaseContactRelation &contact_relation)
        : FreeSurfaceIndication<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(this->contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~FreeSurfaceIndication(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real *> contact_Vol_;
};

/**
 * @class NonWetting
 * @brief Non wetting surface particles include free-surface ones and interfacial ones near the non-wetted structure.
 * @brief Even the position divergence of interfacial fluid particles has satisfied with the threshold of spatial-temporal
   identification approach to be identified as internal ones,they will remain as free-surface ones if without
   any wetted neighboring solid particles.
 */
class NonWetting;

template <>
class FreeSurfaceIndication<Contact<NonWetting>>
    : public FreeSurfaceIndication<DataDelegateContact>
{
  public:
    FreeSurfaceIndication(BaseContactRelation &contact_relation);
    virtual ~FreeSurfaceIndication(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real *> contact_phi_, contact_Vol_;
};

using FreeSurfaceIndicationComplex =
    ComplexInteraction<FreeSurfaceIndication<Inner<>, Contact<>>>;

using SpatialTemporalFreeSurfaceIndicationComplex =
    ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>>>;

using WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex =
    ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<NonWetting>>>;
} // namespace SPH
#endif // SURFACE_INDICATION_H
