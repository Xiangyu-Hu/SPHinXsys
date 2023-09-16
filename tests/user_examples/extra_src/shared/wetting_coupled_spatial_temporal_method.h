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
 * @file 	fluid_surface_complex.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef WETTING_COUPLED_SPATIAL_TEMPORAL_COMPLEX_H
#define WETTING_COUPLED_SPATIAL_TEMPORAL_COMPLEX_H

#include "fluid_surface_complex.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class NonWettingSurfaceIndication
 * @brief Non wetting surface particles include free-surface ones and interfacial ones near the non-wetted structure.
 * @brief Even the position divergence of interfacial fluid pariticles has satisfied with the threshold of spatial-temporal 
   identification approach to be identified as internal ones,they will remain as free-surface ones if without 
   any wetted neighboring solid particles.
 */
class NonWettingSurfaceIndication : public FreeSurfaceIndicationComplex
{
  public:
    NonWettingSurfaceIndication(BaseInnerRelation &inner_relation,
                                               BaseContactRelation &contact_relation, Real threshold = 0.75, Real criterion = 0.99);
    explicit NonWettingSurfaceIndication(ComplexRelation &complex_relation, Real threshold = 0.75, Real criterion = 0.99);
    virtual ~NonWettingSurfaceIndication(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        FreeSurfaceIndicationInner::interaction(index_i, dt);

        Real pos_div = 0.0;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                pos_div -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.r_ij_[n];
            }
        }
        pos_div_[index_i] += pos_div;

        if (pos_div_[index_i] > this->threshold_by_dimensions_)
        {
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t j = contact_neighborhood.j_[n];
                    if ((*(contact_phi_[k]))[j] > wetting_criterion)
                    {
                        pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
                        break;
                    }
                    else
                    {
                        pos_div_[index_i] = 0.5 * this->threshold_by_dimensions_;
                    }
                }
                if (pos_div_[index_i] == 2.0 * this->threshold_by_dimensions_)
                    break;
            }
        }
    };

  protected:
    Real wetting_criterion;
    StdVec<StdLargeVec<Real> *> contact_phi_;
};

using WettingCoupledSpatialTemporalFreeSurfaceIdentificationComplex =
    SpatialTemporalFreeSurfaceIdentification<NonWettingSurfaceIndication>;
using SpatialTemporalFreeSurfaceIdentificationComplex =
    SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationComplex>;


} // namespace fluid_dynamics
} // namespace SPH
#endif // WETTING_COUPLED_SPATIAL_TEMPORAL_COMPLEX_H