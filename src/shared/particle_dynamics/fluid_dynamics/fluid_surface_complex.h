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

#ifndef FLUID_SURFACE_COMPLEX_H
#define FLUID_SURFACE_COMPLEX_H

#include "fluid_dynamics_complex.hpp"
#include "fluid_surface_inner.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FreeSurfaceIndicationComplex
 * @brief indicate the particles near the free fluid surface.
 */
class FreeSurfaceIndicationComplex : public FreeSurfaceIndicationInner, public FluidContactData
{
  public:
    FreeSurfaceIndicationComplex(BaseInnerRelation &inner_relation,
                                 BaseContactRelation &contact_relation, Real threshold = 0.75);
    explicit FreeSurfaceIndicationComplex(ComplexRelation &complex_relation, Real threshold = 0.75);
    virtual ~FreeSurfaceIndicationComplex(){};

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
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
};
using SpatialTemporalFreeSurfaceIdentificationComplex =
    SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationComplex>;

/** the cases with free surface and freestream */
using DensitySummationFreeSurfaceComplex = DensitySummationFreeSurface<DensitySummationComplex>;
using DensitySummationFreeStreamComplex = DensitySummationFreeStream<DensitySummationFreeSurfaceComplex>;
/** the case with variable smoothing length */
using DensitySummationFreeSurfaceComplexAdaptive = DensitySummationFreeSurface<DensitySummationComplexAdaptive>;
using DensitySummationFreeStreamComplexAdaptive = DensitySummationFreeStream<DensitySummationFreeSurfaceComplexAdaptive>;

/**
 * @class ColorFunctionGradientComplex
 * @brief indicate the particles near the free fluid surface.
 */
class ColorFunctionGradientComplex : public ColorFunctionGradientInner, public FluidContactData
{
  public:
    ColorFunctionGradientComplex(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation);
    ColorFunctionGradientComplex(ComplexRelation &complex_relation);
    virtual ~ColorFunctionGradientComplex(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        ColorFunctionGradientInner::interaction(index_i, dt);

        Vecd gradient = Vecd::Zero();
        if (pos_div_[index_i] < threshold_by_dimensions_)
        {
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    gradient -= contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
                }
            }
        }
        color_grad_[index_i] += gradient;
        surface_norm_[index_i] = color_grad_[index_i] / (color_grad_[index_i].norm() + TinyReal);
    };

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
};

/**
 * @class 	SurfaceNormWithWall
 * @brief  Modify surface norm when contact with wall
 */
class SurfaceNormWithWall : public LocalDynamics, public FSIContactData
{
  public:
    SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle);
    virtual ~SurfaceNormWithWall(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real large_dist(1.0e6);
        Vecd n_i = surface_norm_[index_i];
        Real smoothing_factor(1.0);
        Vecd smooth_norm = Vecd::Zero();
        Vecd n_i_w = Vecd::Zero();
        /** Contact interaction. */
        if (surface_indicator_[index_i] == 1)
        {
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
                Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                {
                    size_t index_j = wall_neighborhood.j_[n];
                    if (wall_neighborhood.r_ij_[n] < large_dist)
                    {
                        Vecd n_w_t = n_i - n_i.dot(n_k[index_j]) * n_k[index_j];
                        Vecd n_t = n_w_t / (n_w_t.norm() + TinyReal);
                        n_i_w = n_t * sin(contact_angle_) + cos(contact_angle_) * n_k[index_j];
                        /** No change for multi-resolution. */
                        Real r_ij = wall_neighborhood.r_ij_[n] * n_k[index_j].dot(wall_neighborhood.e_ij_[n]);
                        if (r_ij <= smoothing_length_)
                        {
                            smoothing_factor = 0.0;
                        }
                        else
                        {
                            smoothing_factor = (r_ij - smoothing_length_) / smoothing_length_;
                        }
                        large_dist = wall_neighborhood.r_ij_[n];
                        smooth_norm = smoothing_factor * n_i + (1.0 - smoothing_factor) * n_i_w;
                        surface_norm_[index_i] = smooth_norm / (smooth_norm.norm() + TinyReal);
                    }
                }
            }
        }
    };

  protected:
    Real contact_angle_;
    Real smoothing_length_;
    Real particle_spacing_;
    StdLargeVec<int> &surface_indicator_;
    StdLargeVec<Vecd> &surface_norm_;
    StdLargeVec<Real> &pos_div_;
    StdVec<StdLargeVec<Vecd> *> wall_n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_SURFACE_COMPLEX_H