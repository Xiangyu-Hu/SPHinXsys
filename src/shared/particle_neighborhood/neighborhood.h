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
 * @file 	neighborhood.h
 * @brief 	There are the classes for particle neighborhood.
 * It saves the information for carrying out inter-particle (or pair) interaction,
 * and also considered as the local configuration of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "all_kernels.h"
#include "base_data_type_package.h"
#include "sphinxsys_containers.h"

namespace SPH
{

class SPHBody;
class BodyPart;
class SPHAdaptation;

/**
 * @class Neighborhood
 * @brief A neighborhood around particle i.
 */
class Neighborhood
{
  public:
    size_t current_size_;   /**< the current number of neighbors */
    size_t allocated_size_; /**< the limit of neighbors does not require memory allocation  */

    StdVec<size_t> j_;   /**< index of the neighbor particle. */
    StdVec<Real> W_ij_;  /**< kernel value or particle volume contribution */
    StdVec<Real> dW_ij_; /**< derivative of kernel function or inter-particle surface contribution */
    StdVec<Real> r_ij_;  /**< distance between j and i. */
    StdVec<Vecd> e_ij_;  /**< unit vector pointing from j to i or inter-particle surface direction */

    Neighborhood() : current_size_(0), allocated_size_(0) {};
    ~Neighborhood() {};

    void removeANeighbor(size_t neighbor_n);
};
using ParticleConfiguration = StdVec<Neighborhood>;

/**
 * @class NeighborBuilder
 * @brief Base class for building a neighbor particle j around particles i.
 */
class NeighborBuilder
{
  protected:
    Kernel *kernel_;
    //----------------------------------------------------------------------
    //	Below are for constant smoothing length.
    //----------------------------------------------------------------------
    void createNeighbor(Neighborhood &neighborhood, const Real &distance, const Vecd &displacement, size_t j_index);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance, const Vecd &displacement, size_t j_index);
    //----------------------------------------------------------------------
    //	Below are for variable smoothing length.
    //----------------------------------------------------------------------
    void createNeighbor(Neighborhood &neighborhood, const Real &distance,
                        const Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                            const Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min);
    static Kernel *chooseKernel(SPHBody &body, SPHBody &target_body);

  public:
    NeighborBuilder(Kernel *kernel) : kernel_(kernel) {};
    virtual ~NeighborBuilder() {};
    virtual void operator()(Neighborhood &neighborhood,
                            const Vecd &pos_i, size_t index_i, const ListData &list_data_j) = 0;
};

/**
 * @class NeighborBuilderInner
 * @brief A inner neighbor builder functor.
 */
class NeighborBuilderInner : public NeighborBuilder
{
  public:
    explicit NeighborBuilderInner(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;
};

/**
 * @class NeighborBuilderInnerAdaptive
 * @brief A inner neighbor builder functor when the particles have different smoothing lengths.
 */
class NeighborBuilderInnerAdaptive : public NeighborBuilder
{
  public:
    explicit NeighborBuilderInnerAdaptive(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  protected:
    Real *h_ratio_;
};

/**
 * @class NeighborBuilderSelfContact
 * @brief A self-contact neighbor builder functor.
 */
class NeighborBuilderSelfContact : public NeighborBuilder
{
  public:
    explicit NeighborBuilderSelfContact(SPHBody &body);
    virtual ~NeighborBuilderSelfContact() {};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  protected:
    Vecd *pos0_;
};

/**
 * @class NeighborBuilderContact
 * @brief A contact neighbor builder functor for contact relation.
 */
class NeighborBuilderContact : public NeighborBuilder
{
  public:
    NeighborBuilderContact(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderContact() {};
    virtual void operator()(Neighborhood &neighborhood,
                            const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;
};

/**
 * @class NeighborBuilderSurfaceContact
 * @brief A solid contact neighbor builder functor when bodies having surface contact.
 */
class NeighborBuilderSurfaceContact : public NeighborBuilderContact
{
  private:
    UniquePtrKeeper<Kernel> kernel_keeper_;

  public:
    NeighborBuilderSurfaceContact(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderSurfaceContact() {};
};

/**
 * @class NeighborBuilderContactBodyPart
 * @brief A contact neighbor builder functor when particles j belongs a body part.
 */
class NeighborBuilderContactBodyPart : public NeighborBuilder
{
  public:
    NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part);
    virtual ~NeighborBuilderContactBodyPart() {};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  protected:
    int *part_indicator_; /**< indicator of the body part */
};

/**
 * @class NeighborBuilderContactAdaptive
 * @brief A contact neighbor builder functor when the particles have different smoothing lengths.
 */
class NeighborBuilderContactAdaptive : public NeighborBuilder
{
  public:
    explicit NeighborBuilderContactAdaptive(SPHBody &body, SPHBody &contact_body);
    virtual ~NeighborBuilderContactAdaptive() {};
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  protected:
    SPHAdaptation &adaptation_, &contact_adaptation_;
    Real relative_h_ref_;
};

/**
 * @class BaseNeighborBuilderShell
 * @brief A base neighbor builder functor for contact relation from/to shell.
 */
class BaseNeighborBuilderContactShell : public NeighborBuilder
{
  public:
    explicit BaseNeighborBuilderContactShell(SPHBody &shell_body);

  protected:
    UniquePtrKeeper<Kernel> kernel_keeper_;
    Vecd *n_; // normal direction of contact body
    Real *thickness_;
    Real *k1_ave_;           // 1st principle curvature of contact body
    Real *k2_ave_;           // 2nd principle curvature of contact body
    Real particle_distance_; // reference spacing of contact body

    void createNeighbor(Neighborhood &neighborhood, const Real &distance,
                        size_t index_j, const Real &W_ij,
                        const Real &dW_ij, const Vecd &e_ij);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
                            size_t index_j, const Real &W_ij,
                            const Real &dW_ij, const Vecd &e_ij);
};

/**
 * @class BaseNeighborBuilderFromShell
 * @brief A base neighbor builder functor for solid/shell or fluid/shell contact relation.
 */
class BaseNeighborBuilderContactFromShell : public BaseNeighborBuilderContactShell
{
  public:
    BaseNeighborBuilderContactFromShell(SPHBody &body, SPHBody &contact_body, bool normal_correction);

  protected:
    void update_neighbors(Neighborhood &neighborhood, const Vecd &pos_i, size_t index_i,
                          const ListData &list_data_j);

  private:
    Real direction_corrector_;
};

/**
 * @class NeighborBuilderContactFromShellToFluid
 * @brief A contact neighbor builder functor for contact relation from shell (contact body) to fluid (source body).
 */
class NeighborBuilderContactFromShellToFluid : public BaseNeighborBuilderContactFromShell
{
  public:
    NeighborBuilderContactFromShellToFluid(SPHBody &body, SPHBody &contact_body, bool normal_correction);
    inline void operator()(Neighborhood &neighborhood,
                           const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override
    {
        update_neighbors(neighborhood, pos_i, index_i, list_data_j);
    };
};

/**
 * @class NeighborBuilderContactFromFluidToShell
 * @brief A contact neighbor builder functor for contact relation from shell (contact body) to fluid (source body).
 */
class NeighborBuilderContactFromFluidToShell : public BaseNeighborBuilderContactShell
{
  public:
    NeighborBuilderContactFromFluidToShell(SPHBody &body, SPHBody &contact_body, bool normal_correction);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  private:
    Real direction_corrector_;
};

/**
 * @class ShellNeighborBuilderInnerWithContactKernel
 * @brief A inner neighbor builder functor with reduced kernel.
 * The smoothing length is equal to that of the contact body
 */
class ShellNeighborBuilderInnerWithContactKernel : public NeighborBuilderInner
{
  public:
    explicit ShellNeighborBuilderInnerWithContactKernel(SPHBody &body, SPHBody &contact_body);

  private:
    UniquePtrKeeper<Kernel> kernel_keeper_;
};

/**
 * @class NeighborBuilderShellSelfContact
 * @brief A self-contact neighbor builder functor of shell.
 */
class NeighborBuilderShellSelfContact : public BaseNeighborBuilderContactShell
{
  public:
    explicit NeighborBuilderShellSelfContact(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  private:
    Real *k1_; // 1st principle curvature of contact body
    Real *k2_; // 2nd principle curvature of contact body
    Vecd *pos0_;
};

/**
 * @class NeighborBuilderSurfaceContactFromShell
 * @brief A solid contact neighbor builder functor between solid and a shell
 */
class NeighborBuilderSurfaceContactFromShell : public BaseNeighborBuilderContactFromShell
{
  public:
    NeighborBuilderSurfaceContactFromShell(SPHBody &body, SPHBody &contact_body, bool normal_correction);
    inline void operator()(Neighborhood &neighborhood,
                           const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override
    {
        update_neighbors(neighborhood, pos_i, index_i, list_data_j);
    }
};

/**
 * @class NeighborBuilderSurfaceContactFromSolid
 * @brief A solid contact neighbor builder functor when bodies having surface contact with offset_Wij reduction.
 */
class NeighborBuilderSurfaceContactFromSolid : public NeighborBuilderSurfaceContact
{
  private:
    Real offset_W_ij_;
    void createNeighbor(Neighborhood &neighborhood, const Real &distance, const Vecd &displacement, size_t j_index);
    void initializeNeighbor(Neighborhood &neighborhood, const Real &distance, const Vecd &displacement, size_t j_index);

  public:
    NeighborBuilderSurfaceContactFromSolid(SPHBody &body, SPHBody &contact_body);
    ~NeighborBuilderSurfaceContactFromSolid() override = default;
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;
};

/**
 * @class NeighborBuilderInnerAdaptive
 * @brief A inner neighbor builder functor when the particles have different smoothing lengths for splitting algorithm.
 *        a particle can only see neighbors with ascending ids or higher levels
 */
class NeighborBuilderSplitInnerAdaptive : public NeighborBuilder
{
  public:
    explicit NeighborBuilderSplitInnerAdaptive(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  private:
    Real *h_ratio_;
    int *level_;
};
} // namespace SPH
#endif // NEIGHBORHOOD_H