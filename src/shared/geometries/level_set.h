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
 * @file 	level_set.h
 * @brief 	Level set is a function which is defined as signed distance to a surface or interface.
 * @author	Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */

#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include "adaptation.h"
#include "base_geometry.h"
#include "mesh_with_data_packages.hpp"

namespace SPH
{
/**
 * @class BaseLevelSet
 * @brief A abstract describes a level set field defined on a mesh.
 * Level set is a signed distance function to an interface where the zero level set locates.
 * Here, the region with negative level set is considered as the region enclose by the interface.
 */
class BaseLevelSet : public BaseMeshField
{
  public:
    BaseLevelSet(Shape &shape, SPHAdaptation &sph_adaptation);
    virtual ~BaseLevelSet(){};

    virtual void cleanInterface(Real small_shift_factor) = 0;
    virtual void correctTopology(Real small_shift_factor) = 0;
    virtual bool probeIsWithinMeshBound(const Vecd &position) = 0;
    virtual Real probeSignedDistance(const Vecd &position) = 0;
    virtual Vecd probeNormalDirection(const Vecd &position) = 0;
    virtual Vecd probeLevelSetGradient(const Vecd &position) = 0;
    virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) = 0;
    virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) = 0;

  protected:
    Shape &shape_; /**< the geometry is described by the level set. */
    SPHAdaptation &sph_adaptation_;

    /** a cut cell is a cut by the level set. */
    /** "Multi-scale modeling of compressible multi-fluid flows with conservative interface method."
     * Hu, X. Y., et al., Proceedings of the Summer Program. Vol. 301. Stanford, CA, USA:
     * Center for Turbulence Research, Stanford University, 2010.*/
    Real CutCellVolumeFraction(Real phi, const Vecd &phi_gradient, Real data_spacing);
};

/**
 * @class LevelSet
 * @brief Mesh with level set data as packages.
 * Note that the mesh containing the data packages are cell-based
 * but within the data package, the data is grid-based.
 * Note that the level set data is initialized after the constructor.
 */
class LevelSet : public MeshWithGridDataPackages<GridDataPackage<4, 1>>,
                 public BaseLevelSet
{
  public:
    typedef GridDataPackage<4, 1> LevelSetDataPackage;
    ConcurrentVec<LevelSetDataPackage *> core_data_pkgs_; /**< packages near to zero level set. */
    Real global_h_ratio_;

    /** This constructor only initialize far field. */
    LevelSet(BoundingBox tentative_bounds, Real data_spacing, size_t buffer_size, Shape &shape, SPHAdaptation &sph_adaptation);
    /** This constructor generate inner packages too. */
    LevelSet(BoundingBox tentative_bounds, Real data_spacing, Shape &shape, SPHAdaptation &sph_adaptation);
    virtual ~LevelSet(){};

    virtual void cleanInterface(Real small_shift_factor) override;
    virtual void correctTopology(Real small_shift_factor) override;
    virtual bool probeIsWithinMeshBound(const Vecd &position) override;
    virtual Real probeSignedDistance(const Vecd &position) override;
    virtual Vecd probeNormalDirection(const Vecd &position) override;
    virtual Vecd probeLevelSetGradient(const Vecd &position) override;
    virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) override;
    virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) override;
    virtual void writeMeshFieldToPlt(std::ofstream &output_file) override;
    bool isWithinCorePackage(Vecd position);
    Real computeKernelIntegral(const Vecd &position);
    Vecd computeKernelGradientIntegral(const Vecd &position);

  protected:
    MeshVariable<Real> &phi_;
    MeshVariable<int> &near_interface_id_;
    MeshVariable<Vecd> &phi_gradient_;
    MeshVariable<Real> &kernel_weight_;
    MeshVariable<Vecd> &kernel_gradient_;
    Kernel &kernel_;

    void initializeDataForSingularPackage(LevelSetDataPackage *data_pkg, Real far_field_level_set);
    void initializeBasicDataForAPackage(LevelSetDataPackage *data_pkg, Shape &shape);
    void redistanceInterfaceForAPackage(LevelSetDataPackage *core_data_pkg);

    void finishDataPackages();
    void reinitializeLevelSet();
    void markNearInterface(Real small_shift_factor);
    void redistanceInterface();
    void diffuseLevelSetSign();
    void updateLevelSetGradient();
    void updateKernelIntegrals();
    bool isInnerPackage(const Arrayi &cell_index);
    void initializeDataInACell(const Arrayi &cell_index);
    void initializeAddressesInACell(const Arrayi &cell_index);
    void tagACellIsInnerPackage(const Arrayi &cell_index);

    // upwind algorithm choosing candidate difference by the sign
    Real upwindDifference(Real sign, Real df_p, Real df_n);
};

/**
 * @class RefinedLevelSet
 * @brief level set  which has double resolution of a coarse level set.
 */
class RefinedLevelSet : public RefinedMesh<LevelSet>
{
  public:
    RefinedLevelSet(BoundingBox tentative_bounds, LevelSet &coarse_level_set, Shape &shape, SPHAdaptation &sph_adaptation);
    virtual ~RefinedLevelSet(){};

  protected:
    void initializeDataInACellFromCoarse(const Arrayi &cell_index);
};

/**
 * @class MultilevelLevelSet
 * @brief Defining a multilevel level set for a complex region.
 */
class MultilevelLevelSet : public MultilevelMesh<BaseLevelSet, LevelSet, RefinedLevelSet>
{
  public:
    MultilevelLevelSet(BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation);
    virtual ~MultilevelLevelSet(){};

    virtual void cleanInterface(Real small_shift_factor) override;
    virtual void correctTopology(Real small_shift_factor) override;
    virtual bool probeIsWithinMeshBound(const Vecd &position) override;
    virtual Real probeSignedDistance(const Vecd &position) override;
    virtual Vecd probeNormalDirection(const Vecd &position) override;
    virtual Vecd probeLevelSetGradient(const Vecd &position) override;
    virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) override;
    virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) override;

  protected:
    inline size_t getProbeLevel(const Vecd &position);
    inline size_t getCoarseLevel(Real h_ratio);
};
} // namespace SPH
#endif // LEVEL_SET_H
