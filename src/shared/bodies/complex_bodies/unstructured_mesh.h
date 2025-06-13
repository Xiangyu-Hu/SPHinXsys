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
 * @file 	unstructured_mesh.h
 * @brief 	Here, we define the common shared classes for FVM.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef UNSTRUCTURED_MESH_H
#define UNSTRUCTURED_MESH_H

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class ANSYSMesh
 * @brief ANASYS mesh.file parser class
 */
class ANSYSMesh
{
  public:
    ANSYSMesh(const std::string &full_path);
    virtual ~ANSYSMesh() {};

    StdVec<size_t> types_of_boundary_condition_;
    StdVec<Vecd> node_coordinates_;
    StdVec<Vecd> elements_centroids_;
    StdVec<Real> elements_volumes_;
    StdVec<StdVec<size_t>> elements_nodes_connection_;
    StdVec<StdVec<StdVec<size_t>>> mesh_topology_;
    Real MinMeshEdge() { return min_distance_between_nodes_; }

  protected:
    Real min_distance_between_nodes_;

    void getDataFromMeshFile(const std::string &full_path);
    void getElementCenterCoordinates();
    void getMinimumDistanceBetweenNodes();
};

/**
 * @class BaseInnerRelationInFVM
 * @brief The abstract relation within a SPH body in FVM
 */
class BaseInnerRelationInFVM : public BaseInnerRelation
{
  public:
    RealBody *real_body_;
    StdVec<Vecd> &node_coordinates_;
    StdVec<StdVec<StdVec<size_t>>> &mesh_topology_;

    explicit BaseInnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh);
    virtual ~BaseInnerRelationInFVM() {};

  protected:
    Vecd *pos_;
    Real *Vol_;
    virtual void resetNeighborhoodCurrentSize() override;
};

/**
 * @class NeighborBuilderInFVM
 * @brief Base neighbor relation between particles i and j.
 */
class NeighborBuilderInFVM
{
  protected:
    //----------------------------------------------------------------------
    //	Below are for constant smoothing length.
    //----------------------------------------------------------------------
    void createRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const;
    void initializeRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ij,
                            Vecd &interface_normal_direction, size_t j_index) const;

  public:
    NeighborBuilderInFVM() {};
    virtual ~NeighborBuilderInFVM() {};
};

/**
 * @class NeighborBuilderInnerInFVM
 * @brief A inner neighbor builder functor in FVM.
 */
class NeighborBuilderInnerInFVM : public NeighborBuilderInFVM
{
  public:
    explicit NeighborBuilderInnerInFVM(SPHBody *body) : NeighborBuilderInFVM() {};
    void operator()(Neighborhood &neighborhood, Real &distance,
                    Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createRelation(neighborhood, distance, dW_ij, interface_normal_direction, j_index)
            : initializeRelation(neighborhood, distance, dW_ij, interface_normal_direction, j_index);
        neighborhood.current_size_++;
    };
};

/** a small functor for obtaining particle index for container index */
struct SPHBodyParticlesIndex
{
    size_t operator()(size_t particle_index) const { return particle_index; };
};

/**
 * @class InnerRelationInFVM
 * @brief The first concrete relation within a SPH body
 */
class InnerRelationInFVM : public BaseInnerRelationInFVM
{
  protected:
    SPHBodyParticlesIndex get_particle_index_;
    NeighborBuilderInnerInFVM get_inner_neighbor_;

  public:
    explicit InnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh);
    virtual ~InnerRelationInFVM() {};

    /** generalized particle search algorithm */
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles,
                                    ParticleConfiguration &particle_configuration,
                                    GetParticleIndex &get_particle_index,
                                    GetNeighborRelation &get_neighbor_relation);
    virtual void updateConfiguration() override;
};

} // namespace SPH
#endif // UNSTRUCTURED_MESH_H
