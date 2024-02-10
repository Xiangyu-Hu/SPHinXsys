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

#include "base_fluid_dynamics.h"
#include "base_particle_generator.h"
#include "compressible_fluid.h"
#include "fluid_body.h"
#include "io_vtk.h"
using namespace std;
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
    virtual ~ANSYSMesh(){};

    StdVec<size_t> types_of_boundary_condition_;
    StdLargeVec<Vecd> node_coordinates_;
    StdLargeVec<Vecd> elements_centroids_;
    StdLargeVec<Real> elements_volumes_;
    StdLargeVec<StdVec<size_t>> elements_nodes_connection_;
    vector<vector<vector<size_t>>> mesh_topology_;
    double min_distance_between_nodes_;

  protected:
    void readNodeCoordinate(const std::string &text_line, StdLargeVec<Vec2d> &node_coordinates);
    void readNodeCoordinate(const std::string &text_line, StdLargeVec<Vec3d> &node_coordinates);
    void getDataFromMeshFile(const std::string &full_path);
    void getElementCenterCoordinates();
    void gerMinimumDistanceBetweenNodes();
};

/**
 * @class BaseInnerRelationInFVM
 * @brief The abstract relation within a SPH body in FVM
 */
class BaseInnerRelationInFVM : public BaseInnerRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize() override;

  public:
    RealBody *real_body_;
    StdLargeVec<Vecd> &node_coordinates_;
    vector<vector<vector<size_t>>> &mesh_topology_;

    explicit BaseInnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh);
    virtual ~BaseInnerRelationInFVM(){};
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
    void createRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const;
    void initializeRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j,
                            Vecd &interface_normal_direction, size_t j_index) const;

  public:
    NeighborBuilderInFVM(){};
    virtual ~NeighborBuilderInFVM(){};
};

/**
 * @class NeighborBuilderInnerInFVM
 * @brief A inner neighbor builder functor in FVM.
 */
class NeighborBuilderInnerInFVM : public NeighborBuilderInFVM
{
  public:
    explicit NeighborBuilderInnerInFVM(SPHBody *body) : NeighborBuilderInFVM(){};
    void operator()(Neighborhood &neighborhood, Real &distance,
                    Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createRelation(neighborhood, distance, dW_ijV_j, interface_normal_direction, j_index)
            : initializeRelation(neighborhood, distance, dW_ijV_j, interface_normal_direction, j_index);
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
    virtual ~InnerRelationInFVM(){};

    /** generalized particle search algorithm */
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles,
                                    ParticleConfiguration &particle_configuration,
                                    GetParticleIndex &get_particle_index,
                                    GetNeighborRelation &get_neighbor_relation);
    virtual void updateConfiguration() override;
};

/**
 * @class BaseGhostCreation
 * @brief Base class for the ghost particle
 */
class GhostCreationFromMesh : public GeneralDataDelegateSimple
{
  public:
    GhostCreationFromMesh(RealBody &real_body, Ghost<UnstructuredMesh> &ghost_boundary);
    virtual ~GhostCreationFromMesh(){};

  protected:
    Ghost<UnstructuredMesh> &ghost_boundary_;
    ANSYSMesh &ansys_mesh_;
    std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
    StdLargeVec<Vecd> &node_coordinates_;
    vector<vector<vector<size_t>>> &mesh_topology_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &Vol_;
    void addGhostParticleAndSetInConfiguration();

  public:
    std::pair<size_t, size_t> &ghost_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};

/**
 * @class BodyStatesRecordingInMeshToVtp
 * @brief  Write files for bodies
 * the output file is VTK XML format in FVMcan visualized by ParaView the data type vtkPolyData
 */
class BodyStatesRecordingInMeshToVtp : public BodyStatesRecording
{
  public:
    BodyStatesRecordingInMeshToVtp(SPHBody &body, ANSYSMesh &ansys_mesh);
    virtual ~BodyStatesRecordingInMeshToVtp(){};

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    StdLargeVec<Vecd> &node_coordinates_;
    StdLargeVec<StdVec<size_t>> &elements_nodes_connection_;
};

//----------------------------------------------------------------------
//	BoundaryConditionSetupInFVM
//----------------------------------------------------------------------
class BoundaryConditionSetupInFVM : public fluid_dynamics::FluidDataInner
{
  public:
    BoundaryConditionSetupInFVM(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation);
    virtual ~BoundaryConditionSetupInFVM(){};
    virtual void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyGivenValueInletFlow(size_t ghost_index){};
    virtual void applyOutletBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyTopBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyFarFieldBoundary(size_t ghost_index){};
    // Common functionality for resetting boundary conditions
    void resetBoundaryConditions();

  protected:
    StdLargeVec<Real> &rho_, &Vol_, &mass_, &p_;
    StdLargeVec<Vecd> &vel_, &pos_, &mom_;
    std::pair<size_t, size_t> &ghost_bound_;
    vector<vector<size_t>> &each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> &each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> &each_boundary_type_contact_real_index_;
};
} // namespace SPH
#endif // UNSTRUCTURED_MESH_H