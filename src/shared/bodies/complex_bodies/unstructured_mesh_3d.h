/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2023 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	unstructured_mesh_3d.h
 * @brief 	Here, we define the common shared classes for FVM.
 * @author	Yash Mandaokar, Zhentong Wang and Xiangyu Hu
 */
#ifndef UNSTRUCTURED_MESH_3D_H
#define UNSTRUCTURED_MESH_3D_H

#include "base_particle_generator.h"
#include "compressible_fluid.h"
#include "fluid_body.h"
#include "io_vtk.h"
#include <Eigen/Dense>
using namespace std;
namespace SPH
{
/**
 * @class ANSYSMeshFile
 * @brief ANASYS mesh.file parser class
 */
class ANSYSMesh_3d
{
  public:
    ANSYSMesh_3d(const std::string &full_path);
    virtual ~ANSYSMesh_3d(){};

    string full_path_;
    vector<size_t> types_of_boundary_condition_;
    vector<vector<Real>> node_coordinates_;
    vector<vector<Real>> point_coordinates;
    StdLargeVec<Vecd> elements_centroids_;
    StdLargeVec<Real> elements_volumes_;
    vector<vector<size_t>> elements_nodes_connection_;
    StdLargeVec<Vec3d> elements_neighbors_connection_;
    vector<vector<vector<size_t>>> mesh_topology_;
    double min_distance_between_nodes_;

    void getDataFromMeshFile3d(const std::string &full_path);
    void getElementCenterCoordinates();
    void getMinimumDistanceBetweenNodes();
    
};


class BaseInnerRelationInFVM_3d : public BaseInnerRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize() override;

  public:
    RealBody *real_body_;
    vector<vector<vector<size_t>>> &mesh_topology_;
    vector<vector<Real>> &node_coordinates_;
    explicit BaseInnerRelationInFVM_3d(RealBody &real_body, ANSYSMesh_3d &ansys_mesh);
    virtual ~BaseInnerRelationInFVM_3d(){};

    virtual void resizeConfiguration() override;
};

/**
 * @class NeighborBuilderInFVM
 * @brief Base neighbor relation between particles i and j.*/
 
class NeighborBuilderInFVM_3d
{
  protected:
    //Kernel *kernel_;
    //----------------------------------------------------------------------
    //	Below are for constant smoothing length.
    //----------------------------------------------------------------------
    void createRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const;
    void initializeRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j,
                            Vecd &interface_normal_direction, size_t j_index) const;

  public:
      NeighborBuilderInFVM_3d(){} //: kernel_(nullptr){};
    virtual ~NeighborBuilderInFVM_3d(){};
};

/**
 * @class NeighborBuilderInnerInFVM
 * @brief A inner neighbor builder functor in FVM.*/
 
class NeighborBuilderInnerInFVM_3d : public NeighborBuilderInFVM_3d
{
  public:
    explicit NeighborBuilderInnerInFVM_3d(SPHBody *body) : NeighborBuilderInFVM_3d(){};
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
struct SPHBodyParticlesIndex_3d
{
    size_t operator()(size_t particle_index) const { return particle_index; };
};

/**
 * @class InnerRelationInFVM
 * @brief The first concrete relation within a SPH body*/
 
class InnerRelationInFVM_3d : public BaseInnerRelationInFVM_3d
{
  protected:
    SPHBodyParticlesIndex_3d get_particle_index_;
    NeighborBuilderInnerInFVM_3d get_inner_neighbor_;

  public:
    explicit InnerRelationInFVM_3d(RealBody &real_body, ANSYSMesh_3d &ansys_mesh);
    virtual ~InnerRelationInFVM_3d(){};

    /** generalized particle search algorithm */
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles,
                                    ParticleConfiguration &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation);
    virtual void updateConfiguration() override;
};


/**
 * @class BaseGhostCreation
 * @brief Base class for the ghost particle
 */
class GhostCreationFromMesh_3d : public GeneralDataDelegateSimple
{
  public:
      GhostCreationFromMesh_3d(RealBody &real_body, ANSYSMesh_3d &ansys_mesh);
      virtual ~GhostCreationFromMesh_3d(){};
      vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
      vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
      vector<vector<size_t>> each_boundary_type_contact_real_index_;

      protected:
      std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
      vector<vector<Real>> &node_coordinates_;
      vector<vector<vector<size_t>>> &mesh_topology_;
      StdLargeVec<Vecd> &pos_;
      StdVec<IndexVector> ghost_particles_;
      StdLargeVec<Real> &Vol_;
      size_t &total_ghost_particles_;
      size_t &real_particles_bound_;
      void addGhostParticleAndSetInConfiguration();
};

class BodyStatesRecordingInMeshToVtu : public BodyStatesRecording
{
  public:
    BodyStatesRecordingInMeshToVtu(SPHBody &body, ANSYSMesh_3d &ansys_mesh);
    virtual ~BodyStatesRecordingInMeshToVtu(){};

    protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    vector<vector<Real>> &node_coordinates_;
    vector<vector<size_t>> &elements_nodes_connection_;
};


class BoundaryConditionSetupInFVM_3d : public fluid_dynamics::FluidDataInner
{
  public:
      BoundaryConditionSetupInFVM_3d(BaseInnerRelationInFVM_3d& inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                                vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index);
    virtual ~BoundaryConditionSetupInFVM_3d(){};

    virtual void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyGivenValueInletFlow(size_t ghost_index){};
    virtual void applyOutletBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyTopBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyFarFieldBoundary(size_t ghost_index){};
    virtual void applyPressureOutletBC(size_t ghost_index, size_t index_i){};
    virtual void applysymmetry(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyVelocityInletFlow(size_t ghost_index, size_t index_i) {};
    // Common functionality for resetting boundary conditions
    void resetBoundaryConditions();

  protected:
    StdLargeVec<Real> &rho_, &Vol_, &mass_, &p_;
    StdLargeVec<Vecd> &vel_, &pos_, &mom_;
    size_t &total_ghost_particles_;
    size_t &real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};


} // namespace SPH*/
#endif // UNSTRUCTURED_MESH_3D_H