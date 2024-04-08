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
    StdLargeVec<Vecd> elements_neighbors_connection_;
    vector<vector<vector<size_t>>> mesh_topology_;
    Real MinMeshEdge() { return min_distance_between_nodes_; }
    double min_distance_between_nodes_;

    void getDataFromMeshFile3d(const std::string &full_path);
    void getElementCenterCoordinates();
    void getMinimumDistanceBetweenNodes();
    
};
class UnstructuredMesh_3d;
template <> // Base class for generating particles from mesh centroids
class GeneratingMethod<UnstructuredMesh_3d>
{
  public:
    GeneratingMethod(ANSYSMesh_3d& ansys_mesh);
    virtual ~GeneratingMethod(){};

  protected:
    StdLargeVec<Vecd> &elements_centroids_;
    StdLargeVec<Real> &elements_volumes_;
};
template <>
class ParticleGenerator<UnstructuredMesh_3d>
    : public ParticleGenerator<Base>, public GeneratingMethod<UnstructuredMesh_3d>
{
  public:
    ParticleGenerator(SPHBody& sph_body, ANSYSMesh_3d& ansys_mesh);
    virtual ~ParticleGenerator(){};
    virtual void initializeGeometricVariables() override;
};
using ParticleGeneratorUnstructuredMesh_3d = ParticleGenerator<UnstructuredMesh_3d>;


class BaseInnerRelationInFVM_3d : public BaseInnerRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize() override;

  public:
    RealBody *real_body_;
    vector<vector<Real>>& node_coordinates_;
    vector<vector<vector<size_t>>> &mesh_topology_;
    explicit BaseInnerRelationInFVM_3d(RealBody &real_body, ANSYSMesh_3d &ansys_mesh);
    virtual ~BaseInnerRelationInFVM_3d(){};
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
    void createRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const;
    void initializeRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ij,
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
      GhostCreationFromMesh_3d(RealBody &real_body, ANSYSMesh_3d &ansys_mesh, Ghost<ReserveSizeFactor>& ghost_boundary);
      virtual ~GhostCreationFromMesh_3d(){};

      protected:
      Ghost<ReserveSizeFactor>& ghost_boundary_;
      std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
      vector<vector<Real>>& node_coordinates_;
      vector<vector<vector<size_t>>>& mesh_topology_;
      StdLargeVec<Vecd>& pos_;
      StdLargeVec<Real>& Vol_;
      void addGhostParticleAndSetInConfiguration();

    public:
    std::pair<size_t, size_t>& ghost_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
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
    SPHBody &bounds_;
};


class BoundaryConditionSetupInFVM_3d : public fluid_dynamics::FluidDataInner
{
  public:
      BoundaryConditionSetupInFVM_3d(BaseInnerRelationInFVM_3d& inner_relation, GhostCreationFromMesh_3d& ghost_creation);
    virtual ~BoundaryConditionSetupInFVM_3d(){};

    virtual void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyGivenValueInletFlow(size_t ghost_index){};
    virtual void applyOutletBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyTopBoundary(size_t ghost_index, size_t index_i){};
    virtual void applyFarFieldBoundary(size_t ghost_index){};
    virtual void applyPressureOutletBC(size_t ghost_index, size_t index_i){};
    virtual void applySymmetryBoundary(size_t ghost_index, size_t index_i, Vecd e_ij){};
    virtual void applyVelocityInletFlow(size_t ghost_index, size_t index_i) {};
    // Common functionality for resetting boundary conditions
    void resetBoundaryConditions();

  protected:
      StdLargeVec<Real>& rho_, & Vol_, & mass_, & p_;
      StdLargeVec<Vecd>& vel_, & pos_, & mom_;
      std::pair<size_t, size_t>& ghost_bound_;
      vector<vector<size_t>>& each_boundary_type_with_all_ghosts_index_;
      vector<vector<Vecd>>& each_boundary_type_with_all_ghosts_eij_;
      vector<vector<size_t>>& each_boundary_type_contact_real_index_;
};


} // namespace SPH*/
#endif // UNSTRUCTURED_MESH_3D_H