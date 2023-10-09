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
 * @file 	common_shared_FVM_classes.h
 * @brief 	Here, we define the common shared classes for FVM.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef COMMON_SHARED_FVM_CLASSES_H
#define COMMON_SHARED_FVM_CLASSES_H

#include "base_particle_generator.h"
#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "io_vtk.h"
using namespace std;
namespace SPH
{
/**
 * @class readMeshFile
 * @brief ANASYS mesh.file parser class
 */
class readMeshFile
{
  public:
    readMeshFile(std::string full_path)
    {
        full_path_ = full_path;
        getDataFromMeshFile();
        getElementCenterCoordinates();
        gerMinimumDistanceBetweenNodes();
    };
    virtual ~readMeshFile(){};

    void getDataFromMeshFile();
    void getElementCenterCoordinates();
    void gerMinimumDistanceBetweenNodes();
    string full_path_;
    vector<size_t> types_of_boundary_condition_;
    vector<vector<Real>> point_coordinates_2D_;
    vector<vector<Real>> point_coordinates;
    StdLargeVec<Vecd> elements_center_coordinates_;
    StdLargeVec<Real> elements_volumes_;
    vector<vector<size_t>> elements_nodes_connection_;
    StdLargeVec<Vec3d> elements_neighbors_connection_;
    vector<vector<vector<size_t>>> cell_lists_;
    double min_distance_between_nodes_;
};

/**
 * @class BaseInnerRelationInFVM
 * @brief The abstract relation within a SPH body in FVM
 */
class BaseInnerRelationInFVM : public BaseInnerRelation
{
  protected:
    virtual void resetNeighborhoodCurrentSize();

  public:
    RealBody *real_body_;
    vector<vector<vector<size_t>>> all_needed_data_from_mesh_file_;
    vector<vector<Real>> nodes_coordinates_;
    explicit BaseInnerRelationInFVM(RealBody &real_body, vector<vector<vector<size_t>>> data_inpute, vector<vector<Real>> nodes_coordinates);
    virtual ~BaseInnerRelationInFVM(){};

    virtual void resizeConfiguration() override;
};

/**
 * @class ParticleGeneratorInFVM
 * @brief Generate particle directly from position-and-volume data.
 */
class ParticleGeneratorInFVM : public ParticleGenerator
{
  public:
    explicit ParticleGeneratorInFVM(SPHBody &sph_body)
        : ParticleGenerator(sph_body){};
    ParticleGeneratorInFVM(SPHBody &sph_body, const StdLargeVec<Vecd> &positions, const StdLargeVec<Real> &elements_volumes);
    virtual ~ParticleGeneratorInFVM(){};
    /** Initialize geometrical variable for observe particles. */
    virtual void initializeGeometricVariables() override;

  protected:
    StdLargeVec<Vecd> elements_center_coordinates_;
    StdLargeVec<Real> elements_volumes_;
};

/**
 * @class NeighborBuilderInFVM
 * @brief Base neighbor relation between particles i and j.
 */
class NeighborBuilderInFVM
{
  protected:
    Kernel *kernel_;
    //----------------------------------------------------------------------
    //	Below are for constant smoothing length.
    //----------------------------------------------------------------------
    void createRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const;
    void initializeRelation(Neighborhood &neighborhood, Real &distance, Real &dW_ijV_j,
                            Vecd &interface_normal_direction, size_t j_index) const;

  public:
    NeighborBuilderInFVM() : kernel_(nullptr){};
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
    explicit InnerRelationInFVM(RealBody &real_body, vector<vector<vector<size_t>>> data_inpute, vector<vector<Real>> nodes_coordinates);
    virtual ~InnerRelationInFVM(){};

    /** generalized particle search algorithm */
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles, ParticleConfiguration &particle_configuration, 
        GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation);
    virtual void updateConfiguration() override;
};

/**
 * @class BaseGhostCreation
 * @brief Base class for the ghost particle
 */
class GhostCreationFromMesh : public GeneralDataDelegateSimple
{
  public:
    GhostCreationFromMesh(RealBody &real_body, vector<vector<vector<size_t>>> &data_inpute, vector<vector<Real>> nodes_coordinates);
    virtual ~GhostCreationFromMesh(){};
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;

  protected:
    std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
    vector<vector<vector<size_t>>> &all_needed_data_from_mesh_file_;
    vector<vector<Real>> nodes_coordinates_;
    StdLargeVec<Vecd> &pos_;
    StdVec<IndexVector> ghost_particles_;
    StdLargeVec<Real> &Vol_;
    size_t &total_ghost_particles_;
    size_t &real_particles_bound_;
    void addGhostParticleAndSetInConfiguration();
};

/**
 * @class BodyStatesRecordingInMeshToVtp
 * @brief  Write files for bodies
 * the output file is VTK XML format in FVMcan visualized by ParaView the data type vtkPolyData
 */
class BodyStatesRecordingInMeshToVtp : public BodyStatesRecording
{
  public:
    BodyStatesRecordingInMeshToVtp(IOEnvironment &io_environment, SPHBody &body, 
        vector<vector<size_t>> elements_nodes_connection, vector<vector<Real>> nodes_coordinates)
        : BodyStatesRecording(io_environment, body), elements_nodes_connection_(elements_nodes_connection), nodes_coordinates_(nodes_coordinates){};
    BodyStatesRecordingInMeshToVtp(IOEnvironment &io_environment, SPHBodyVector bodies, 
        vector<vector<size_t>> elements_nodes_connection, vector<vector<Real>> nodes_coordinates)
        : BodyStatesRecording(io_environment, bodies), elements_nodes_connection_(elements_nodes_connection), nodes_coordinates_(nodes_coordinates){};
    virtual ~BodyStatesRecordingInMeshToVtp(){};

  protected:
    virtual void writeWithFileName(const std::string &sequence) override;
    vector<vector<size_t>> elements_nodes_connection_;
    vector<vector<Real>> nodes_coordinates_;
};

//----------------------------------------------------------------------
//	BoundaryConditionSetupInFVM
//----------------------------------------------------------------------
class BoundaryConditionSetupInFVM : public fluid_dynamics::FluidDataInner
{
  public:
    BoundaryConditionSetupInFVM(BaseInnerRelationInFVM &inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
                                vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index);
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
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_, &pos_;
    size_t &total_ghost_particles_;
    size_t &real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};
} // namespace SPH
#endif // COMMON_SHARED_FVM_CLASSES_H