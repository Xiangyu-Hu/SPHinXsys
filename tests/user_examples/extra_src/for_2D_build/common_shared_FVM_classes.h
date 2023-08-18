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
 *  HU1527/12-1 and HU1527/12-4													*
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
        gerMaximumDistanceBetweenNodes();
    };
    virtual ~readMeshFile(){};

    void getDataFromMeshFile();
    void getElementCenterCoordinates();
    void gerMaximumDistanceBetweenNodes();
    string full_path_;
    vector<size_t> types_of_boundary_condition_;
    vector<vector<Real>> point_coordinates_2D_;
    vector<vector<Real>> point_coordinates;
    StdLargeVec<Vecd> elements_center_coordinates_;
    StdLargeVec<Real> elements_volumes_;
    vector<vector<size_t>> elements_nodes_connection_;
    StdLargeVec<Vec3d> elements_neighbors_connection_;
    vector<vector<vector<size_t>>> cell_lists_;
    double max_distance_between_nodes_;
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
    void searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles,
                                    ParticleConfiguration &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation);
    virtual void updateConfiguration() override;
};

/**
 * @class BaseGhostCreation
 * @brief Base class for the ghost particle
 */
class GhostCreationFromMesh : public GeneralDataDelegateSimple
{
  public:
    GhostCreationFromMesh(RealBody &real_body, vector<vector<vector<size_t>>> &data_inpute, vector<vector<Real>> nodes_coordinates)
        : GeneralDataDelegateSimple(real_body), all_needed_data_from_mesh_file_(data_inpute), nodes_coordinates_(nodes_coordinates),
          pos_(particles_->pos_), Vol_(particles_->Vol_), total_ghost_particles_(particles_->total_ghost_particles_),
          real_particles_bound_(particles_->real_particles_bound_)
    {
        each_boundary_type_with_all_ghosts_index_.resize(50);
        each_boundary_type_with_all_ghosts_eij_.resize(50);
        each_boundary_type_contact_real_index_.resize(50);
        ghost_particles_.resize(1);
        addGhostParticleAndSetInConfiguration();
    }
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

    void addGhostParticleAndSetInConfiguration()
    {
        for (size_t i = 0; i != ghost_particles_.size(); ++i)
            ghost_particles_[i].clear();

        for (size_t index_i = 0; index_i != real_particles_bound_; ++index_i)
        {
            for (size_t neighbor_index = 0; neighbor_index != all_needed_data_from_mesh_file_[index_i].size(); ++neighbor_index)
            {
                size_t boundary_type = all_needed_data_from_mesh_file_[index_i][neighbor_index][1];
                if (all_needed_data_from_mesh_file_[index_i][neighbor_index][1] != 2)
                {
                    mutex_create_ghost_particle_.lock();
                    size_t ghost_particle_index = particles_->insertAGhostParticle(index_i);
                    size_t node1_index = all_needed_data_from_mesh_file_[index_i][neighbor_index][2];
                    size_t node2_index = all_needed_data_from_mesh_file_[index_i][neighbor_index][3];
                    Vecd node1_position = Vecd(nodes_coordinates_[node1_index][0], nodes_coordinates_[node1_index][1]);
                    Vecd node2_position = Vecd(nodes_coordinates_[node2_index][0], nodes_coordinates_[node2_index][1]);
                    Vecd ghost_particle_position = 0.5 * (node1_position + node2_position);

                    all_needed_data_from_mesh_file_[index_i][neighbor_index][0] = ghost_particle_index + 1;
                    ghost_particles_[0].push_back(ghost_particle_index);
                    pos_[ghost_particle_index] = ghost_particle_position;
                    mutex_create_ghost_particle_.unlock();

                    all_needed_data_from_mesh_file_.resize(ghost_particle_index);
                    std::vector<std::vector<size_t>> new_element;

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element1 = {index_i + 1, boundary_type, node1_index, node2_index};
                    new_element.push_back(sub_element1);

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element2 = {index_i + 1, boundary_type, node1_index, node2_index};
                    new_element.push_back(sub_element2);

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element3 = {index_i + 1, boundary_type, node1_index, node2_index};
                    new_element.push_back(sub_element3);

                    // Add the new element to all_needed_data_from_mesh_file_
                    all_needed_data_from_mesh_file_.push_back(new_element);
                    // all_needed_data_from_mesh_file_[ghost_particle_index][0][0].push_back(size_t(0);

                    // creating the boundary files with ghost particle index
                    each_boundary_type_with_all_ghosts_index_[boundary_type].push_back(ghost_particle_index);

                    // creating the boundary files with contact real particle index
                    each_boundary_type_contact_real_index_[boundary_type].push_back(index_i);

                    // creating the boundary files with ghost eij
                    Vecd interface_area_vector = node1_position - node2_position;
                    Real interface_area_size = interface_area_vector.norm();
                    Vecd unit_vector = interface_area_vector / interface_area_size;
                    // normal unit vector
                    Vecd normal_vector = Vecd(unit_vector[1], -unit_vector[0]);
                    // judge the direction
                    Vecd particle_position = pos_[index_i];
                    Vecd node1_to_center_direction = particle_position - node1_position;
                    if (node1_to_center_direction.dot(normal_vector) < 0)
                    {
                        normal_vector = -normal_vector;
                    };
                    each_boundary_type_with_all_ghosts_eij_[boundary_type].push_back(normal_vector);
                }
            }
        }
    };
};

} // namespace SPH
#endif // COMMON_SHARED_FVM_CLASSES_H