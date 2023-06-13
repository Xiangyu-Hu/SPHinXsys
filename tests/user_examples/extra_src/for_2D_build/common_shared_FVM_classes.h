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
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
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

#include "fluid_body.h"
#include "general_dynamics.h"
#include "fluid_dynamics_inner.h"
#include "compressible_fluid.h"
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
        };
        virtual ~readMeshFile(){};

        void getDataFromMeshFile();
        //void getNodesOfElementFromCellLists();
        void getElementCenterCoordinates();
        string full_path_;
        vector<size_t> types_of_boundary_condition_;
        vector<vector<Real>> point_coordinates_2D_;
        //StdLargeVec<Vec3d> point_coordinates_3D_;
        vector<vector<Real>> point_coordinates;
        //StdLargeVec<Vec3d> elements_nodes_connection_;
        StdLargeVec<Vecd> elements_center_coordinates_;
        StdLargeVec<Real> elements_volumes_;
        vector<vector<size_t>> elements_nodes_connection_;
        StdLargeVec<Vec3d> elements_neighbors_connection_;
        vector<vector<vector<size_t>>> cell_lists_;
    };

    /**
    * @class NeighborhoodInFVM
    * @brief A neighborhood around particle i in FVM.
    */
    class NeighborhoodInFVM : public Neighborhood
    {
    public:
	    StdLargeVec<size_t> boundary_type_;	  /**< index of the neighbor particle. */
	    StdLargeVec<Real> interface_size_; /**< interface area size */

	    NeighborhoodInFVM() : Neighborhood() {};
	    ~NeighborhoodInFVM() {};

	    void removeANeighbor(size_t neighbor_n);
    };

    /** Neighborhoods for all particles in a body for a inner body relation. */
    using ParticleConfigurationInFVM = StdLargeVec<NeighborhoodInFVM>;

    /**
    * @class BaseInnerRelationInFVM
    * @brief The abstract relation within a SPH body in FVM
    */
    class BaseInnerRelationInFVM : public SPHRelation
    {
      protected:
        virtual void resetNeighborhoodCurrentSize();

      public:
        RealBody *real_body_;
        vector<vector<vector<size_t>>> all_needed_data_from_mesh_file_;
        vector<vector<Real>> nodes_coordinates_;
        ParticleConfigurationInFVM inner_configuration_in_FVM_; /**< inner configuration for the neighbor relations. */
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
        void createRelation(NeighborhoodInFVM &neighborhood, Real &distance, Real &dW_ijV_j, Real &interface_size,
                            Vecd &interface_normal_direction, size_t bc_type, size_t j_index) const;
        void initializeRelation(NeighborhoodInFVM &neighborhood, Real &distance, Real &dW_ijV_j, Real &interface_size,
                                Vecd &interface_normal_direction, size_t bc_type, size_t j_index) const;
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
        void operator()(NeighborhoodInFVM &neighborhood, Real &distance,
                        Real &dW_ijV_j, Real &interface_size, Vecd &interface_normal_direction, size_t bc_type, size_t j_index) const
        {
            neighborhood.current_size_ >= neighborhood.allocated_size_
                ? createRelation(neighborhood, distance, dW_ijV_j, interface_size, interface_normal_direction, bc_type, j_index)
                : initializeRelation(neighborhood, distance, dW_ijV_j, interface_size, interface_normal_direction, bc_type, j_index);
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
                                        ParticleConfigurationInFVM &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation);
        virtual void updateConfiguration() override;
    };

    /**
     * @class DataDelegateInnerInFVM
     * @brief prepare data for inner particle dynamics
     */
    template <class ParticlesType = BaseParticles,
              class BaseDataDelegateType = DataDelegateSimple<ParticlesType>>
    class DataDelegateInnerInFVM : public BaseDataDelegateType
    {
      public:
        explicit DataDelegateInnerInFVM(BaseInnerRelationInFVM &inner_relation)
            : BaseDataDelegateType(inner_relation.getSPHBody()),
              inner_configuration_in_FVM_(inner_relation.inner_configuration_in_FVM_){};
        virtual ~DataDelegateInnerInFVM(){};

      protected:
        /** inner configuration of the designated body */
        ParticleConfigurationInFVM &inner_configuration_in_FVM_;
    };

}
#endif // COMMON_SHARED_FVM_CLASSES_H