/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	level_set.h
* @brief 	This is the base classes of mesh, which describe ordered and indexed
*			data sets.  Depending on application, there are different data 
* 			saved on the mesh. The intersection points of mesh lines are called 
*			grid points, the element enclosed by mesh lines (2D) or faces (3D) called 
*			cells. The mesh line or face are also called cell faces. Grid points are
*			also called cell corners.
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef LEVEL_SET_H
#define LEVEL_SET_H

#include "mesh_with_data_packages.h"
#include "mesh_with_data_packages.hpp"
#include "geometry.h"

namespace SPH
{
	class LevelSet;
	class Kernel;
	/**
	 * @class LevelSetDataPackage
	 * @brief Fixed memory level set data packed in a package.
	 * Level set is the signed distance to an interface, 
	 * here, the surface of a body.
	 */
	class LevelSetDataPackage : public BaseDataPackage<4, 6>
	{
	public:
		bool is_core_pkg_;					 /**< If true, the package is near to zero level set. */
		PackageData<Real> phi_;				 /**< the level set or signed distance. */
		PackageDataAddress<Real> phi_addrs_; /**< address for the level set. */
		PackageData<Vecd> n_;				 /**< level set normalized gradient, to approximate interface normal direction */
		PackageDataAddress<Vecd> n_addrs_;
		PackageData<Real> kernel_weight_;
		PackageDataAddress<Real> kernel_weight_addrs_;
		PackageData<Vecd> kernel_gradient_;
		PackageDataAddress<Vecd> kernel_gradient_addrs_;
		/** mark the near interface cells. 0 for zero level set cut cells,
		  * -1 and 1 for negative and positive cut cells,  
		  * 0 can also be for other cells in the region closed 
		  * by negative and positive cut cells 
		  */
		PackageData<int> near_interface_id_;
		PackageDataAddress<int> near_interface_id_addrs_;

		LevelSetDataPackage();
		virtual ~LevelSetDataPackage(){};

		void assignAllPackageDataAddress(Vecu data_index, LevelSetDataPackage *src_pkg, Vecu addrs_index);
		void initializeBasicData(ComplexShape &complex_shape);
		void initializeWithUniformData(Real level_set);
		void computeKernelIntegrals(LevelSet &level_set);
		void computeNormalDirection();
		void stepReinitialization();
		void markNearInterface(Real small_shift_factor);
	};

	/**
	  * @class BaseLevelSet
	  * @brief A abstract describes a level set field defined on a mesh.
	  */
	class BaseLevelSet : public BaseMeshField
	{
	public:
		BaseLevelSet(ComplexShape &complex_shape, ParticleAdaptation &particle_adaptation);
		virtual ~BaseLevelSet(){};

		virtual bool probeIsWithinMeshBound(const Vecd &position) = 0;
		virtual Real probeSignedDistance(const Vecd &position) = 0;
		virtual Vecd probeNormalDirection(const Vecd &position) = 0;
		virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) = 0;
		virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) = 0;
		virtual void cleanInterface(bool isSmoothed = false) = 0;

	protected:
		ComplexShape &complex_shape_; /**< the geometry is described by the level set. */
		ParticleAdaptation &particle_adaptation_;

		/** for computing volume fraction occupied by a shape.*/
		Real computeHeaviside(Real phi, Real half_width);
	};

	/**
	 * @class LevelSet
	 * @brief Mesh with level set data as packages.
	 * Note that the mesh containing the data packages are cell-based 
	 * but within the data package, the data is grid-based.
	 */
	class LevelSet
		: public MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>
	{
	public:
		ConcurrentVector<LevelSetDataPackage *> core_data_pkgs_; /**< packages near to zero level set. */
		Real global_h_ratio_;
		Real small_shift_factor_;

		LevelSet(BoundingBox tentative_bounds, Real data_spacing,
				 ComplexShape &complex_shape, ParticleAdaptation &particle_adaptation, Real small_shift_factor = 0.75);
		virtual ~LevelSet(){};

		virtual bool probeIsWithinMeshBound(const Vecd &position) override;
		virtual Real probeSignedDistance(const Vecd &position) override;
		virtual Vecd probeNormalDirection(const Vecd &position) override;
		virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) override;
		virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) override;
		virtual void cleanInterface(bool isSmoothed = false) override;
		virtual void writeMeshFieldToPlt(std::ofstream &output_file) override;
		bool isWithinCorePackage(Vecd position);
		Real computeKernelIntegral(const Vecd &position);
		Vecd computeKernelGradientIntegral(const Vecd &position);

	protected:
		Kernel &kernel_;

		void reinitializeLevelSet();
		void markNearInterface();
		void redistanceInterface();
		void updateNormalDirection();
		void updateNormalDirectionForAPackage(LevelSetDataPackage *inner_data_pkg, Real dt = 0.0);
		void updateKernelIntegrals();
		void updateKernelIntegralsForAPackage(LevelSetDataPackage *inner_data_pkg, Real dt = 0.0);
		void stepReinitializationForAPackage(LevelSetDataPackage *inner_data_pkg, Real dt = 0.0);
		void markNearInterfaceForAPackage(LevelSetDataPackage *core_data_pkg, Real dt = 0.0);
		void redistanceInterfaceForAPackage(LevelSetDataPackage *core_data_pkg, Real dt = 0.0);
		virtual void initializeDataInACell(const Vecu &cell_index, Real dt) override;
		virtual void initializeAddressesInACell(const Vecu &cell_index, Real dt) override;
		virtual void tagACellIsInnerPackage(const Vecu &cell_index, Real dt) override;
		virtual void initializeDataPackages() override;
	};

	/**
	  * @class MultilevelCellLinkedList
	  * @brief Defining a multilevel level set for a complex region.
	  */
	class MultilevelLevelSet : public MultilevelMesh<BaseLevelSet, LevelSet>
	{
	public:
		MultilevelLevelSet(BoundingBox tentative_bounds, Real reference_data_spacing,
						   size_t total_levels, Real maximum_spacing_ratio,
						   ComplexShape &complex_shape, ParticleAdaptation &particle_adaptation);
		virtual ~MultilevelLevelSet(){};

		virtual bool probeIsWithinMeshBound(const Vecd &position) override;
		virtual Real probeSignedDistance(const Vecd &position) override;
		virtual Vecd probeNormalDirection(const Vecd &position) override;
		virtual Real probeKernelIntegral(const Vecd &position, Real h_ratio = 1.0) override;
		virtual Vecd probeKernelGradientIntegral(const Vecd &position, Real h_ratio = 1.0) override;
		virtual void cleanInterface(bool isSmoothed = false) override;

	protected:
		inline size_t getProbeLevel(const Vecd &position);
		inline size_t getMeshLevel(Real h_ratio);
	};
}
#endif //LEVEL_SET_H