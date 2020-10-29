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
* @version	0.1
* @version  0.2.0
* 			Now narrow bounded levelset mesh is added to replace the whole domain background levelset mesh. 
*/

#pragma once

#include "mesh_with_data_packages.h"
#include "mesh_with_data_packages.hpp"
#include "geometry.h"

namespace SPH
{
	/**
	 * @class LevelSetDataPackage
	 * @brief Fixed memory level set data packed in a package.
	 */
	class LevelSetDataPackage : public BaseDataPackage<4, 6>
	{
	public:
		bool is_core_pkg_;	/**< If true, the package is near to zero level set. */

		/** level set is the signed distance to an interface, 
		  * here, the surface of a body */
		PackageData<Real> phi_;
		PackageDataAddress<Real> phi_addrs_;
		/** level set normalized gradient, to approximate interface normal direction */
		PackageData<Vecd> n_;
		PackageDataAddress<Vecd> n_addrs_;
		/** level set curvature, to approximate interface curvature */
		PackageData<Real> kappa_;
		PackageDataAddress<Real> kappa_addrs_;
		/** mark the near interface cells. 0 for zero level set cut cells,
		  * -1 and 1 for negative and positive cut cells,  
		  * 0 can also be for other cells in the region closed 
		  * by negative and positive cut cells 
		  */
		PackageData<int> near_interface_id_;
		PackageDataAddress<int> near_interface_id_addrs_;

		/** default constructor */
		LevelSetDataPackage();
		virtual ~LevelSetDataPackage() {};
	
		/** This function assigns a set of package data address for all variables. */
		void assignAllPackageDataAddress(Vecu data_index, LevelSetDataPackage* src_pkg, Vecu addrs_index);
		/** This function initialize the level set package constructed by default. */
		void initializeDataPackage(ComplexShape& complex_shape);
		/**
		 *@brief This function initialize with uniform level set field
		 *@param[in] level_set(Real) Level set value
		 *@param[in] normal_direction(Vecd) Normal direction of levelset field.
		 */
		void initializeWithUniformData(Real level_set, Vecd normal_direction);
		/** This function compute normal direction for all level set in the package */
		void computeNormalDirection();
		/** This function applies one step reinitialization */
		void stepReinitialization();
		/** This function marks the near interface ids */
		void markNearInterface();
	};

	/**
	  * @class BaseLevelSet
	  * @brief A abstract describes a mesh with level set data packages.
	  */
	class BaseLevelSet : public BaseMeshWithDataPackages
	{
	public:
		/** Constructor using domain information. */
		BaseLevelSet(Vecd lower_bound, 		/**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spcaing. */
			size_t buffer_width = 0	/**< BUffer size. */
		);
		/** Constructor using mesh information directly. */
		BaseLevelSet(Vecd mesh_lower_bound, /**< Lower bound. */
			Vecu number_of_cells,  /**< Upper bound. */
			Real cell_spacing		/**< Cell spcaing. */
		);
		virtual ~BaseLevelSet() {};

		/**
		 *@brief This function probe the level set at a off-grid position.
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Real probeLevelSet(Vecd position) = 0;
		/**
		 *@brief This function probe the normal direction at a off-grid position
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Vecd probeNormalDirection(Vecd position) = 0;
		/**update the normal direction */
		virtual void updateNormalDirection() = 0;
		/**
		*@brief This function: 
		* 1) mark the cutcell in whole domain include zero/positive/negative cutcell;
		* 2) clean the non-consistency 0 levelset in whole domain;
		* 3) reinitialize the levelset vaule in whole domain after clean opeartion;
		* 4) smoothing the levelset according to curvature (optional);
		*/
		virtual void cleanInterface(bool isSmoothed = false) = 0;
		/**
		 *@brief This function reinitialize levelset value in a whole domain
		 */
		virtual void reinitializeLevelSet() = 0;
		/**
		 *@brief This function marks the near interface
		 */
		virtual void markNearInterface() = 0;
		/**
		 *@brief This function delete unresolved interface by redistance
		 */
		virtual void redistanceInterface() = 0;
		/**
		*@brief This function calculate the integration of kernel function outside the surafce
		*/
		virtual Vecd computeKernelIntegral(Vecd position, Kernel* kernel) = 0;
	protected:
		Real computeHeaviside(Real phi, Real half_width);
	};

	/**
	 * @class LevelSet
	 * @brief Mesh with level set data as packages.
	 */
	class LevelSet
		: public MeshWithDataPackages<BaseLevelSet, LevelSetDataPackage>
	{
	public:
		/** Core packages which are near to zero level set. */
		ConcurrentVector<LevelSetDataPackage*> core_data_pkgs_;

		/** Constructor using domain and sph body information. */
		LevelSet(ComplexShape& complex_shape,  	/**< Link to geomentry. */
			Vecd lower_bound,      /**< Lower bound. */
			Vecd upper_bound, 		/**< Upper bound. */
			Real grid_spacing, 	/**< Grid spcaing. */
			size_t buffer_width = 0 /**< Buffer size. */
		);
		virtual ~LevelSet() {};

		/**
		 *@brief This function initialize the Levelset data package.
		 */
		virtual void initializeDataPackages() override;
		/**
		 *@brief This function probe the level set at a off-grid position.
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Real probeLevelSet(Vecd position) override;
		/**
		 *@brief This function probe the normal direction at a off-grid position
		 *@param[in] position(Vecd) The enquiry postion
		 */
		virtual Vecd probeNormalDirection(Vecd position) override;
		/**
		 *@brief This function update the norm of levelset field using central difference scheme.
		 */
		virtual void updateNormalDirection() override;
		/**
		 *@brief This function clean the zero level set.
		 */
		virtual void cleanInterface(bool isSmoothed = false) override;
		/**
		 *@brief This function reinitialize levelset value in a whole domain
		 */
		virtual void reinitializeLevelSet() override;
		/**
		 *@brief This function marks the near interface
		 */
		virtual void markNearInterface() override;
		/**
		 *@brief This function redistance unresolved interface
		 */
		virtual void redistanceInterface() override;
		/**
		 *@brief This function output mesh data for Paraview
		 *@param[out] output_file(ofstream) output ofstream.
		 */
		virtual void writeMeshToVtuFile(ofstream& output_file) override {};
		/**
		 *@brief This function output mesh data for Tecplot visualization
		 *@param[out] output_file(ofstream) output ofstream.
		 */
		virtual void writeMeshToPltFile(ofstream& output_file) override;
		/*test below*/
		/**
		*@brief This function calculate the integration of kernel function outside the surafce
		*/
		virtual Vecd computeKernelIntegral(Vecd position, Kernel* kernel) override;

	protected:
		/**the geometry is described by the level set. */
		ComplexShape& complex_shape_;
		/**
		 *@brief This function initialize level set in a cell.
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void initializeDataInACell(Vecu cell_index, Real dt) override;
		/**
		 *@brief This function initialize the addresses in a data package
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void initializeAddressesInACell(Vecu cell_index, Real dt) override;
		/**
		 *@brief This function tag if a data package is inner package.
		 *@param[in] cell_index(Vecu) Index of cell
		 *@param[in] dt(Real) not used herein.
		 */
		virtual void tagACellIsInnerPackage(Vecu cell_index, Real dt) override;
		void updateNormalDirectionForAPackage(LevelSetDataPackage* inner_data_pkg, Real dt = 0.0);
		void stepReinitializationForAPackage(LevelSetDataPackage* inner_data_pkg, Real dt = 0.0);
		void markNearInterfaceForAPackage(LevelSetDataPackage* core_data_pkg, Real dt = 0.0);
		void redistanceInterfaceForAPackage(LevelSetDataPackage* core_data_pkg, Real dt = 0.0);

	};
}
