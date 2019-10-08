/**
* @file 	base_mesh.h
* @brief 	This is the base classes of mesh, which describe ordered and indexed
*			data sets.  Depending on application, there are different data 
* 			saved on the mesh. The intersection points of mesh lines are called 
*			grid points, the element enclosed by mesh lines (2D) or faces (3D) called 
*			cells. The mesh line or face are also called cell faces. Grid points are
*			also called cell corners.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_data_package.h"

#include <fstream>

using namespace std;

namespace SPH {

	/**
	 * @class Mesh
	 * @brief Abstract base class for defining basic mesh prpoerties.
	 * The mesh is proposed for several functions.
	 * First, it is used in cell linked list for neighbor search.
	 * Second, it is used for background maps such as level sets.
	 * This class is the counterpart of the class particles.
	 */
	class Mesh
	{

	protected:
		/** buffer to avoid bound check */
		size_t buffer_size_;
		/** bounds */
		Vecd mesh_lower_bound_, mesh_upper_bound_;
		/** grid spcing */
		Real grid_spacing_;
		/** number grid points by dimension */
		Vecu number_of_grid_points_;

		/** computing number of total lattices */
		Vecu CalcNumberOfGridPoints();
	public:
		/** Constructor */
		Mesh(Vecd lower_bound, Vecd upper_bound, Real grid_spacing, size_t buffer_size = 0);
		virtual ~Mesh() {};

		/** accesss protected variables */
		Vecd GetLowerBound() { return mesh_lower_bound_; };
		Vecd GetUpperBound() { return mesh_upper_bound_; };
		Real GetGridSpacing() { return grid_spacing_; };
		Vecu GetNumberOfGridPoints() { return number_of_grid_points_; };

		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() = 0;

		/** find grid indexes from point poistion */
		Vecu GridIndexesFromPosition(Vecd position);
		/** find grid poistion from indexes */
		Vecd GridPositionFromIndexes(Vecu gird_indexes);
		/** Find cell position from indexes.
		  * It is the position shift in the upper-right direction half grid size */
		Vecd CellPositionFromIndexes(Vecu cell_indexes);

		/** output mesh data for Paraview visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) = 0;
		/** output mesh data for Tecplot visuallization */
		virtual void WriteMeshToPltFile(ofstream &output_file) = 0;

	};

	/**
	 * @class BackgroundData
	 * @brief Data for background mesh.
	 * Describing complex geometrics using level set, 
	 * which is the distance to the surface of the geometry 
	 * and the direction leads to the nearest point on the surface    
	 */
	class BackgroundData
	{
	public:
		/** Empty constructor */
		BackgroundData() {};
		BackgroundData(Real level_set, Vecd normal_direction);
		virtual ~BackgroundData() {};

		/** level set is the signed distance to
		  * an interface, here, the surface of a body */
		Real phi_;
		/** displacment to the surface */
		Vecd n_;
		/** curvature */
		Real kappa_;
	};

	/**
	 * @class MeshBackground
	 * @brief Background mesh for inital particle relaxation.
	 */
	class MeshBackground : public Mesh
	{
		/** base mesh data */
		matrix_grid mesh_background_data_;

	public:
		MeshBackground(Vecd lower_bound, Vecd upper_bound, 
			Real grid_spacing, size_t buffer_size = 0);
		virtual ~MeshBackground() { DeleteMeshDataMatrix(); };

		/** allocate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() override;

		/** initialize level set and displacement to surface
		  * for body region geometry */
		void InitializeLevelSetData(SPHBody &body);
		void ComputeCurvatureFromLevelSet(SPHBody &body);
		/** probe the base mesh data */
		Vecd ProbeNormalDirection(Vecd Point);
		Real ProbeLevelSet(Vecd Point);
		Real ProbeCurvature(Vecd Point);

		/** output mesh data for visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) override;
		virtual void WriteMeshToPltFile(ofstream &output_file) override;
	};
}
