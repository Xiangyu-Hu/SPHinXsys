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
* @file mesh_cell_linked_list.h
* @brief Here gives the classes for managing cell linked lists. This is the basic class 
* for building the particle configurations.
* @details  The cell linked list saves for each body a list of particles
* located within the cell.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_mesh.h"

/** this is a reformulation of tbb parallel_sort for particle data */
namespace tbb {
	namespace interafce9 {
		namespace internal {

			using tbb::internal::no_assign;

			/** sorting particle */
			template<typename RandomAccessIterator, typename Compare, typename SwapType>
			class QuickSortParticleRange : private no_assign {

				inline size_t median_of_three(const RandomAccessIterator& array, size_t l, size_t m, size_t r) const {
					return comp_(array[l], array[m]) ? (comp_(array[m], array[r]) ? m : (comp_(array[l], array[r]) ? r : l))
						: (comp_(array[r], array[m]) ? m : (comp_(array[r], array[l]) ? r : l));
				}

				inline size_t PseudoMedianOfNine(const RandomAccessIterator& array, const QuickSortParticleRange& range) const {
					size_t offset = range.size_ / 8u;
					return median_of_three(array,
						median_of_three(array, 0, offset, offset * 2),
						median_of_three(array, offset * 3, offset * 4, offset * 5),
						median_of_three(array, offset * 6, offset * 7, range.size_ - 1));

				}

				size_t splitRange(QuickSortParticleRange& range) {
					RandomAccessIterator array = range.begin_;
					RandomAccessIterator key0 = range.begin_;
					size_t m = PseudoMedianOfNine(array, range);
					if (m) swap_particle_data_(array, array + m);

					size_t i = 0;
					size_t j = range.size_;
					// Partition interval [i+1,j-1] with key *key0.
					for (;;) {
						__TBB_ASSERT(i < j, NULL);
						// Loop must terminate since array[l]==*key0.
						do {
							--j;
							__TBB_ASSERT(i <= j, "bad ordering relation?");
						} while (comp_(*key0, array[j]));
						do {
							__TBB_ASSERT(i <= j, NULL);
							if (i == j) goto quick_sort_particle_partition;
							++i;
						} while (comp_(array[i], *key0));
						if (i == j) goto quick_sort_particle_partition;
						swap_particle_data_(array + i, array + j);
					}
				quick_sort_particle_partition:
					// Put the partition key were it belongs
					swap_particle_data_(array + j, key0);
					// array[l..j) is less or equal to key.
					// array(j..r) is greater or equal to key.
					// array[j] is equal to key
					i = j + 1;
					size_t new_range_size = range.size_ - i;
					range.size_ = j;
					return new_range_size;
				}

			public:

				static const size_t grainsize_ = 500;
				const Compare& comp_;
				SwapType& swap_particle_data_;
				size_t size_;
				RandomAccessIterator begin_;

				QuickSortParticleRange(RandomAccessIterator begin,
					size_t size, const Compare& compare, SwapType& swap_particle_data) :
					comp_(compare), swap_particle_data_(swap_particle_data),
					size_(size), begin_(begin) {}

				bool empty() const { return size_ == 0; }
				bool is_divisible() const { return size_ >= grainsize_; }

				QuickSortParticleRange(QuickSortParticleRange& range, split)
					: comp_(range.comp_), swap_particle_data_(range.swap_particle_data_)
					, size_(splitRange(range))
					// +1 accounts for the pivot element, which is at its correct place
					// already and, therefore, is not included into subranges.
					, begin_(range.begin_ + range.size_ + 1) {}
			};

			/*
			Description : QuickSort in Iterator format
			Link        : https://stackoverflow.com/a/54976413/3547485
			Ref			: http://www.cs.fsu.edu/~lacher/courses/COP4531/lectures/sorts/slide09.html
			*/
			template <typename RandomAccessIterator, typename Compare, typename SwapType>
			RandomAccessIterator Partition(RandomAccessIterator first, RandomAccessIterator last, Compare& compare, SwapType& swap_particle_data)
			{
				auto pivot = std::prev(last, 1);
				auto i = first;
				for (auto j = first; j != pivot; ++j) {
					// bool format 
					if (compare(*j, *pivot)) {
						swap_particle_data(i++, j);
					}
				}
				swap_particle_data(i, pivot);
				return i;
			}

			template <typename RandomAccessIterator, typename Compare, typename SwapType>
			void SerialQuickSort(RandomAccessIterator first, RandomAccessIterator last, Compare& compare, SwapType& swap_particle_data)
			{
				if (std::distance(first, last) > 1) {
					RandomAccessIterator bound = Partition(first, last, compare, swap_particle_data);
					SerialQuickSort(first, bound, compare, swap_particle_data);
					SerialQuickSort(bound + 1, last, compare, swap_particle_data);
				}
			}

			/*
			Description : Insertsort in Iterator format
			Link        : http://www.codecodex.com/wiki/Insertion_sort
			*/

			template< typename RandomAccessIterator, typename Compare, typename SwapType>
			void InsertionSort(RandomAccessIterator First, RandomAccessIterator Last, Compare& compare, SwapType& swap_particle_data)
			{
				RandomAccessIterator min = First;
				for (RandomAccessIterator i = First + 1; i < Last; ++i)
					if (compare(*i, *min)) min = i;

				swap_particle_data(First, min);
				while (++First < Last)
					for (RandomAccessIterator j = First; compare(*j, *(j - 1)); --j)
						swap_particle_data((j - 1), j);
			}

			/** Body class used to sort elements in a range that is smaller than the grainsize. */
			template<typename RandomAccessIterator, typename Compare, typename SwapType>
			struct QuickSortParticleBody {
				void operator()(const QuickSortParticleRange<RandomAccessIterator, Compare, SwapType>& range) const {
					SerialQuickSort(range.begin_, range.begin_ + range.size_, range.comp_, range.swap_particle_data_);
				}
			};

		}
	}
}

namespace SPH {

	class SPHSystem;
	class SPHBody;
	class BaseParticles;
	class Kernel;

	/**
	 * @class CompareParticleSequence
	 * @brief compare the sequence of two particles
	 */
	struct CompareParticleSequence 
	{ 
		bool operator () (const size_t& x, const size_t& y) const 
		{ return x < y; }; 
	};

	/**
	 * @class SwapParticleData
	 * @brief swap sortable particle data according to a sequence
	 */
	class SwapParticleData
	{
	protected:
		StdLargeVec<size_t>& sequence_;
		StdLargeVec<size_t>& unsorted_id_;
		StdVec<StdLargeVec<Matd>*>& sortable_matrices_;
		StdVec<StdLargeVec<Vecd>*>& sortable_vectors_;
		StdVec<StdLargeVec<Real>*>& sortable_scalars_;

	public:
		SwapParticleData(BaseParticles* base_particles);
		~SwapParticleData() {};

		/** the operater overload for swapping particle data.
		 *  the arguments are the same with std::iter_swap
		 */
		void operator () (size_t* a, size_t* b);
	};

	/**
	 * @class CellList
	 * @brief The linked list for one cell
	 */
	class CellList
	{
	public:
		/** using concurrent vectors due to writting conflicts when building the list */
		ConcurrentIndexVector concurrent_particle_indexes_;
		/** non-concurrent cell linked list rewritten for building neighbor list */
		CellListDataVector cell_list_data_;
		/** the index vector for real particles. */
		IndexVector real_particle_indexes_;

		CellList();
		~CellList() {};
	};

	/**
	 * @class BaseMeshCellLinkedList
	 * @brief Abstract class for mesh cell linked list.
	 */
	class BaseMeshCellLinkedList : public Mesh
	{
	protected:
		SPHBody* body_;
		BaseParticles* base_particles_;
		Kernel* kernel_;
		SwapParticleData* swap_particle_data_;
		CompareParticleSequence compare_;
		tbb::interafce9::internal::
			QuickSortParticleRange<size_t*, CompareParticleSequence, SwapParticleData>* 
			quick_sort_particle_range_;
		tbb::interafce9::internal::
			QuickSortParticleBody<size_t*, CompareParticleSequence, SwapParticleData> 
			quick_sort_particle_body_;

		/** clear the cell lists */
		void ClearCellLists(Vecu& number_of_cells, matrix_cell cell_linked_lists);
		/** clear split cell lists in this mesh*/
		void ClearSplitCellLists(SplitCellLists& split_cell_lists);
		/** update split particle list in this mesh */
		void UpdateSplitCellLists(SplitCellLists& split_cell_lists,
			Vecu& number_of_cells, matrix_cell cell_linked_lists);
		/** update cell linked list data in this mesh */
		void UpdateCellListData(matrix_cell cell_linked_lists);
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		BaseMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound, 
			Real cell_spacing, size_t buffer_width = 2);
		/** Constructor with the direct information of the mesh. */
		BaseMeshCellLinkedList(SPHBody* body, 
			Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~BaseMeshCellLinkedList() {};

		/** computing search range for building contact configuration */
		int computeSearchRange(int origin_refinement_level, int target_refinement_level);
		/** choose a kernel for building up inter refinement level configuration */
		Kernel& ChoosingKernel(Kernel* original_kernel, Kernel* target_kernel);
		/** get the address of cell list */
		virtual CellList* CellListFromIndex(Vecu cell_index) = 0;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() = 0;

		/** Assign base particles to the mesh cell linked list. */
		void assignBaseParticles(BaseParticles* base_particles);
		/** Assign kernel to the mesh cell linked list. */
		void reassignKernel(Kernel* kernel);

		/** update the cell lists */
		virtual void UpdateCellLists() = 0;

		/** Insert a cell-linked_list entry. */
		virtual void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) = 0;
		virtual void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) = 0;

		/** find the nearest list data entry */
		virtual ListData findNearestListDataEntry(Vecd& position) = 0;

		/** sorting particle data according to the cell location of particles */
		virtual void sortingParticleData();
		/** computing the sequence which indicate the order of sorted particle data */
		virtual void computingSequence(StdLargeVec<size_t>& sequence) = 0;
		/** update the reference of sorted data from unsorted data */
		virtual void updateSortedId();

	};

	/**
	 * @class MeshCellLinkedList
	 * @brief Defining a mesh cell linked list for a body.
	 * The meshes for all bodies share the same global coordinates.
	 */
	class MeshCellLinkedList : public BaseMeshCellLinkedList
	{
	protected:
		/** cut_off radius */
		Real cutoff_radius_;
		/** The array for of mesh cells, i.e. mesh data.
		 * Within each cell, a list is saved with the indexes of particles.*/
		matrix_cell cell_linked_lists_;
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		MeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real cell_spacing, size_t buffer_width = 2);
		/** direct construct with mesh information. */
		MeshCellLinkedList(SPHBody* body, Vecd mesh_lower_bound,
			Vecu number_of_cells, Real cell_spacing);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~MeshCellLinkedList() { deleteMeshDataMatrix(); };

		/** access protected members */
		virtual CellList* CellListFromIndex(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() override { return cell_linked_lists_; };

		/** allcate memories for mesh data */
		virtual void allocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;

		/** output mesh data for visualization */
		virtual void writeMeshToVtuFile(ofstream &output_file) override {};
		virtual void writeMeshToPltFile(ofstream &output_file) override {};

		/** Insert a cell-linked_list entry. */
		void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override;

		/** find the nearest list data entry */
		virtual ListData findNearestListDataEntry(Vecd& position) override;
		/** computing the sequence which indicate the order of sorted particle data */
		virtual void computingSequence(StdLargeVec<size_t>& sequence) override;
	};

	/**
	  * @class MultilevelMeshCellLinkedList
	  * @brief Defining a multimesh cell linked list for a body
	  * for multiresolution particle configuration.
	  */
	class MultilevelMeshCellLinkedList : public BaseMeshCellLinkedList
	{
	protected:
		/** total levels of the storage pyramid.
		  * it is an odd value so that the reference level is the middle level.*/
		size_t total_levels_;
		/**cell spacing for the mesh levels*/
		StdVec<Real> cell_spacing_levels_;
		/** number of cells by dimension */
		StdVec<Vecu> number_of_cells_levels_;
		/** split cell list for building configuration.*/
		StdVec<SplitCellLists> split_cell_lists_levels_;
		/** Projected cell linked lists for building configuration.*/
		StdVec<MeshDataMatrix<CellList>> cell_linked_lists_levels_;
		/** point to every mesh level. */
		StdVec<MeshCellLinkedList*> mesh_cell_linked_list_levels_;

		/** determine mesh level of a particle. */
		size_t getLevelFromCutOffRadius(Real smoothing_length);
		/** find cell indexes from point position in a level */
		Vecu getLevelCellIndexesFromPosition(Vecd& position, size_t level);
		/** Insert a cell-linked_list entry to the projected particle list. */
		void InsertACellLinkedListEntryAtALevel(size_t particle_index, Vecd& position, Vecu& cell_index, size_t level);
	public:
		/** Constructor to achieve alignment of all mesh levels. */
		MultilevelMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real reference_cell_spacing, size_t total_levels = 1, size_t buffer_width = 2);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~MultilevelMeshCellLinkedList() {};

		/** access protected members */
		virtual CellList* CellListFromIndex(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() override { return cell_linked_lists_levels_[0]; };

		/** allcate memories for mesh data */
		virtual void allocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;

		/** Insert a cell-linked_list entry to the projected particle list. */
		void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override {};

		/** find the nearest list data entry */
		virtual ListData findNearestListDataEntry(Vecd& position) override { return ListData(0, Vecd(0)); };
		/** computing the sequence which indicate the order of sorted particle data */
		virtual void computingSequence(StdLargeVec<size_t>& sequence) override {};
	};
}
