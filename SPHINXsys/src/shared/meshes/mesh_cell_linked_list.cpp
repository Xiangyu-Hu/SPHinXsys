/**
 * @file 	mesh_cell_linked_list.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "mesh_cell_linked_list.h"
#include "base_kernel.h"
#include "base_body.h"
#include "base_particles.h"


namespace SPH {
	//=================================================================================================//
	BaseMeshCellLinkedList
		::BaseMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real cell_spacing, size_t buffer_size)
		: Mesh(lower_bound, upper_bound, cell_spacing, buffer_size), 
		body_(body), contact_map_(),
		base_particles_(NULL), kernel_(body->kernel_) {}
	//=================================================================================================//
	BaseMeshCellLinkedList
		::BaseMeshCellLinkedList(SPHBody* body, 
			Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing)
		: Mesh(mesh_lower_bound, number_of_cells, cell_spacing),
		body_(body), contact_map_(),
		base_particles_(NULL), kernel_(body->kernel_) {}
	//=================================================================================================//
	int BaseMeshCellLinkedList::ComputingSearchRage(int orign_refinement_level,
		int target_refinement_level)
	{
		return orign_refinement_level >= target_refinement_level
			? 1 : powern(2, target_refinement_level - orign_refinement_level);
	}
	//=================================================================================================//
	Kernel& BaseMeshCellLinkedList
		::ChoosingKernel(Kernel* orignal_kernel, Kernel* target_kernel)
	{
		return orignal_kernel->GetSmoothingLength()
				> target_kernel->GetSmoothingLength()
			? *orignal_kernel : *target_kernel;
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::assignParticles(BaseParticles* base_particles) 
	{ 
		base_particles_ = base_particles; 
	};
	//=================================================================================================//
	void BaseMeshCellLinkedList::assignContactMap(SPHBodyContactMap contact_map)
	{
		contact_map_ = contact_map;
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::reassignKernel(Kernel* kernel)
	{
		kernel_ = kernel;
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::ClearSplitCellLists(SplitCellLists& split_cell_lists)
	{
		for (size_t i = 0; i < split_cell_lists.size(); i++)
			split_cell_lists[i].clear();
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList
		::BuildInnerConfiguration(ParticleConfiguration& inner_configuration)
	{
		UpdateInnerConfiguration(inner_configuration);
	}
	void BaseMeshCellLinkedList::UpdateContactConfiguration()
	{
		UpdateInteractionConfiguration(contact_map_.second);
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::BuildContactConfiguration()
	{
		UpdateContactConfiguration();
	}
	//=================================================================================================//
	MeshCellLinkedList::MeshCellLinkedList(SPHBody* body, Vecd lower_bound,
		Vecd upper_bound, Real cell_spacing, size_t buffer_size)
		: BaseMeshCellLinkedList(body, lower_bound, upper_bound, cell_spacing, buffer_size),
		cutoff_radius_(cell_spacing) {}
	//=================================================================================================//
	MeshCellLinkedList::MeshCellLinkedList(SPHBody* body, Vecd mesh_lower_bound,
		Vecu number_of_cells, Real cell_spacing)
		: BaseMeshCellLinkedList(body, mesh_lower_bound, number_of_cells, cell_spacing),
		cutoff_radius_(cell_spacing) {}
	//=================================================================================================//
	void MeshCellLinkedList::UpdateCellLists()
	{
		ClearCellLists(number_of_cells_, cell_linked_lists_);
		StdLargeVec<BaseParticleData>& base_particle_data = base_particles_->base_particle_data_;
		size_t number_of_particles = body_->number_of_particles_;
		//rebuild the corresponding particle list.
		parallel_for(blocked_range<size_t>(0, number_of_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					InsertACellLinkedListEntry(i, base_particle_data[i].pos_n_);
				}
			}, ap);
		UpdateSplitCellLists(body_->split_cell_lists_, number_of_cells_, cell_linked_lists_);
	}
	//=================================================================================================//
	MultilevelMeshCellLinkedList
		::MultilevelMeshCellLinkedList(SPHBody* body, Vecd lower_bound,
		Vecd upper_bound, Real reference_cell_spacing, size_t total_levels, size_t buffer_size)
		: BaseMeshCellLinkedList(body, lower_bound, upper_bound, reference_cell_spacing, buffer_size),
		total_levels_(total_levels)
	{
		/**intialize lists. */
		split_cell_lists_levels_.resize(total_levels);
		cell_linked_lists_levels_.resize(total_levels);

		/** bulid the zero level mesh first.*/
		size_t middle_level = (total_levels - 1) / 2;
		Real zero_level_cell_spacing = reference_cell_spacing * powern(2.0, middle_level);
		cell_spacing_levels_.push_back(zero_level_cell_spacing);
		MeshCellLinkedList* zero_level_mesh
			= new MeshCellLinkedList(body, lower_bound,	upper_bound, zero_level_cell_spacing, buffer_size);
		mesh_cell_linked_list_levels_.push_back(zero_level_mesh);
		size_t number_of_split_cell_lists = powern(3, Vecd(0).size());
		split_cell_lists_levels_[0].resize(number_of_split_cell_lists);
		Vecu zero_level_number_of_cells = zero_level_mesh->getNumberOfCells();
		number_of_cells_levels_.push_back(zero_level_number_of_cells);
		/** copy zero evel mesh perperties to this. */
		copyMeshProperties(zero_level_mesh);

		/** other levels. */
		for (size_t level = 1; level != total_levels; ++level) {
			Real cell_spacing = zero_level_cell_spacing * powern(0.5, level);
			cell_spacing_levels_.push_back(cell_spacing);
			Vecu number_of_cells
				= zero_level_number_of_cells * powern(2, level);
			MeshCellLinkedList* mesh_cell_linked_list_level
				= new MeshCellLinkedList(body, mesh_lower_bound_, number_of_cells, cell_spacing);
			mesh_cell_linked_list_levels_.push_back(mesh_cell_linked_list_level);
			split_cell_lists_levels_[level].resize(number_of_split_cell_lists);
			number_of_cells_levels_.push_back(number_of_cells);
		}
	}
	//=================================================================================================//
	size_t MultilevelMeshCellLinkedList
		::getLevelFromCutOffRadius(Real smoothing_length)
	{
		size_t current_level = 0;
		Real cut_off_radius = kernel_->GetCutOffRadius(smoothing_length);
		for (size_t level = 1; level != cell_spacing_levels_.size(); ++level)
		{
			if (cut_off_radius < 0.5 * cell_spacing_levels_[level - 1]) current_level = level;
		}
		return current_level;
	}
	//=================================================================================================//
	Vecu MultilevelMeshCellLinkedList::getLevelCellIndexesFromPosition(Vecd& position, size_t level)
	{
		Real cell_spacing = cell_spacing_levels_[level];
		Vecu& number_of_cells = number_of_cells_levels_[level];
		Vecd rltpos = position - mesh_lower_bound_;
		Vecu cell_pos(0);
		for (int n = 0; n < rltpos.size(); n++)
		{
			cell_pos[n] = clamp((int)floor(rltpos[n] / cell_spacing),
				0, int(number_of_cells[n]) - 1);
		}
		return cell_pos;
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::
		InsertACellLinkedListEntry(size_t index_i, Vecd particle_position)
	{
		BaseParticleData& base_particle_data_i
			= base_particles_->base_particle_data_[index_i];
		size_t current_level = 0;
		Real cut_off_radius = kernel_->GetCutOffRadius(base_particle_data_i.smoothing_length_);
		Vecd& position = base_particle_data_i.pos_n_;
		Vecu curretn_cell_index(0);
		for (size_t level = 1; level != cell_spacing_levels_.size(); ++level)
		{
			if (cut_off_radius < 0.5 * cell_spacing_levels_[level - 1]) {
				current_level = level;
				curretn_cell_index = getLevelCellIndexesFromPosition(position, level);
				InsertACellLinkedListEntryAtALevel(index_i, position, curretn_cell_index, level);
			}
		}
		mesh_cell_linked_list_levels_[current_level]->InsertACellLinkedListEntry(index_i, particle_position);
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::AllocateMeshDataMatrix()
	{
		for (size_t l = 0; l != total_levels_; ++l) {
			mesh_cell_linked_list_levels_[l]->AllocateMeshDataMatrix();
		}
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::DeleteMeshDataMatrix()
	{
		for (size_t l = 0; l != total_levels_; ++l) {
			mesh_cell_linked_list_levels_[l]->DeleteMeshDataMatrix();
		}
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::UpdateCellLists()
	{
		for (size_t level = 0; level != total_levels_; ++level) {
			matrix_cell cell_linked_list
				= mesh_cell_linked_list_levels_[level]->getCellLinkedLists();
			ClearCellLists(number_of_cells_levels_[level], cell_linked_list);
			ClearCellLists(number_of_cells_levels_[level], cell_linked_lists_levels_[level]);

		}
		StdLargeVec<BaseParticleData>& base_particle_data 
			= base_particles_->base_particle_data_;
		size_t number_of_particles = body_->number_of_particles_;
		//rebuild the corresponding particle list.
		parallel_for(blocked_range<size_t>(0, number_of_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					InsertACellLinkedListEntry(i, base_particle_data[i].pos_n_);
				}
			}, ap);

		for (size_t level = 0; level != total_levels_; ++level) {
			matrix_cell cell_linked_list
				= mesh_cell_linked_list_levels_[level]->getCellLinkedLists();
			UpdateSplitCellLists(split_cell_lists_levels_[level],
				number_of_cells_levels_[level], cell_linked_list);
		}
		UpdateSplitCellLists(body_->split_cell_lists_, 
			number_of_cells_levels_[0], cell_linked_lists_levels_[0]);
	}
	//=================================================================================================//
}