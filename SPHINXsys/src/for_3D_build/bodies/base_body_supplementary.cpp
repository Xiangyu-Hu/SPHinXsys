/**
 * @file 	base_body_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_body.h"
#include "base_particles.h"
#include "mesh_cell_linked_list.h"

namespace SPH 
{
	//=================================================================================================//
	void BodyPartByCell::TagBodyPart()
	{
		BaseMeshCellLinkedList *mesh_cell_linked_list
			= body_->base_mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->NumberOfCells();

		for (int i = 0; i < number_of_cells[0]; ++i)
			for (int j = 0; j < number_of_cells[1]; ++j)
				for (int k = 0; k < number_of_cells[2]; ++k)
				{
					bool is_contained = false;
					for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++m)
							for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells[2]) - 1); ++n)
							{
								Vecd cell_position = mesh_cell_linked_list
									->CellPositionFromIndexes(Vecu(l, m, n));
								if (body_part_shape_.checkContain(cell_position))
									is_contained = true;
							}
					if (is_contained == true) 
						body_part_cells_.push_back(mesh_cell_linked_list->CellListFormIndex(Vecu(i, j, k)));
				}
	}
	//=================================================================================================//
	void NearBodySurface::TagBodyPart()
	{
		BaseMeshCellLinkedList* mesh_cell_linked_list
			= body_->base_mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->NumberOfCells();

		for (int i = 0; i < number_of_cells[0]; ++i)
			for (int j = 0; j < number_of_cells[1]; ++j)
				for (int k = 0; k < number_of_cells[2]; ++k)
				{
					bool is_near = false;
					for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++m)
							for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells[2]) - 1); ++n)
							{
								Vecd cell_position = mesh_cell_linked_list->CellPositionFromIndexes(Vecu(l, m, n));
								if (body_->levelset_mesh_->isWithinMeshBound(cell_position)) 
								{
									Real phii =body_->levelset_mesh_->probeLevelSet(cell_position);
									if (fabs(phii) <= mesh_cell_linked_list->GridSpacing()) is_near = true;
								}
							}
					if (is_near == true)
						body_part_cells_.push_back(mesh_cell_linked_list->CellListFormIndex(Vecu(i, j, k)));
				}
	}
	//=================================================================================================//
	void SolidBodyPartForSimbody::TagBodyPart()
	{
		BodyPartByParticle::TagBodyPart();

		Real body_part_volume(0);
		initial_mass_center_ = Vec3d(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];

			initial_mass_center_ += base_particle_data_i.Vol_ * base_particle_data_i.pos_0_;
			body_part_volume += base_particle_data_i.Vol_;
		}

		initial_mass_center_ /= body_part_volume;

		//computing unit intertia
		Vec3d intertia_moments(0);
		Vec3d intertia_products(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];

			Vec3d displacement = (base_particle_data_i.pos_0_ - initial_mass_center_);
			intertia_moments[0] += base_particle_data_i.Vol_
				* (displacement[1] * displacement[1] + displacement[2] * displacement[2]);
			intertia_moments[1] += base_particle_data_i.Vol_
				* (displacement[0] * displacement[0] + displacement[2] * displacement[2]);
			intertia_moments[2] += base_particle_data_i.Vol_
				* (displacement[0] * displacement[0] + displacement[1] * displacement[1]);
			intertia_products[0] -= base_particle_data_i.Vol_ * displacement[0] * displacement[1];
			intertia_products[1] -= base_particle_data_i.Vol_ * displacement[0] * displacement[2];
			intertia_products[2] -= base_particle_data_i.Vol_ * displacement[1] * displacement[2];

		}
		intertia_moments /= body_part_volume;
		intertia_products /= body_part_volume;

		body_part_mass_properties_
			= new SimTK::MassProperties(body_part_volume * solid_body_density_,
				Vec3d(0), SimTK::UnitInertia(intertia_moments, intertia_products));
	}
	//=================================================================================================//
}
