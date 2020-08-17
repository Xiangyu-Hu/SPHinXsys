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
			for (int j = 0; j < number_of_cells[1]; ++j) {
				bool is_contained = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++l)
					{
						Vecd cell_position = mesh_cell_linked_list
							->CellPositionFromIndexes(Vecu(k, l));
						if (body_part_shape_.checkContain(cell_position))
							is_contained = true;
					}
				if (is_contained == true) 
					body_part_cells_.push_back(mesh_cell_linked_list->CellListFormIndex(Vecu(i, j)));
			}
	}
	//=================================================================================================//
	void NearBodySurface::TagBodyPart()
	{
		BaseMeshCellLinkedList* mesh_cell_linked_list
			= body_->base_mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->NumberOfCells();

		for (int i = 0; i < number_of_cells[0]; ++i)
			for (int j = 0; j < number_of_cells[1]; ++j) {

				bool is_near = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++l)
					{
						Vecd cell_position = mesh_cell_linked_list->CellPositionFromIndexes(Vecu(k, l));
						if (body_->levelset_mesh_->isWithinMeshBound(cell_position)) 
						{
							Real phii = body_->levelset_mesh_->probeLevelSet(cell_position);
							if (fabs(phii) <= mesh_cell_linked_list->GridSpacing()) is_near = true;
						}
					}
				if (is_near == true)
					body_part_cells_.push_back(mesh_cell_linked_list->CellListFormIndex(Vecu(i, j)));
			}
	}
	//=================================================================================================//
	void SolidBodyPartForSimbody::TagBodyPart()
	{
		BodyPartByParticle::TagBodyPart();

		Real body_part_volume(0);
		Vecd mass_center = Vecd(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];

			mass_center += base_particle_data_i.Vol_ * base_particle_data_i.pos_0_;
			body_part_volume += base_particle_data_i.Vol_;
		}

		mass_center /= body_part_volume;
		initial_mass_center_ = Vec3d(mass_center[0], mass_center[1], 0.0);

		//computing unit intertia
		Real Ix = 0.0;
		Real Iy = 0.0;
		Real Iz = 0.0;
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];

			Vecd displacement = (base_particle_data_i.pos_0_ - mass_center);
			Real r_x = (base_particle_data_i.pos_0_[1] - mass_center[1]);
			Ix += base_particle_data_i.Vol_ * r_x * r_x;
			Real r_y = (base_particle_data_i.pos_0_[0] - mass_center[0]);
			Iy += base_particle_data_i.Vol_ * r_y * r_y;
			Iz += base_particle_data_i.Vol_
				* (base_particle_data_i.pos_0_ - mass_center).normSqr();
		}
		Ix /= body_part_volume;
		Iy /= body_part_volume;
		Iz /= body_part_volume;

		body_part_mass_properties_
			= new SimTK::MassProperties(body_part_volume * solid_body_density_,
				Vec3d(0), SimTK::UnitInertia(Ix, Iy, Iz));
	}
	//=================================================================================================//
}
