
#include "outflow_boundary.h"

namespace SPH
{
    void ModifiedDoNothingConditionWithFace::UpdateOutletRegionParticles(size_t index_i, Real dt)
    {
        const Vecd& postion = pos_n_[index_i];
        Vecu target_cell_index = cell_linked_list_->CellIndexFromPosition(postion);
        int ii = (int)target_cell_index[0];
        int jj = (int)target_cell_index[1];
        int kk = (int)target_cell_index[2];

        Neighborhood neighborhood;
        for (int l = SMAX(ii - 1, 0); l <= SMIN(ii + 1, int(cell_linked_list_->NumberOfCells()[0]) - 1); ++l)
            for (int m = SMAX(jj - 1, 0); m <= SMIN(jj + 1, int(cell_linked_list_->NumberOfCells()[1]) - 1); ++m)
                for (int n = SMAX(kk - 1, 0); n <= SMIN(jj + 1, int(cell_linked_list_->NumberOfCells()[2] - 1)); ++n)
                {
                    ListDataVector& target_particles = cell_linked_list_->getCellLists()[l][m][n].cell_list_data_;
                    for (const ListData& list_data : target_particles)
                    {
                        if (body_region_.getSignedDistance(list_data.second) < body_region_.getRegionWidth()) // outlet region particles
                            continue;

                        // displacement pointing from neighboring particle to origin particle
                        Vecd displacement = postion - list_data.second;
                        relation_inner_(neighborhood, displacement, index_i, list_data.first);
                    }
                }

        Vecd sum_v(0);
        Real denominator = 0.0;

        for (size_t n = 0; n < neighborhood.current_size_; ++n)
        {
            size_t index_j = neighborhood.j_[n];
            Real W_ij_V_j = neighborhood.W_ij_[n] * particles_->Vol_[index_j];
            sum_v += W_ij_V_j * vel_n_[index_j];
            denominator += W_ij_V_j;
        }

        if (denominator > Eps)
        {
            vel_n_[index_i] = sum_v / denominator;
        }
    }

} // namespace SPH