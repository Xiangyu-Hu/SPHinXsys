/**
 * @file 	particle_iterators.cpp
 * @brief 	This is the implementation of the template class particle particle iterators
 * @author	Xiangyu Hu
 */

#include "particle_iterators.h"

#include "cell_linked_list.h"

//=================================================================================================//
namespace SPH
{
    //=================================================================================================//
    size_t SizeOfLoopRange(const size_t &all_real_particles)
    {
        return all_real_particles;
    };
    //=============================================================================================//
    size_t SizeOfLoopRange(const IndexVector &body_part_particles)
    {
        return body_part_particles.size();
    };
    //=============================================================================================//
    size_t SizeOfLoopRange(const CellLists &body_part_cells)
    {
        size_t size_of_loop_range = 0;
        for (size_t i = 0; i != body_part_cells.size(); ++i)
        {
            size_of_loop_range += body_part_cells[i]->cell_list_data_.size();
        }
        return body_part_cells.size();
    };
    //=============================================================================================//
}
//=================================================================================================//