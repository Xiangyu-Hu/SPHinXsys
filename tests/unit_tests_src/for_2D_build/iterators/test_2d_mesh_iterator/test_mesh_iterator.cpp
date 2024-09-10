#include "base_data_type.h"
#include "mesh_iterators.hpp"
#include <gtest/gtest.h>

using namespace SPH;

size_t get_cell_index(const Array2i &index)
{
    return 3 * (index[0] % 3) + (index[1] % 3);
}

size_t get_1d_index(const Array2i &all_cells, const Array2i &index)
{
    return all_cells[1] * index[0] + index[1];
}

void function(const Array2i &index, std::vector<size_t> &vec, const Array2i &all_cells)
{
    size_t index_cell = get_cell_index(index);
    size_t index_1d = get_1d_index(all_cells, index);
    vec[index_1d] += index_cell;
}

TEST(test_mesh_iterator, split_for)
{
    Array2i all_cells(9, 9);
    MeshRange mesh_range(Array2i::Zero(), all_cells);
    Array2i stride = 3 * Array2i::Ones();

    const size_t total_cells = all_cells[0] * all_cells[1];

    std::vector<size_t> vec_1(total_cells);
    std::vector<size_t> vec_2(total_cells);

    mesh_split_for(
        mesh_range, stride,
        [&](const Array2i &index)
        {
            function(index, vec_1, all_cells);
        });

    mesh_split_parallel_for(
        mesh_range, stride,
        [&](const Array2i &index)
        {
            function(index, vec_2, all_cells);
        });

    for (size_t i = 0; i < total_cells; ++i)
    {
        ASSERT_EQ(vec_1[i], vec_2[i]);
    }
}