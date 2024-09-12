#include "base_data_type.h"
#include "mesh_iterators.hpp"
#include <gtest/gtest.h>

using namespace SPH;

size_t get_cell_index(const Array3i &index)
{
    return 9 * (index[0] % 3) + 3 * (index[1] % 3) + (index[2] % 3);
}

size_t get_1d_index(const Array3i &all_cells, const Array3i &index)
{
    return all_cells[1] * all_cells[2] * index[0] + all_cells[2] * index[1] + index[2];
}

void function(const Array3i &index, std::vector<size_t> &vec, const Array3i &all_cells)
{
    size_t index_cell = get_cell_index(index);
    size_t index_1d = get_1d_index(all_cells, index);
    vec[index_1d] += index_cell;
}

TEST(test_mesh_iterator, split_for)
{
    Array3i all_cells(10, 10, 10);
    MeshRange mesh_range(Array3i::Zero(), all_cells);
    Array3i stride = 3 * Array3i::Ones();

    const size_t total_cells = all_cells[0] * all_cells[1] * all_cells[2];

    std::vector<size_t> vec_1(total_cells);
    std::vector<size_t> vec_2(total_cells);

    // forward test
    mesh_stride_forward_for(
        mesh_range, stride,
        [&](const Array3i &index)
        {
            function(index, vec_2, all_cells);
        });

    mesh_stride_forward_parallel_for(
        mesh_range, stride,
        [&](const Array3i &index)
        {
            function(index, vec_1, all_cells);
        });

    for (size_t i = 0; i < total_cells; ++i)
        ASSERT_EQ(vec_1[i], vec_2[i]);

    // backward test
    mesh_stride_backward_for(
        mesh_range, stride,
        [&](const Array3i &index)
        {
            function(index, vec_2, all_cells);
        });

    mesh_stride_backward_parallel_for(
        mesh_range, stride,
        [&](const Array3i &index)
        {
            function(index, vec_1, all_cells);
        });

    for (size_t i = 0; i < total_cells; ++i)
        ASSERT_EQ(vec_1[i], vec_2[i]);
}