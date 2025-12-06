#include "algorithm_primitive.h"

namespace SPH
{
//=================================================================================================//
QuickSort::SwapIndex::SwapIndex(UnsignedInt *sequence, UnsignedInt *index_permutation)
    : sequence_(sequence), index_permutation_(index_permutation) {}
//=================================================================================================//
void QuickSort::SwapIndex::operator()(UnsignedInt *a, UnsignedInt *b)
{
    std::swap(*a, *b);

    UnsignedInt index_a = a - sequence_;
    UnsignedInt index_b = b - sequence_;
    std::swap(index_permutation_[index_a], index_permutation_[index_b]);
}
//=================================================================================================//
void QuickSort::sort(const ParallelPolicy &ex_policy, UnsignedInt size, UnsignedInt start_index)
{
    quick_sort_range_.begin_ = sequence_ + start_index;
    quick_sort_range_.size_ = size;
    tbb::parallel_for(quick_sort_range_, quick_sort_body_);
}
//=================================================================================================//
} // namespace SPH
