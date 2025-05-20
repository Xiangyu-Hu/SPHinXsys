
#ifndef IO_PLT_HPP
#define IO_PLT_HPP

#include "io_plt.h"

namespace SPH
{
//=============================================================================================//
template <int N>
void PltEngine::writeAQuantityHeader(
    std::ofstream &out_file, const Eigen::Matrix<Real, N, 1> &quantity, const std::string &quantity_name)
{
    for (int i = 0; i != N; ++i)
        out_file << "\"" << quantity_name << "[" << i << "]\"" << "   ";
}
//=============================================================================================//
template <int N, int M>
void PltEngine::writeAQuantityHeader(
    std::ofstream &out_file, const Eigen::Matrix<Real, N, M> &quantity, const std::string &quantity_name)
{
    for (int i = 0; i != N; ++i)
        for (int j = 0; j != M; ++j)
            out_file << "\"" << quantity_name << "[" << i << "][" << j << "]\"" << "   ";
}
//=============================================================================================//
template <int N>
void PltEngine::writeAQuantity(std::ofstream &out_file, const Eigen::Matrix<Real, N, 1> &quantity)
{
    for (int i = 0; i < N; ++i)
        out_file << std::fixed << std::setprecision(9) << quantity[i] << "   ";
}
//=============================================================================================//
template <int N, int M>
void PltEngine::writeAQuantity(std::ofstream &out_file, const Eigen::Matrix<Real, N, M> &quantity)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            out_file << std::fixed << std::setprecision(9) << quantity(i, j) << "   ";
}
//=================================================================================================//
} // namespace SPH
#endif // IO_PLT_HPP