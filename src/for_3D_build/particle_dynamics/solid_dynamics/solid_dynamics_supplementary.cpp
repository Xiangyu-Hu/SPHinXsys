#include "all_solid_dynamics.h"
#include "base_data_package.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "neighborhood.h"
#include "polar_decomposition_3x3.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "weakly_compressible_fluid.h"

using namespace polar;

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
void UpdateElasticNormalDirection::update(size_t index_i, Real dt)
{
    Matd &F = F_[index_i];
    Mat3d R;
    Real Q[9], H[9], A[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i * 3 + j] = F(i, j);

    polar::polar_decomposition(Q, H, A);
    // this decomposition has the form A = Q*H, where Q is orthogonal and H is symmetric positive semi-definite.
    // Ref. "An algorithm to compute the polar decomposition of a 3*3 matrix, Nicholas J. Higham et al. Numer Algor(2016) "
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R(i, j) = Q[i * 3 + j];
    n_[index_i] = R * n0_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
  //=====================================================================================================//
} // namespace SPH
