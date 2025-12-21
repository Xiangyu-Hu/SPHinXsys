#include "level_set_transformation.hpp"

namespace SPH
{
//=============================================================================================//
UpdateLevelSetGradient::UpdateLevelSetGradient(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerPackageVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
UpdateKernelIntegrals::UpdateKernelIntegrals(
    SparseMeshField<4> &data_mesh, UnsignedInt resolution_level,
    NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method)
    : BaseMeshLocalDynamics(data_mesh, resolution_level),
      neighbor_method_(neighbor_method),
      mv_phi_(*data_mesh.getPackageVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getPackageVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      mv_kernel_weight_(*data_mesh.registerPackageVariable<Real>("KernelWeight")),
      mv_kernel_gradient_(*data_mesh.registerPackageVariable<Vecd>("KernelGradient")),
      mv_kernel_second_gradient_(*data_mesh.registerPackageVariable<Matd>("KernelSecondGradient"))
{
    data_mesh.setBoundaryData(&mv_kernel_weight_, resolution_level, [&](UnsignedInt k)
                              { return PackageVariableData<Real>::Constant(k == 0 ? 0.0 : 1.0); });
    data_mesh.setBoundaryData(&mv_kernel_gradient_, resolution_level, [&](UnsignedInt k)
                              { return PackageVariableData<Vecd>::Constant(Vecd::Zero()); });
    data_mesh.setBoundaryData(&mv_kernel_second_gradient_, resolution_level, [&](UnsignedInt k)
                              { return PackageVariableData<Matd>::Constant(Matd::Zero()); });
}
//=================================================================================================//
} // namespace SPH
