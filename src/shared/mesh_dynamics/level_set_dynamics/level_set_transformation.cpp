#include "level_set_transformation.hpp"

namespace SPH
{
//=============================================================================================//
UpdateLevelSetGradient::UpdateLevelSetGradient(MeshWithGridDataPackagesType &data_mesh)
    : BaseMeshLocalDynamics(data_mesh),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.registerMeshVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()) {}
//=============================================================================================//
UpdateKernelIntegrals::UpdateKernelIntegrals(
    MeshWithGridDataPackagesType &data_mesh, NeighborMethod<SPHAdaptation, SPHAdaptation> &neighbor_method)
    : BaseMeshLocalDynamics(data_mesh), neighbor_method_(neighbor_method),
      mv_phi_(*data_mesh.getMeshVariable<Real>("LevelSet")),
      mv_phi_gradient_(*data_mesh.getMeshVariable<Vecd>("LevelSetGradient")),
      dv_cell_neighborhood_(data_mesh.getCellNeighborhood()),
      mv_kernel_weight_(*data_mesh.registerMeshVariable<Real>("KernelWeight")),
      mv_kernel_gradient_(*data_mesh.registerMeshVariable<Vecd>("KernelGradient")),
      mv_kernel_second_gradient_(*data_mesh.registerMeshVariable<Matd>("KernelSecondGradient"))
{
    IndexHandler &index_handler = data_mesh.getIndexHandler();
    Real far_field_distance = index_handler.GridSpacing() * (Real)index_handler.BufferWidth();
    initializeSingularPackages(0, -far_field_distance);
    initializeSingularPackages(1, far_field_distance);
}
//=================================================================================================//
void UpdateKernelIntegrals::initializeSingularPackages(
    UnsignedInt package_index, Real far_field_level_set)
{
    auto &kernel_weight = mv_kernel_weight_.Data()[package_index];
    auto &kernel_gradient = mv_kernel_gradient_.Data()[package_index];
    auto &kernel_second_gradient = mv_kernel_second_gradient_.Data()[package_index];

    mesh_for_each(Arrayi::Zero(), Arrayi::Constant(pkg_size),
                  [&](const Arrayi &data_index)
                  {
                      kernel_weight(data_index) = far_field_level_set < 0.0 ? 0 : 1.0;
                      kernel_gradient(data_index) = Vecd::Zero();
                      kernel_second_gradient(data_index) = Matd::Zero();
                  });
}
//=================================================================================================//
} // namespace SPH
