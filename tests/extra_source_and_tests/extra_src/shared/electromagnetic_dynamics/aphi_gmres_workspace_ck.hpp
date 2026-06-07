#ifndef APHI_GMRES_WORKSPACE_CK_HPP
#define APHI_GMRES_WORKSPACE_CK_HPP

#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline RegisterAphiGMRESWorkspaceCK::RegisterAphiGMRESWorkspaceCK(SPHBody &sph_body, UnsignedInt restart_dimension)
    : LocalDynamics(sph_body), restart_dimension_(restart_dimension),
      workspace_(buildAphiGMRESWorkspaceNames(restart_dimension))
{
    for (const AphiBlockNames &block_names : workspace_.v_basis)
    {
        registerBlock(block_names);
    }
    for (const AphiBlockNames &block_names : workspace_.z_basis)
    {
        registerBlock(block_names);
    }
}

inline void RegisterAphiGMRESWorkspaceCK::registerBlock(const AphiBlockNames &block_names)
{
    particles_->registerStateVariable<Vecd>(block_names.a_real);
    particles_->registerStateVariable<Vecd>(block_names.a_imag);
    particles_->registerStateVariable<Real>(block_names.phi_real);
    particles_->registerStateVariable<Real>(block_names.phi_imag);
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_WORKSPACE_CK_HPP
