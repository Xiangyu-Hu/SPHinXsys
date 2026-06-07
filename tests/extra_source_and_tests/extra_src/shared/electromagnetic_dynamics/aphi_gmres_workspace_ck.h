#ifndef APHI_GMRES_WORKSPACE_CK_H
#define APHI_GMRES_WORKSPACE_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{

/** Maximum supported GMRES restart dimension m (registers GMRESV0..GMRESV{m}). */
static constexpr UnsignedInt AphiGMRESMaxRestartDimension = 80;

inline bool aphiGMRESRestartDimensionValid(UnsignedInt restart_dimension)
{
    return restart_dimension > 0 && restart_dimension <= AphiGMRESMaxRestartDimension;
}

inline AphiBlockNames aphiGMRESVBlock(UnsignedInt index)
{
    const std::string prefix = "GMRESV" + std::to_string(index);
    return AphiBlockNames{prefix + "AReal", prefix + "AImag", prefix + "PhiReal", prefix + "PhiImag"};
}

inline AphiBlockNames aphiGMRESZBlock(UnsignedInt index)
{
    const std::string prefix = "GMRESZ" + std::to_string(index);
    return AphiBlockNames{prefix + "AReal", prefix + "AImag", prefix + "PhiReal", prefix + "PhiImag"};
}

struct AphiGMRESWorkspaceNames
{
    std::vector<AphiBlockNames> v_basis;
    std::vector<AphiBlockNames> z_basis;
};

inline AphiGMRESWorkspaceNames buildAphiGMRESWorkspaceNames(UnsignedInt restart_dimension)
{
    AphiGMRESWorkspaceNames workspace;
    for (UnsignedInt index = 0; index <= restart_dimension; ++index)
    {
        workspace.v_basis.push_back(aphiGMRESVBlock(index));
    }
    for (UnsignedInt index = 0; index < restart_dimension; ++index)
    {
        workspace.z_basis.push_back(aphiGMRESZBlock(index));
    }
    return workspace;
}

class RegisterAphiGMRESWorkspaceCK : public LocalDynamics
{
  public:
    RegisterAphiGMRESWorkspaceCK(SPHBody &sph_body, UnsignedInt restart_dimension);
    virtual ~RegisterAphiGMRESWorkspaceCK() = default;

  protected:
    void registerBlock(const AphiBlockNames &block_names);

    UnsignedInt restart_dimension_;
    AphiGMRESWorkspaceNames workspace_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_WORKSPACE_CK_H
