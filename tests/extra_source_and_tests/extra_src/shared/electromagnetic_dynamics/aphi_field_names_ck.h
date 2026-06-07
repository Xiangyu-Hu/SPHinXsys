#ifndef APHI_FIELD_NAMES_CK_H
#define APHI_FIELD_NAMES_CK_H

#include <string>

namespace SPH
{
namespace electromagnetics
{

struct AphiBlockNames
{
    std::string a_real;
    std::string a_imag;
    std::string phi_real;
    std::string phi_imag;
};

struct AphiMaterialNames
{
    std::string sigma = "Sigma";
    std::string nu = "Nu";
};

struct AphiDiagnosticNames
{
    std::string div_a_real = "DivAReal";
    std::string div_a_imag = "DivAImag";
    AphiBlockNames laplace = {"LapAReal", "LapAImag", "LapPhiReal", "LapPhiImag"};
};

struct AphiVariableNames
{
    AphiBlockNames solution = {"AReal", "AImag", "PhiReal", "PhiImag"};
    AphiBlockNames rhs = {"RhsAReal", "RhsAImag", "RhsPhiReal", "RhsPhiImag"};
    AphiBlockNames lhs = {"LhsAReal", "LhsAImag", "LhsPhiReal", "LhsPhiImag"};
    AphiBlockNames residual = {"ResidualAReal", "ResidualAImag", "ResidualPhiReal", "ResidualPhiImag"};
    AphiBlockNames true_residual = {"TrueResAReal", "TrueResAImag", "TrueResPhiReal", "TrueResPhiImag"};
    AphiBlockNames r_hat = {"RHatAReal", "RHatAImag", "RHatPhiReal", "RHatPhiImag"};
    AphiBlockNames search = {"SearchAReal", "SearchAImag", "SearchPhiReal", "SearchPhiImag"};
    AphiBlockNames v = {"VAReal", "VAImag", "VPhiReal", "VPhiImag"};
    AphiBlockNames v_old = {"VOldAReal", "VOldAImag", "VOldPhiReal", "VOldPhiImag"};
    AphiBlockNames s = {"SAReal", "SAImag", "SPhiReal", "SPhiImag"};
    AphiBlockNames t = {"TAReal", "TAImag", "TPhiReal", "TPhiImag"};
    /** z = M^{-1} p in right-preconditioned BiCGStab. */
    AphiBlockNames z = {"ZAReal", "ZAImag", "ZPhiReal", "ZPhiImag"};
    /** y = M^{-1} s in right-preconditioned BiCGStab. */
    AphiBlockNames y = {"YAReal", "YAImag", "YPhiReal", "YPhiImag"};
    AphiMaterialNames material = {};
    AphiDiagnosticNames diagnostic = {};
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_FIELD_NAMES_CK_H
