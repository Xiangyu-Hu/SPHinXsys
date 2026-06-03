#ifndef ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H
#define ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H

#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieCoilFieldNames
{
    std::string j_src_real = "JSrcReal";
    std::string j_src_imag = "JSrcImag";
};

struct OphelieGlassFieldNames
{
    std::string sigma = "Sigma";

    std::string a_coil_real = "ACoilReal";
    std::string a_coil_imag = "ACoilImag";
    std::string b_coil_real = "BCoilReal";
    std::string b_coil_imag = "BCoilImag";

    std::string a_ind_real = "AIndReal";
    std::string a_ind_imag = "AIndImag";
    std::string b_ind_real = "BIndReal";
    std::string b_ind_imag = "BIndImag";

    /** Working total vector potential: A_coil + A_ind. */
    std::string a_src_real = "ASrcReal";
    std::string a_src_imag = "ASrcImag";
    std::string b_src_real = "BSrcReal";
    std::string b_src_imag = "BSrcImag";

    std::string e_real = "EReal";
    std::string e_imag = "EImag";
    std::string j_real = "JReal";
    std::string j_imag = "JImag";
    std::string joule_heat = "JouleHeat";

    std::string phi_imag = "PhiImag";
    std::string grad_phi_imag = "GradPhiImag";
    std::string phi_rhs_imag = "PhiRhsImag";
    std::string phi_lhs_imag = "PhiLhsImag";
    std::string phi_laplace_diag = "PhiLaplaceDiag";
    /** Diagnostic SPH divergence of JImag (not a governing unknown). */
    std::string div_j_imag = "DivJImag";
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FIELD_NAMES_H
