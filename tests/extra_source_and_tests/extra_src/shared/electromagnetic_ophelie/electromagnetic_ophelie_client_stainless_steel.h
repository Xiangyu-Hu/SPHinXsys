/**
 * @file electromagnetic_ophelie_client_stainless_steel.h
 * @brief Doc228 development defaults for 304 / 1.4301 stainless (room-T).
 *
 * source_status=development_assumption — not a high-T paddle strength certificate.
 * Reference: Outokumpu Core 304/4301 datasheet (see Doc228).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_CLIENT_STAINLESS_STEEL_H
#define ELECTROMAGNETIC_OPHELIE_CLIENT_STAINLESS_STEEL_H

#include "sphinxsys.h"

#include <fstream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieClientStainless304Development
{
    std::string material_id = "stainless_304_1p4301_room_T_dev";
    std::string source_status = "development_assumption";
    std::string reference =
        "Outokumpu Core range datasheet Core 304/4301 (Doc228); room-T engineering defaults";
    Real rho_kg_m3 = Real(7900);
    Real youngs_modulus_pa = Real(2.0e11); // 200 GPa
    Real poisson = Real(0.30);
    Real k_w_m_k = Real(15);
    Real cp_j_kg_k = Real(500);
    Real alpha_1_k = Real(16.0e-6); // reserved; thermal strain not enabled by default
    Real t_ref_k = Real(293.15);
};

inline OphelieClientStainless304Development makeClientStainless304Development()
{
    return {};
}

inline bool writeClientStainless304Json(const OphelieClientStainless304Development &m, const std::string &path)
{
    std::ofstream out(path);
    if (!out)
        return false;
    out << "{\n";
    out << "  \"material_id\": \"" << m.material_id << "\",\n";
    out << "  \"source_status\": \"" << m.source_status << "\",\n";
    out << "  \"reference\": \"" << m.reference << "\",\n";
    out << "  \"rho_kg_m3\": " << m.rho_kg_m3 << ",\n";
    out << "  \"youngs_modulus_Pa\": " << m.youngs_modulus_pa << ",\n";
    out << "  \"poisson\": " << m.poisson << ",\n";
    out << "  \"k_W_per_mK\": " << m.k_w_m_k << ",\n";
    out << "  \"cp_J_per_kgK\": " << m.cp_j_kg_k << ",\n";
    out << "  \"alpha_1_per_K\": " << m.alpha_1_k << ",\n";
    out << "  \"T_ref_K\": " << m.t_ref_k << ",\n";
    out << "  \"note\": \"S0 may use softer rubber-like E for analytic SPH path gate; "
           "paddle E0/E1 must use this stainless preset without secretly lowering E\"\n";
    out << "}\n";
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_CLIENT_STAINLESS_STEEL_H
