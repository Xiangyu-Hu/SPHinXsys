#ifndef ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H
#define ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

class RegisterOphelieCoilFields : public LocalDynamics
{
  public:
    explicit RegisterOphelieCoilFields(SPHBody &sph_body, const OphelieCoilFieldNames &names)
        : LocalDynamics(sph_body)
    {
        auto &particles = particles_;
        particles->template registerStateVariable<Vecd>(names.j_src_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_src_imag, ZeroData<Vecd>::value);
        particles->addVariableToWrite<Vecd>(names.j_src_real);
    }
};

class RegisterOphelieGlassFields : public LocalDynamics
{
  public:
    explicit RegisterOphelieGlassFields(SPHBody &sph_body, const OphelieGlassFieldNames &names)
        : LocalDynamics(sph_body)
    {
        auto &particles = particles_;
        particles->template registerStateVariable<Real>(names.sigma, Real(0));
        particles->template registerStateVariable<Vecd>(names.a_coil_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.a_coil_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_coil_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_coil_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.a_ind_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.a_ind_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_ind_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_ind_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.a_src_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.a_src_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_src_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.b_src_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.e_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.e_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.joule_heat, Real(0));
        particles->template registerStateVariable<Real>(names.phi_real, Real(0));
        particles->template registerStateVariable<Vecd>(names.grad_phi_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.phi_rhs_real, Real(0));
        particles->template registerStateVariable<Real>(names.phi_lhs_real, Real(0));
        particles->template registerStateVariable<Real>(names.phi_imag, Real(0));
        particles->template registerStateVariable<Vecd>(names.grad_phi_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.phi_rhs_imag, Real(0));
        particles->template registerStateVariable<Real>(names.phi_lhs_imag, Real(0));
        particles->template registerStateVariable<Real>(names.phi_laplace_diag, Real(0));
        particles->template registerStateVariable<Real>(names.div_j_imag, Real(0));
        particles->template registerStateVariable<Real>(names.edge_flux_residual_imag, Real(0));
        particles->template registerStateVariable<Real>(names.edge_flux_residual_real, Real(0));
        particles->template registerStateVariable<Real>(names.joule_heat_edge, Real(0));
        particles->template registerStateVariable<Real>(names.power_edge_particle, Real(0));
        particles->template registerStateVariable<Real>(names.joule_heat_edge_recon_imag, Real(0));
        particles->template registerStateVariable<Real>(names.joule_heat_edge_recon_real, Real(0));
        particles->template registerStateVariable<Real>(names.joule_heat_edge_recon_complex, Real(0));
        particles->template registerStateVariable<Vecd>(names.e_imag_particle_diag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_imag_particle_diag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.joule_heat_particle_diag, Real(0));
        particles->template registerStateVariable<Vecd>(names.e_edge_recon_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_edge_recon_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.e_edge_recon_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(names.j_edge_recon_real, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.edge_recon_condition, Real(0));
        particles->template registerStateVariable<Real>(names.edge_recon_fallback, Real(0));
        particles->template registerStateVariable<Real>(names.edge_drop_abs_max, Real(0));
        particles->template registerStateVariable<Real>(names.edge_drop_sq_mean, Real(0));
        particles->template registerStateVariable<Real>(names.edge_q_antisym_max, Real(0));
        particles->template registerStateVariable<Real>(names.edge_q_antisym_sq_sum, Real(0));
        particles->template registerStateVariable<Real>(names.edge_q_scale_sq_sum, Real(0));
        particles->template registerStateVariable<Real>(names.edge_q_neighbor_count, Real(0));
        particles->template registerStateVariable<Real>(names.edge_q_nonfinite_count, Real(0));
        particles->template registerStateVariable<Real>(names.phi_boundary_mask, Real(0));
        particles->template registerStateVariable<Vecd>(names.phi_boundary_normal, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.phi_boundary_gn, Real(0));
        particles->template registerStateVariable<Matd>(names.phi_grad_linear_correction,
                                                        IdentityMatrix<Matd>::value);
        particles->template registerStateVariable<Real>(names.krylov_workspace, Real(0));
        particles->template registerStateVariable<Real>(names.krylov_basis, Real(0));
        particles->template registerStateVariable<Real>(names.krylov_scratch_a, Real(0));
        particles->template registerStateVariable<Real>(names.krylov_scratch_b, Real(0));
        particles->template registerStateVariable<Vecd>(kTeam7EEdgeTangentLsDiag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Vecd>(kTeam7JEdgeTangentLsDiag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(kTeam7EdgeTangentLsFallback, Real(0));

        particles->addVariableToWrite<Vecd>(names.a_coil_real);
        particles->addVariableToWrite<Vecd>(names.a_coil_imag);
        particles->addVariableToWrite<Vecd>(names.a_ind_real);
        particles->addVariableToWrite<Vecd>(names.a_ind_imag);
        particles->addVariableToWrite<Vecd>(names.a_src_real);
        particles->addVariableToWrite<Vecd>(names.a_src_imag);
        particles->addVariableToWrite<Vecd>(names.b_src_real);
        particles->addVariableToWrite<Vecd>(names.e_real);
        particles->addVariableToWrite<Vecd>(names.e_imag);
        particles->addVariableToWrite<Vecd>(names.j_real);
        particles->addVariableToWrite<Vecd>(names.j_imag);
        particles->addVariableToWrite<Real>(names.joule_heat);
        particles->addVariableToWrite<Real>(names.joule_heat_edge);
        particles->addVariableToWrite<Real>(names.joule_heat_edge_recon_imag);
        particles->addVariableToWrite<Real>(names.joule_heat_edge_recon_real);
        particles->addVariableToWrite<Real>(names.joule_heat_edge_recon_complex);
        particles->addVariableToWrite<Vecd>(names.e_edge_recon_imag);
        particles->addVariableToWrite<Vecd>(names.e_edge_recon_real);
        particles->addVariableToWrite<Vecd>(names.j_edge_recon_imag);
        particles->addVariableToWrite<Vecd>(names.j_edge_recon_real);
        particles->addVariableToWrite<Real>(names.edge_flux_residual_imag);
        particles->addVariableToWrite<Real>(names.edge_flux_residual_real);
        particles->addVariableToWrite<Real>(names.sigma);
        particles->addVariableToWrite<Real>(names.phi_imag);
        particles->addVariableToWrite<Real>(names.phi_real);
        particles->addVariableToWrite<Vecd>(names.grad_phi_imag);
        particles->addVariableToWrite<Vecd>(names.grad_phi_real);
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H
