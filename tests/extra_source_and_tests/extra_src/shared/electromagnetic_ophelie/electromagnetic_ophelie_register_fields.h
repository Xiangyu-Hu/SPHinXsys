#ifndef ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H
#define ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_field_names.h"

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
        particles->template registerStateVariable<Real>(names.phi_imag, Real(0));
        particles->template registerStateVariable<Vecd>(names.grad_phi_imag, ZeroData<Vecd>::value);
        particles->template registerStateVariable<Real>(names.phi_rhs_imag, Real(0));
        particles->template registerStateVariable<Real>(names.phi_lhs_imag, Real(0));
        particles->template registerStateVariable<Real>(names.phi_laplace_diag, Real(0));
        particles->template registerStateVariable<Real>(names.div_j_imag, Real(0));

        particles->addVariableToWrite<Vecd>(names.a_coil_real);
        particles->addVariableToWrite<Vecd>(names.a_ind_real);
        particles->addVariableToWrite<Vecd>(names.a_src_real);
        particles->addVariableToWrite<Vecd>(names.b_src_real);
        particles->addVariableToWrite<Vecd>(names.e_imag);
        particles->addVariableToWrite<Vecd>(names.j_imag);
        particles->addVariableToWrite<Real>(names.joule_heat);
        particles->addVariableToWrite<Real>(names.sigma);
        particles->addVariableToWrite<Real>(names.phi_imag);
        particles->addVariableToWrite<Vecd>(names.grad_phi_imag);
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_REGISTER_FIELDS_H
