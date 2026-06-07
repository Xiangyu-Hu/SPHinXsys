#ifndef APHI_JOULE_HEATING_CK_H
#define APHI_JOULE_HEATING_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

struct AphiJouleHeatingFieldNames
{
    std::string grad_phi_real = "GradPhiReal";
    std::string grad_phi_imag = "GradPhiImag";
    std::string electric_field_a_real = "ElectricFieldAReal";
    std::string electric_field_a_imag = "ElectricFieldAImag";
    std::string current_density_real = "CurrentDensityReal";
    std::string current_density_imag = "CurrentDensityImag";
    std::string joule_heat_source = "JouleHeatSource";
};

/** Uncorrected pairwise scalar gradient of phi (matches grad(phi) operator stencil). */
template <typename... RelationTypes>
class AphiComputeScalarPhiGradientCK;

template <typename... Parameters>
class AphiComputeScalarPhiGradientCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiComputeScalarPhiGradientCK(Inner<Parameters...> &inner_relation, const AphiBlockNames &phi_block,
                                            const AphiJouleHeatingFieldNames &field_names);
    template <typename FirstArg, typename SecondArg>
    explicit AphiComputeScalarPhiGradientCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiComputeScalarPhiGradientCK(parameters.identifier_, std::get<0>(parameters.others_),
                                         std::get<1>(parameters.others_)){};
    virtual ~AphiComputeScalarPhiGradientCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *phi_real_;
        Real *phi_imag_;
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
};

/** Contact contribution to uncorrected pairwise scalar grad(phi). Inner assigns; Contact adds. */
template <typename... Parameters>
class AphiComputeScalarPhiGradientCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiComputeScalarPhiGradientCK(Contact<Parameters...> &contact_relation, const AphiBlockNames &phi_block,
                                            const AphiJouleHeatingFieldNames &field_names);
    template <typename FirstArg, typename SecondArg>
    explicit AphiComputeScalarPhiGradientCK(DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiComputeScalarPhiGradientCK(parameters.identifier_, std::get<0>(parameters.others_),
                                         std::get<1>(parameters.others_)){};
    virtual ~AphiComputeScalarPhiGradientCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        Real *phi_real_;
        Real *phi_imag_;
        Real *contact_phi_real_;
        Real *contact_phi_imag_;
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    StdVec<DiscreteVariable<Real> *> dv_contact_phi_real_;
    StdVec<DiscreteVariable<Real> *> dv_contact_phi_imag_;
};

/** E = omega * A_imag - grad(phi_real), E_imag = -omega * A_real - grad(phi_imag). */
class AphiComputeFrequencyElectricFieldCK : public LocalDynamics
{
  public:
    AphiComputeFrequencyElectricFieldCK(SPHBody &sph_body, Real omega, const AphiBlockNames &solution_block,
                                      const AphiJouleHeatingFieldNames &field_names);
    virtual ~AphiComputeFrequencyElectricFieldCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real omega_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
        Vecd *electric_field_a_real_;
        Vecd *electric_field_a_imag_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    DiscreteVariable<Vecd> *dv_electric_field_a_real_;
    DiscreteVariable<Vecd> *dv_electric_field_a_imag_;
};

/** J = sigma E, joule = 0.5 Re(J* · E) on each particle (9D-lite plumbing). */
class AphiComputeJouleHeatSourceCK : public LocalDynamics
{
  public:
    AphiComputeJouleHeatSourceCK(SPHBody &sph_body, const AphiMaterialNames &material_names,
                                 const AphiJouleHeatingFieldNames &field_names);
    virtual ~AphiComputeJouleHeatSourceCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *sigma_;
        Vecd *electric_field_a_real_;
        Vecd *electric_field_a_imag_;
        Vecd *current_density_real_;
        Vecd *current_density_imag_;
        Real *joule_heat_source_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_electric_field_a_real_;
    DiscreteVariable<Vecd> *dv_electric_field_a_imag_;
    DiscreteVariable<Vecd> *dv_current_density_real_;
    DiscreteVariable<Vecd> *dv_current_density_imag_;
    DiscreteVariable<Real> *dv_joule_heat_source_;
};

class RegisterAphiJouleHeatingFieldsCK : public LocalDynamics
{
  public:
    explicit RegisterAphiJouleHeatingFieldsCK(SPHBody &sph_body, const AphiJouleHeatingFieldNames &field_names);
    virtual ~RegisterAphiJouleHeatingFieldsCK() = default;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_JOULE_HEATING_CK_H
