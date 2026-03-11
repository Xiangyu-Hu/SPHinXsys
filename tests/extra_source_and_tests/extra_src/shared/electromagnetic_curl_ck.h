#ifndef ELECTROMAGNETIC_CURL_CK_H
#define ELECTROMAGNETIC_CURL_CK_H

#include "interaction_ck.h"
#include "general_gradient.h"

namespace SPH
{
namespace electromagnetics
{
/**
 * @brief Curl operator for vector potential A in A-phi formulation.
 * Result variable name: "<variable_name>Curl".
 */
template <typename... RelationTypes>
class AphiCurlCK;

template <template <typename...> class RelationType, typename... Parameters>
class AphiCurlCK<Base, RelationType<Parameters...>>
    : public Gradient<Base, Vecd, RelationType<Parameters...>>
{
    using BaseDynamicsType = Gradient<Base, Vecd, RelationType<Parameters...>>;

  public:
    explicit AphiCurlCK(RelationType<Parameters...> &relation,
                        const std::string &variable_name = "VectorPotential");
    template <typename FirstArg>
    explicit AphiCurlCK(DynamicsArgs<RelationType<Parameters...>, FirstArg> parameters)
        : AphiCurlCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AphiCurlCK() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        AngularVecd *curl_;
    };

  protected:
    DiscreteVariable<AngularVecd> *dv_curl_;
};

template <typename... Parameters>
class AphiCurlCK<Inner<Parameters...>>
    : public AphiCurlCK<Base, Inner<Parameters...>>
{
    using BaseDynamicsType = AphiCurlCK<Base, Inner<Parameters...>>;

  public:
    explicit AphiCurlCK(Inner<Parameters...> &inner_relation,
                        const std::string &variable_name = "VectorPotential")
        : BaseDynamicsType(inner_relation, variable_name) {};
    template <typename FirstArg>
    explicit AphiCurlCK(DynamicsArgs<Inner<Parameters...>, FirstArg> parameters)
        : AphiCurlCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AphiCurlCK() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseDynamicsType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };
};

template <typename... Parameters>
class AphiCurlCK<Contact<Parameters...>>
    : public AphiCurlCK<Base, Contact<Parameters...>>
{
    using BaseDynamicsType = AphiCurlCK<Base, Contact<Parameters...>>;

  public:
    explicit AphiCurlCK(Contact<Parameters...> &contact_relation,
                        const std::string &variable_name = "VectorPotential");
    template <typename FirstArg>
    explicit AphiCurlCK(DynamicsArgs<Contact<Parameters...>, FirstArg> parameters)
        : AphiCurlCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AphiCurlCK() {};

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        Vecd *contact_variable_;
    };

  protected:
    StdVec<DiscreteVariable<Vecd> *> dv_contact_variable_;
};

using AphiCurlCKInner = AphiCurlCK<Inner<>>;
using AphiCurlCKComplex = ComplexInteraction<AphiCurlCK<Inner<>, Contact<>>>;
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_CURL_CK_H
