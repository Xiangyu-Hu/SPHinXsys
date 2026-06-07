#ifndef APHI_PAIRWISE_DIV_A_CK_H
#define APHI_PAIRWISE_DIV_A_CK_H

#include "base_general_dynamics.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... RelationTypes>
class AphiPairwiseVectorDivergenceCK;

/** OPERATOR default: divA_i = sum_j g_ij · (A_i - A_j), g_ij = uncorrected pairwise gradient weight. */
template <typename... Parameters>
class AphiPairwiseVectorDivergenceCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiPairwiseVectorDivergenceCK(Inner<Parameters...> &inner_relation, const std::string &input_name,
                                            const std::string &output_name);
    template <typename FirstArg, typename SecondArg>
    explicit AphiPairwiseVectorDivergenceCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiPairwiseVectorDivergenceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                         std::get<1>(parameters.others_)){};
    virtual ~AphiPairwiseVectorDivergenceCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Vecd *input_;
        Real *output_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Real> *dv_output_;
};

/** Contact pairwise divA: gather cross-body g_ij·(A_i - A_j) to owner output (no atomic). */
template <typename... Parameters>
class AphiPairwiseVectorDivergenceCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiPairwiseVectorDivergenceCK(Contact<Parameters...> &contact_relation, const std::string &input_name,
                                            const std::string &output_name);
    template <typename FirstArg, typename SecondArg>
    explicit AphiPairwiseVectorDivergenceCK(DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiPairwiseVectorDivergenceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                         std::get<1>(parameters.others_)){};
    virtual ~AphiPairwiseVectorDivergenceCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Vecd *input_;
        Real *output_;
        Real *contact_Vol_;
        Vecd *contact_input_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Real> *dv_output_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_input_;
};

template <typename... RelationTypes>
class AphiPairwiseScalarGradientCK;

template <typename... RelationTypes>
class AphiPairwiseVectorGradientCK;

/** OPERATOR default: grad f_i = sum_j g_ij (f_i - f_j) for scalar f. */
template <typename... Parameters>
class AphiPairwiseScalarGradientCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiPairwiseScalarGradientCK(Inner<Parameters...> &inner_relation, const std::string &input_name,
                                          const std::string &output_name);
    template <typename FirstArg, typename SecondArg>
    explicit AphiPairwiseScalarGradientCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiPairwiseScalarGradientCK(parameters.identifier_, std::get<0>(parameters.others_),
                                       std::get<1>(parameters.others_)){};
    virtual ~AphiPairwiseScalarGradientCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *input_;
        Vecd *output_;
    };

  protected:
    DiscreteVariable<Real> *dv_input_;
    DiscreteVariable<Vecd> *dv_output_;
};

/** Contact pairwise scalar grad: gather cross-body g_ij (f_i - f_j) to owner output (no atomic). */
template <typename... Parameters>
class AphiPairwiseScalarGradientCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiPairwiseScalarGradientCK(Contact<Parameters...> &contact_relation, const std::string &input_name,
                                          const std::string &output_name);
    template <typename FirstArg, typename SecondArg>
    explicit AphiPairwiseScalarGradientCK(DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiPairwiseScalarGradientCK(parameters.identifier_, std::get<0>(parameters.others_),
                                       std::get<1>(parameters.others_)){};
    virtual ~AphiPairwiseScalarGradientCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *input_;
        Vecd *output_;
        Real *contact_Vol_;
        Real *contact_input_;
    };

  protected:
    DiscreteVariable<Real> *dv_input_;
    DiscreteVariable<Vecd> *dv_output_;
    StdVec<DiscreteVariable<Real> *> dv_contact_input_;
};

template <typename... RelationTypes>
class AphiPairwiseVectorGradientCK;

/** OPERATOR default: grad A_i = sum_j g_ij (A_i - A_j)^T as row-wise outer product (uncorrected pairwise). */
template <typename... Parameters>
class AphiPairwiseVectorGradientCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiPairwiseVectorGradientCK(Inner<Parameters...> &inner_relation, const std::string &input_name,
                                          const std::string &output_name);
    template <typename FirstArg, typename SecondArg>
    explicit AphiPairwiseVectorGradientCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiPairwiseVectorGradientCK(parameters.identifier_, std::get<0>(parameters.others_),
                                       std::get<1>(parameters.others_)){};
    virtual ~AphiPairwiseVectorGradientCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Vecd *input_;
        Matd *output_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Matd> *dv_output_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PAIRWISE_DIV_A_CK_H
