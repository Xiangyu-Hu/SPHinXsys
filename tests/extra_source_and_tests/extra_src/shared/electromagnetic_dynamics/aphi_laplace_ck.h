#ifndef APHI_LAPLACE_CK_H
#define APHI_LAPLACE_CK_H

#include "base_general_dynamics.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename DataType>
struct AphiZeroValue;

template <>
struct AphiZeroValue<Real>
{
    static Real value() { return Real(0); }
};

template <>
struct AphiZeroValue<Vecd>
{
    static Vecd value() { return Vecd::Zero(); }
};

template <typename... RelationTypes>
class AphiPairwiseLaplaceCK;

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
class AphiPairwiseLaplaceCK<Base, DataType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    explicit AphiPairwiseLaplaceCK(RelationType<Parameters...> &relation, const std::string &input_name,
                                   const std::string &coefficient_name, const std::string &output_name,
                                   Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiPairwiseLaplaceCK(DynamicsArgs<RelationType<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiPairwiseLaplaceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiPairwiseLaplaceCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        DataType *input_;
        DataType *output_;
        Real *coefficient_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    std::string input_name_;
    std::string coefficient_name_;
    std::string output_name_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<DataType> *dv_input_;
    DiscreteVariable<DataType> *dv_output_;
    DiscreteVariable<Real> *dv_coefficient_;
};

template <typename DataType, typename... Parameters>
class AphiPairwiseLaplaceCK<Inner<DataType, Parameters...>>
    : public AphiPairwiseLaplaceCK<Base, DataType, Inner<Parameters...>>
{
    using BaseDynamicsType = AphiPairwiseLaplaceCK<Base, DataType, Inner<Parameters...>>;

  public:
    explicit AphiPairwiseLaplaceCK(Inner<Parameters...> &inner_relation, const std::string &input_name,
                                   const std::string &coefficient_name, const std::string &output_name,
                                   Real pair_weight_regularization = Real(0.01))
        : BaseDynamicsType(inner_relation, input_name, coefficient_name, output_name, pair_weight_regularization) {}
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiPairwiseLaplaceCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiPairwiseLaplaceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiPairwiseLaplaceCK() = default;
};

/** Contact pairwise Laplace: harmonic-mean coefficient, gather to owner i (no atomic). */
template <typename DataType, typename... Parameters>
class AphiPairwiseLaplaceCK<Contact<DataType, Parameters...>>
    : public AphiPairwiseLaplaceCK<Base, DataType, Contact<Parameters...>>
{
    using BaseDynamicsType = AphiPairwiseLaplaceCK<Base, DataType, Contact<Parameters...>>;

  public:
    explicit AphiPairwiseLaplaceCK(Contact<Parameters...> &contact_relation, const std::string &input_name,
                                   const std::string &coefficient_name, const std::string &output_name,
                                   Real pair_weight_regularization = Real(0.01))
        : BaseDynamicsType(contact_relation, input_name, coefficient_name, output_name, pair_weight_regularization)
    {
        for (auto *contact_particles : this->contact_particles_)
        {
            dv_contact_input_.push_back(contact_particles->template getVariableByName<DataType>(input_name));
            dv_contact_coefficient_.push_back(contact_particles->template getVariableByName<Real>(coefficient_name));
        }
    }
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiPairwiseLaplaceCK(DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiPairwiseLaplaceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiPairwiseLaplaceCK() = default;

    class InteractKernel : public BaseDynamicsType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataType *contact_input_;
        Real *contact_coefficient_;
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> dv_contact_input_;
    StdVec<DiscreteVariable<Real> *> dv_contact_coefficient_;
};

inline Real AphiHarmonicMean(Real lhs, Real rhs)
{
    return Real(2) * lhs * rhs / (lhs + rhs + TinyReal);
}

inline Real AphiPairwiseNegativeLaplaceWeight(Real dW_ijV_j, Real distance, Real distance_sq,
                                              Real pair_weight_regularization, Real reference_smoothing_length)
{
    const Real h_sq = reference_smoothing_length * reference_smoothing_length;
    return -Real(2) * distance * dW_ijV_j / (distance_sq + pair_weight_regularization * h_sq + TinyReal);
}

/** Uncorrected pairwise gradient weight; matches LinearGradient sign convention when B = I. */
inline Vecd AphiPairwiseGradientWeightUncorrected(Real dW_ijV_j, const Vecd &e_ij)
{
    return -dW_ijV_j * e_ij;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_LAPLACE_CK_H
