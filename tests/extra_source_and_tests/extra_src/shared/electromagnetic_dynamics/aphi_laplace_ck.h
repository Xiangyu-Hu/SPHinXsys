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
                                   Real pair_weight_regularization = Real(0));
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiPairwiseLaplaceCK(DynamicsArgs<RelationType<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiPairwiseLaplaceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiPairwiseLaplaceCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

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
                                   Real pair_weight_regularization = Real(0))
        : BaseDynamicsType(inner_relation, input_name, coefficient_name, output_name, pair_weight_regularization) {}
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiPairwiseLaplaceCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiPairwiseLaplaceCK(parameters.identifier_, std::get<0>(parameters.others_),
                                std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiPairwiseLaplaceCK() = default;
};

inline Real AphiHarmonicMean(Real lhs, Real rhs)
{
    return Real(2) * lhs * rhs / (lhs + rhs + TinyReal);
}

inline Real AphiPairwiseNegativeLaplaceWeight(Real dW_ijV_j, Real distance,
                                              Real pair_weight_regularization, Real reference_smoothing_length)
{
    return -Real(2) * dW_ijV_j /
           (distance + pair_weight_regularization * reference_smoothing_length + TinyReal);
}

/** Uncorrected pairwise gradient weight; matches LinearGradient sign convention when B = I. */
inline Vecd AphiPairwiseGradientWeightUncorrected(Real dW_ijV_j, const Vecd &e_ij)
{
    return -dW_ijV_j * e_ij;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_LAPLACE_CK_H
