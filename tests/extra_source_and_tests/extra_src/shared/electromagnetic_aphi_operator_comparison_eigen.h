#ifndef ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_H
#define ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_H

#include "inner_body_relation.h"
#include <Eigen/Sparse>
#include <complex>

namespace SPH
{
namespace electromagnetics
{
using Complex = std::complex<Real>;
using VecC = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using SpMatC = Eigen::SparseMatrix<Complex>;
using TripletC = Eigen::Triplet<Complex>;

struct DiscreteOperatorComparisonParameters
{
    Real pair_weight_regularization = 0.01;
    Real gauge_penalty = 0.0;
    bool enable_gauge_penalty = false;
};

struct DiscreteModeEnergySummary
{
    Real laplace_energy = 0.0;
    Real weak_curlcurl_energy = 0.0;
    Real divergence_l2 = 0.0;
    Real curl_l2 = 0.0;
    Real max_abs_a = 0.0;
};

class DiscreteEMOperatorComparator
{
  public:
    explicit DiscreteEMOperatorComparator(BaseInnerRelation &inner_relation,
                                          const DiscreteOperatorComparisonParameters &parameters);

    void buildOperators(const StdVec<Real> &magnetic_reluctivity);

    VecC gatherUnknownVector(const StdVec<Complex> &ax,
                             const StdVec<Complex> &ay,
                             const StdVec<Complex> &az) const;

    DiscreteModeEnergySummary evaluateMode(const StdVec<Complex> &ax,
                                           const StdVec<Complex> &ay,
                                           const StdVec<Complex> &az) const;

    int idxAx(size_t particle_index) const { return static_cast<int>(particle_index); }
    int idxAy(size_t particle_index) const { return static_cast<int>(number_of_particles_ + particle_index); }
    int idxAz(size_t particle_index) const { return static_cast<int>(2 * number_of_particles_ + particle_index); }
    size_t NumberOfParticles() const { return number_of_particles_; }

  protected:
    Real pairwiseDiffusionWeight(const Neighborhood &inner_neighborhood,
                                 size_t neighbor_index) const;
    static Real harmonicMean(Real value_i, Real value_j);

    BaseInnerRelation &inner_relation_;
    DiscreteOperatorComparisonParameters parameters_;
    size_t number_of_particles_;
    Real reference_smoothing_length_;

    Real *volumetric_measure_;

    SpMatC grad_x_;
    SpMatC grad_y_;
    SpMatC grad_z_;
    SpMatC laplace_nu_;
    SpMatC laplace_operator_;
    SpMatC weak_curlcurl_operator_;
    SpMatC divergence_matrix_;
    SpMatC curl_matrix_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_OPERATOR_COMPARISON_EIGEN_H
