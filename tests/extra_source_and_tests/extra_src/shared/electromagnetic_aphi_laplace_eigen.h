#ifndef ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_H
#define ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_H

#include "inner_body_relation.h"
#include <Eigen/Sparse>
#include <complex>
#include <string>

namespace SPH
{
namespace electromagnetics
{
using Complex = std::complex<Real>;
using VecC = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using SpMatC = Eigen::SparseMatrix<Complex>;
using TripletC = Eigen::Triplet<Complex>;
using Vec3c = Eigen::Matrix<Complex, Dimensions, 1>;

struct LaplaceStructuredAPhiParameters
{
    Real angular_frequency = 0.0;
    Real geom_length_to_m = 1.0;
    Real gauge_penalty = 0.0;
    bool enable_gauge_penalty = false;
    bool enable_block_scaling = false;
    Real phi_block_scale = 0.0;
    bool enable_diagonal_equilibration = false;
    int diagonal_equilibration_iterations = 2;
    std::string solver_backend = "sparse_lu";
    int iterative_max_iterations = 5000;
    Real iterative_tolerance = 1.0e-8;
    bool fix_phi_reference = true;
    size_t phi_reference_index = 0;
    Real phi_reference_value = 0.0;
    Real pair_weight_regularization = 0.01;
    Real contact_gradient_scale = 1.0;
    Real contact_diffusion_scale = 1.0;
};

struct LaplaceStructuredAPhiBoundaryCondition
{
    StdVec<bool> is_dirichlet;
    StdVec<Complex> ax;
    StdVec<Complex> ay;
    StdVec<Complex> az;
    StdVec<Complex> phi;
};

struct LaplaceStructuredAPhiFields
{
    StdVec<Complex> ax;
    StdVec<Complex> ay;
    StdVec<Complex> az;
    StdVec<Complex> phi;

    StdVec<Vec3c> electric_field;
    StdVec<Vec3c> magnetic_flux_density;
    StdVec<Vec3c> current_density;
    StdVec<Real> joule_heat;
};

struct LaplaceStructuredAPhiDiagnostics
{
    Real linear_residual_norm = 0.0;
    Real divergence_a_l2 = 0.0;
    Real current_continuity_l2 = 0.0;
    Real phi_equation_residual_l2 = 0.0;
    Real physical_current_continuity_l2 = 0.0;
    Real relative_current_continuity_l2 = 0.0;
    Real phi_diffusion_term_l2 = 0.0;
    Real sigma_a_coupling_term_l2 = 0.0;
    Real rhs_phi_l2 = 0.0;
    Real max_abs_a = 0.0;
    Real max_abs_e = 0.0;
    Real max_abs_j = 0.0;
    Real total_joule_power = 0.0;
    Real gauge_penalty_energy = 0.0;
    Real magnetic_energy = 0.0;
};

struct LaplaceStructuredAPhiDiscreteView
{
    size_t number_of_particles = 0;
    Real reference_smoothing_length = 0.0;
    Vecd *positions = nullptr;
    Real *volumetric_measure = nullptr;
    Real *smoothing_length_ratio = nullptr;
    const ParticleConfiguration *particle_configuration = nullptr;
    const StdVec<StdVec<uint8_t>> *neighbor_is_contact = nullptr;
};

class LaplaceStructuredAPhiEigenSolver
{
  public:
    explicit LaplaceStructuredAPhiEigenSolver(BaseInnerRelation &inner_relation,
                                              const LaplaceStructuredAPhiParameters &parameters);
    explicit LaplaceStructuredAPhiEigenSolver(const LaplaceStructuredAPhiDiscreteView &discrete_view,
                                              const LaplaceStructuredAPhiParameters &parameters);

    void assembleSystem(const StdVec<Real> &electrical_conductivity,
                        const StdVec<Real> &magnetic_reluctivity,
                        const VecC &rhs_unconstrained,
                        const LaplaceStructuredAPhiBoundaryCondition &boundary_condition);

    VecC assembleRightHandSideFromExactSolution(const StdVec<Real> &electrical_conductivity,
                                                const StdVec<Real> &magnetic_reluctivity,
                                                const LaplaceStructuredAPhiFields &exact_fields);

    bool solve(LaplaceStructuredAPhiFields &solution_fields,
               std::string &solver_message);

    bool solveDiscreteCompatiblePhi(const StdVec<Real> &electrical_conductivity,
                                    const StdVec<Real> &magnetic_reluctivity,
                                    const LaplaceStructuredAPhiBoundaryCondition &boundary_condition,
                                    LaplaceStructuredAPhiFields &fields,
                                    std::string &solver_message);

    LaplaceStructuredAPhiDiagnostics computeDiagnostics(const StdVec<Real> &electrical_conductivity,
                                                        const StdVec<Real> &magnetic_reluctivity,
                                                        LaplaceStructuredAPhiFields &fields) const;

    VecC gatherUnknownVector(const LaplaceStructuredAPhiFields &fields) const;
    VecC gatherScaledUnknownVector(const LaplaceStructuredAPhiFields &fields) const;
    void scatterUnknownVector(const VecC &unknown_vector,
                              LaplaceStructuredAPhiFields &fields) const;

    size_t NumberOfParticles() const { return number_of_particles_; }
    size_t NumberOfUnknowns() const { return 4 * number_of_particles_; }

    int idxAx(size_t particle_index) const { return static_cast<int>(particle_index); }
    int idxAy(size_t particle_index) const { return static_cast<int>(number_of_particles_ + particle_index); }
    int idxAz(size_t particle_index) const { return static_cast<int>(2 * number_of_particles_ + particle_index); }
    int idxPhi(size_t particle_index) const { return static_cast<int>(3 * number_of_particles_ + particle_index); }

    const SpMatC &systemMatrix() const { return system_matrix_; }
    const VecC &rightHandSide() const { return right_hand_side_; }
    Real UsedPhiBlockScale() const { return used_phi_block_scale_; }
    int UsedDiagonalEquilibrationIterations() const { return used_diagonal_equilibration_iterations_; }
    const std::string &LastSolverBackend() const { return last_solver_backend_; }
    int LastSolverIterations() const { return last_solver_iterations_; }
    Real LastSolverEstimatedError() const { return last_solver_estimated_error_; }

  protected:
    void buildOperatorMatrices(const StdVec<Real> &electrical_conductivity,
                               const StdVec<Real> &magnetic_reluctivity);
    void buildSystemMatrix();
    void initializeBlockScaling();
    void applyBlockScalingToSystem();
    void applyDiagonalEquilibrationToSystem();
    void applyBoundaryCondition(const LaplaceStructuredAPhiBoundaryCondition &boundary_condition);
    VecC physicalToScaledUnknownVector(const VecC &physical_unknown_vector) const;
    VecC scaledToPhysicalUnknownVector(const VecC &scaled_unknown_vector) const;
    Complex scaleBoundaryUnknown(int index, const Complex &physical_value) const;

    Real pairwiseDiffusionWeight(const Neighborhood &inner_neighborhood,
                                 size_t index_i,
                                 size_t neighbor_index,
                                 Real reference_smoothing_length) const;

    static Real harmonicMean(Real value_i, Real value_j);
    static Real complexAbs(const Complex &value);
    static Real choosePhiBlockScale(Real angular_frequency, Real requested_scale);

    LaplaceStructuredAPhiParameters parameters_;
    size_t number_of_particles_;
    Real reference_smoothing_length_;

    Vecd *positions_;
    Real *volumetric_measure_;
    Real *smoothing_length_ratio_;
    const ParticleConfiguration *particle_configuration_;
    const StdVec<StdVec<uint8_t>> *neighbor_is_contact_;

    SpMatC grad_x_;
    SpMatC grad_y_;
    SpMatC grad_z_;
    SpMatC laplace_sigma_;
    SpMatC laplace_nu_;
    SpMatC system_matrix_;
    VecC right_hand_side_;
    StdVec<Real> left_block_scaling_;
    StdVec<Real> right_block_scaling_;
    Real used_phi_block_scale_;
    int used_diagonal_equilibration_iterations_;
    StdVec<Real> last_electrical_conductivity_;
    StdVec<Real> last_magnetic_reluctivity_;
    VecC rhs_unconstrained_physical_;
    std::string last_solver_backend_;
    int last_solver_iterations_;
    Real last_solver_estimated_error_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_LAPLACE_EIGEN_H
