#ifndef ELECTROMAGNETIC_OPHELIE_PARAMETERS_H
#define ELECTROMAGNETIC_OPHELIE_PARAMETERS_H

#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OpheliePhiSolverKind
{
    Jacobi,
    GMRES,
    PCG
};

struct OphelieParameters
{
    /** Glass melt box half-sizes (m). */
    Real glass_halfsize_x_ = 0.325;
    Real glass_halfsize_y_ = 0.25;
    Real glass_halfsize_z_ = 0.325;

    /** Coil shell box: outer half-sizes and inner cutout half-sizes (m). */
    Real coil_outer_halfsize_x_ = 0.475;
    Real coil_outer_halfsize_y_ = 0.275;
    Real coil_outer_halfsize_z_ = 0.475;
    Real coil_inner_halfsize_x_ = 0.38;
    Real coil_inner_halfsize_y_ = 0.26;
    Real coil_inner_halfsize_z_ = 0.38;

    Vecd glass_center_ = Vecd(0.6, 0.5, 0.5);
    Vecd coil_center_ = Vecd(0.6, 0.5, 0.5);

    Real current_amplitude_ = 1.0;
    Real number_of_turns_ = 8.0;
    /** If positive, overrides equivalentCurrentDensity() for coil J0. */
    Real coil_j0_override_ = 0.0;
    /** Native STL: assign J only on outer xy shell (fraction of max radius); rescales J0 via smaller A_cross. */
    bool coil_j_outer_shell_only_ = false;
    Real coil_j_outer_shell_radius_fraction_ = 0.85;
    Real coil_max_xy_radius_m_ = 0.0;
    Real frequency_ = 300.0e3;

    Real sigma_glass_ = 16.0;
    Real target_joule_power_ = 50.0e3;
    /** When false or target_joule_power_<=0, field_scale=power_scale=1 (required for TEAM7 reference compare). */
    bool enable_power_scaling_ = true;

    Real softening_length_ = 1.0e-6;
    Real mu0_ = 4.0 * Pi * 1.0e-7;

    bool enable_phi_correction_ = true;
    OpheliePhiSolverKind phi_solver_kind_ = OpheliePhiSolverKind::PCG;
    Real phi_gauge_penalty_ = 1.0;
    Real pair_weight_regularization_ = 0.01;
    Real phi_jacobi_relaxation_ = 0.85;
    size_t phi_jacobi_max_iterations_ = 4000;
    Real phi_jacobi_tolerance_ = 1.0e-4;
    UnsignedInt phi_gmres_restart_dimension_ = 40;
    UnsignedInt phi_gmres_max_outer_iterations_ = 80;
    Real phi_gmres_tolerance_ = 1.0e-4;
    size_t phi_pcg_max_iterations_ = 6000;
    Real phi_pcg_tolerance_ = 5.0e-4;

    bool enable_self_induction_ = false;
    size_t self_induction_max_iterations_ = 8;
    Real self_induction_j_tolerance_ = 0.05;
    Real self_induction_relaxation_factor_ = 0.15;

    Real omega() const { return 2.0 * Pi * frequency_; }

    Real coilCrossSectionArea() const
    {
        return (coil_outer_halfsize_x_ - coil_inner_halfsize_x_) * 2.0 *
               (coil_outer_halfsize_y_ - coil_inner_halfsize_y_) * 2.0;
    }

    Real equivalentCurrentDensity() const
    {
        const Real coil_height = 2.0 * coil_outer_halfsize_y_;
        const Real radial_gap = 2.0 * (coil_outer_halfsize_x_ - coil_inner_halfsize_x_);
        return number_of_turns_ * current_amplitude_ / (radial_gap * coil_height + TinyReal);
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PARAMETERS_H
