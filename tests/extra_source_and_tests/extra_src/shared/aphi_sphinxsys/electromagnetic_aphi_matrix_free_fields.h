#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_FIELDS_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_FIELDS_H

#include "inner_body_relation.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_types.h"

namespace SPH
{
namespace electromagnetics
{
struct MatrixFreeAPhiParameters
{
    Real angular_frequency = 0.0;
    Real geom_length_to_m = 1.0;
    Real pair_weight_regularization = 0.01;
    Real interface_contrast_threshold = 10.0;
    Real contact_gradient_scale = 1.0;
    Real contact_diffusion_scale = 1.0;
    Real reference_pseudo_time_step = 1.0;
    Real relaxation_scaling = 1.0;
    Real max_change_rate = MaxReal;
    bool fix_phi_reference = true;
    size_t phi_reference_index = 0;
    Complex phi_reference_value = Complex(0.0, 0.0);
};

struct MatrixFreeAPhiFields
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

struct MatrixFreeAPhiDiagnostics
{
    Real residual_ax_l2 = 0.0;
    Real residual_ay_l2 = 0.0;
    Real residual_az_l2 = 0.0;
    Real residual_phi_l2 = 0.0;
    Real divergence_a_l2 = 0.0;
    Real current_continuity_l2 = 0.0;
    Real max_abs_a = 0.0;
    Real max_abs_e = 0.0;
    Real max_abs_j = 0.0;
    Real total_joule_power = 0.0;
};

struct MatrixFreeAPhiDiscreteView
{
    size_t number_of_particles = 0;
    Real reference_smoothing_length = 0.0;
    Vecd *positions = nullptr;
    Real *volumetric_measure = nullptr;
    Real *smoothing_length_ratio = nullptr;
    const ParticleConfiguration *particle_configuration = nullptr;
    const StdVec<StdVec<uint8_t>> *neighbor_is_contact = nullptr;
};
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_FIELDS_H
