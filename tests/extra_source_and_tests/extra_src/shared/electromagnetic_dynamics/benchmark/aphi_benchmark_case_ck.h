#ifndef APHI_BENCHMARK_CASE_CK_H
#define APHI_BENCHMARK_CASE_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace benchmark
{

struct AphiBoxRegion
{
    Real xmin = 0.0;
    Real xmax = 0.0;
    Real ymin = 0.0;
    Real ymax = 0.0;
    Real zmin = 0.0;
    Real zmax = 0.0;
};

struct AphiRegionMaterialProperties
{
    Real sigma = 0.0;
    Real nu = 1.0;
};

struct AphiTeam7LikeUnitBoxLayout
{
    AphiBoxRegion conductor{0.52, 0.68, 0.30, 0.70, 0.20, 0.80};
    AphiBoxRegion coil{0.24, 0.38, 0.15, 0.85, 0.20, 0.80};
    AphiRegionMaterialProperties air{1.0e-4, 1.0};
    AphiRegionMaterialProperties conductor_material{1.0, 1.0};
    AphiRegionMaterialProperties coil_material{1.0e-4, 1.0};
};

/** TEAM7 canonical physical dimensions [m] (fractions mapped from unit box). */
struct AphiTeam7PhysicalDimensions
{
    static constexpr Real length = 1.2;
    static constexpr Real height = 1.0;
    static constexpr Real width = 0.3;
};

inline AphiBoxRegion scaleBoxRegion(const AphiBoxRegion &region, Real body_length, Real body_height, Real body_width)
{
    return AphiBoxRegion{region.xmin * body_length, region.xmax * body_length, region.ymin * body_height,
                         region.ymax * body_height, region.zmin * body_width, region.zmax * body_width};
}

inline AphiTeam7LikeUnitBoxLayout buildTeam7LayoutForBox(Real body_length, Real body_height, Real body_width)
{
    const AphiTeam7LikeUnitBoxLayout unit_layout;
    AphiTeam7LikeUnitBoxLayout scaled_layout = unit_layout;
    scaled_layout.conductor = scaleBoxRegion(unit_layout.conductor, body_length, body_height, body_width);
    scaled_layout.coil = scaleBoxRegion(unit_layout.coil, body_length, body_height, body_width);
    return scaled_layout;
}

inline bool insideBoxRegion(const Vecd &position, const AphiBoxRegion &region)
{
    return position[0] >= region.xmin && position[0] <= region.xmax && position[1] >= region.ymin &&
           position[1] <= region.ymax && position[2] >= region.zmin && position[2] <= region.zmax;
}

inline Real smoothBoxProfile(const Vecd &position, const AphiBoxRegion &region, Real amplitude)
{
    if (!insideBoxRegion(position, region))
    {
        return 0.0;
    }
    const Real xi = (position[0] - region.xmin) / (region.xmax - region.xmin + TinyReal);
    const Real eta = (position[1] - region.ymin) / (region.ymax - region.ymin + TinyReal);
    const Real zeta = (position[2] - region.zmin) / (region.zmax - region.zmin + TinyReal);
    return amplitude * std::sin(Pi * xi) * std::sin(Pi * eta) * std::sin(Pi * zeta);
}

inline void smoothBoxProfileWithGradients(const Vecd &position, const AphiBoxRegion &region, Real amplitude, Real &g,
                                          Real &dgdx, Real &dgdy, Real &dgdz)
{
    g = 0.0;
    dgdx = 0.0;
    dgdy = 0.0;
    dgdz = 0.0;
    if (!insideBoxRegion(position, region))
    {
        return;
    }
    const Real length_x = region.xmax - region.xmin + TinyReal;
    const Real length_y = region.ymax - region.ymin + TinyReal;
    const Real length_z = region.zmax - region.zmin + TinyReal;
    const Real xi = (position[0] - region.xmin) / length_x;
    const Real eta = (position[1] - region.ymin) / length_y;
    const Real zeta = (position[2] - region.zmin) / length_z;
    const Real sin_xi = std::sin(Pi * xi);
    const Real sin_eta = std::sin(Pi * eta);
    const Real sin_zeta = std::sin(Pi * zeta);
    const Real cos_xi = std::cos(Pi * xi);
    const Real cos_eta = std::cos(Pi * eta);
    const Real cos_zeta = std::cos(Pi * zeta);
    g = amplitude * sin_xi * sin_eta * sin_zeta;
    dgdx = amplitude * (Pi / length_x) * cos_xi * sin_eta * sin_zeta;
    dgdy = amplitude * (Pi / length_y) * sin_xi * cos_eta * sin_zeta;
    dgdz = amplitude * (Pi / length_z) * sin_xi * sin_eta * cos_zeta;
}

/** sigma = sigma_left for x <= x_interface, else sigma_right. */
class AssignPiecewiseSigmaHalfSpaceCK : public LocalDynamics
{
  public:
    AssignPiecewiseSigmaHalfSpaceCK(SPHBody &sph_body, Real x_interface, Real sigma_left, Real sigma_right, Real nu,
                                    const AphiMaterialNames &material_names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Real *sigma_;
        Real *nu_;
        Real x_interface_;
        Real sigma_left_;
        Real sigma_right_;
        Real nu_value_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    Real x_interface_;
    Real sigma_left_;
    Real sigma_right_;
    Real nu_value_;
};

/** TEAM7-like conductor / coil / air material tagging on a unit box. */
class AssignTeam7LikeRegionMaterialsCK : public LocalDynamics
{
  public:
    AssignTeam7LikeRegionMaterialsCK(SPHBody &sph_body, const AphiTeam7LikeUnitBoxLayout &layout,
                                     const AphiMaterialNames &material_names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Real *sigma_;
        Real *nu_;
        AphiTeam7LikeUnitBoxLayout layout_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    AphiTeam7LikeUnitBoxLayout layout_;
};

/** Impressed current density on the A-equation RHS block (real/imag Vecd components). */
class AssignImpressedCurrentRhsCK : public LocalDynamics
{
  public:
    AssignImpressedCurrentRhsCK(SPHBody &sph_body, const AphiBlockNames &rhs_block, const AphiBoxRegion &source_region,
                                const Vecd &current_real, const Vecd &current_imag, Real peak_amplitude);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        AphiBoxRegion source_region_;
        Vecd current_real_;
        Vecd current_imag_;
        Real peak_amplitude_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    AphiBoxRegion source_region_;
    Vecd current_real_;
    Vecd current_imag_;
    Real peak_amplitude_;
};

/**
 * Solenoidal current density J = curl C on the A-equation RHS, with C = (0, 0, 0.5 g(x) (x^2 + y^2))
 * and g = smoothBoxProfile (zero with zero gradient outside the source box).
 */
class AssignSolenoidalCurlCurrentRhsCK : public LocalDynamics
{
  public:
    AssignSolenoidalCurlCurrentRhsCK(SPHBody &sph_body, const AphiBlockNames &rhs_block,
                                     const AphiBoxRegion &source_region, Real peak_amplitude,
                                     Real imag_to_real_ratio = 0.1);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        AphiBoxRegion source_region_;
        Real peak_amplitude_;
        Real imag_to_real_ratio_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    AphiBoxRegion source_region_;
    Real peak_amplitude_;
    Real imag_to_real_ratio_;
};

/**
 * Piecewise-linear phi_real with sigma_left * dphi/dx matched across x = x_interface.
 * A block is zero; phi_imag is zero. Used for interface-aware discrete MMS.
 */
class AssignInterfaceFluxMatchedPhiFieldsCK : public LocalDynamics
{
  public:
    AssignInterfaceFluxMatchedPhiFieldsCK(SPHBody &sph_body, Real x_interface, Real sigma_left, Real sigma_right,
                                        const AphiBlockNames &block_names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
        Real x_interface_;
        Real slope_left_;
        Real slope_right_;
        Real offset_right_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    Real x_interface_;
    Real slope_left_;
    Real slope_right_;
    Real offset_right_;
};

/** Cylindrical annulus source region in physical coordinates (xy plane, limited z slab). */
struct AphiAnnularSourceGeometry
{
    Real center_x = 0.372;
    Real center_y = 0.5;
    Real center_z = 0.15;
    Real inner_radius = 0.04;
    Real outer_radius = 0.09;
    Real z_half_height = 0.06;
};

inline bool insideAnnularSourceRegion(const Vecd &position, const AphiAnnularSourceGeometry &geometry)
{
    const Real dx = position[0] - geometry.center_x;
    const Real dy = position[1] - geometry.center_y;
    const Real radial_distance = std::sqrt(dx * dx + dy * dy);
    return radial_distance >= geometry.inner_radius && radial_distance <= geometry.outer_radius &&
           std::abs(position[2] - geometry.center_z) <= geometry.z_half_height;
}

/** Prescribed-current source: sigma = 0 inside annular source region. */
class AssignZeroSigmaInAnnularRegionCK : public LocalDynamics
{
  public:
    AssignZeroSigmaInAnnularRegionCK(SPHBody &sph_body, const AphiAnnularSourceGeometry &geometry,
                                     const AphiMaterialNames &material_names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Real *sigma_;
        AphiAnnularSourceGeometry geometry_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
    AphiAnnularSourceGeometry geometry_;
};

/** Impressed current density on the A-equation RHS inside an annular source region. */
class AssignImpressedCurrentRhsAnnularCK : public LocalDynamics
{
  public:
    AssignImpressedCurrentRhsAnnularCK(SPHBody &sph_body, const AphiBlockNames &rhs_block,
                                       const AphiAnnularSourceGeometry &geometry, const Vecd &current_real,
                                       const Vecd &current_imag, Real peak_amplitude);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        AphiAnnularSourceGeometry geometry_;
        Vecd current_real_;
        Vecd current_imag_;
        Real peak_amplitude_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    AphiAnnularSourceGeometry geometry_;
    Vecd current_real_;
    Vecd current_imag_;
    Real peak_amplitude_;
};

/** Uniform phi_real = -E0 * x, A = 0, phi_imag = 0 — analytic Joule verification (Stage 10A). */
class AssignUniformLinearPhiFieldsCK : public LocalDynamics
{
  public:
    AssignUniformLinearPhiFieldsCK(SPHBody &sph_body, Real electric_field_x, const AphiBlockNames &block_names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
        Real electric_field_x_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    Real electric_field_x_;
};

} // namespace benchmark
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BENCHMARK_CASE_CK_H
