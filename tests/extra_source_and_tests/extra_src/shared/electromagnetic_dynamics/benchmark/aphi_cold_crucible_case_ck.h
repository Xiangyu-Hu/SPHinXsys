#ifndef APHI_COLD_CRUCIBLE_CASE_CK_H
#define APHI_COLD_CRUCIBLE_CASE_CK_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace benchmark
{

/**
 * Simplified cold-crucible layout on a unit box (box-approximated regions, no Contact).
 *
 * Schematic (x-y mid-plane):
 *   [coil_L]  [crucible shell + melt pool]  [coil_R]
 *
 * Material priority: melt > crucible_wall > coil > air.
 */
struct AphiColdCrucibleUnitBoxLayout
{
    AphiBoxRegion melt{0.38, 0.62, 0.38, 0.62, 0.25, 0.75};
    /** Crucible shell (wall only; melt box overrides interior). */
    AphiBoxRegion crucible_wall{0.32, 0.68, 0.32, 0.68, 0.20, 0.80};
    AphiBoxRegion coil_left{0.10, 0.26, 0.20, 0.80, 0.20, 0.80};
    AphiBoxRegion coil_right{0.74, 0.90, 0.20, 0.80, 0.20, 0.80};

    AphiRegionMaterialProperties air{1.0e-4, 1.0};
    /** High-conductivity wall (copper-like, scaled for numerical stability). */
    AphiRegionMaterialProperties crucible_wall_material{1.0e4, 1.0};
    /** Liquid melt (conductive, TEAM7-like magnitude). */
    AphiRegionMaterialProperties melt_material{1.0, 1.0};
    /** Coil form / air-equivalent sigma; current enters via impressed RHS. */
    AphiRegionMaterialProperties coil_material{1.0e-4, 1.0};
};

/** Default physical box for cold-crucible scaffold (fraction-mapped). */
struct AphiColdCruciblePhysicalDimensions
{
    static constexpr Real length = 1.0;
    static constexpr Real height = 1.0;
    static constexpr Real width = 1.0;
};

/** Canonical solver / source settings for Stage 10C scaffold. */
struct AphiColdCrucibleCaseSpec
{
    static constexpr Real body_length = AphiColdCruciblePhysicalDimensions::length;
    static constexpr Real body_height = AphiColdCruciblePhysicalDimensions::height;
    static constexpr Real body_width = AphiColdCruciblePhysicalDimensions::width;

    static constexpr Real dp = 0.1;
    static constexpr Real omega = 1.25;
    static constexpr Real phi_gauge_penalty = 100.0;
    static constexpr Real impressed_current_amplitude = 8.0;

    static constexpr Real tolerance = 5.0e-4;
    static constexpr UnsignedInt restart_dimension = 50;
    static constexpr UnsignedInt max_outer_iterations = 150;

    static constexpr Real min_melt_solution_block_max = 1.0e-4;
};

inline AphiColdCrucibleUnitBoxLayout buildColdCrucibleLayoutForBox(Real body_length, Real body_height, Real body_width)
{
    const AphiColdCrucibleUnitBoxLayout unit_layout;
    AphiColdCrucibleUnitBoxLayout scaled = unit_layout;
    scaled.melt = scaleBoxRegion(unit_layout.melt, body_length, body_height, body_width);
    scaled.crucible_wall = scaleBoxRegion(unit_layout.crucible_wall, body_length, body_height, body_width);
    scaled.coil_left = scaleBoxRegion(unit_layout.coil_left, body_length, body_height, body_width);
    scaled.coil_right = scaleBoxRegion(unit_layout.coil_right, body_length, body_height, body_width);
    return scaled;
}

inline bool insideAnyCoilRegion(const Vecd &position, const AphiColdCrucibleUnitBoxLayout &layout)
{
    return insideBoxRegion(position, layout.coil_left) || insideBoxRegion(position, layout.coil_right);
}

/** Melt / crucible shell / coil / air material tagging. */
class AssignColdCrucibleRegionMaterialsCK : public LocalDynamics
{
  public:
    AssignColdCrucibleRegionMaterialsCK(SPHBody &sph_body, const AphiColdCrucibleUnitBoxLayout &layout,
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
        AphiColdCrucibleUnitBoxLayout layout_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    AphiColdCrucibleUnitBoxLayout layout_;
};

/** Impressed current in left + right coil regions (additive smooth box profiles). */
class AssignColdCrucibleCoilSourceCK : public LocalDynamics
{
  public:
    AssignColdCrucibleCoilSourceCK(SPHBody &sph_body, const AphiBlockNames &rhs_block,
                                   const AphiColdCrucibleUnitBoxLayout &layout, const Vecd &current_real,
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
        AphiColdCrucibleUnitBoxLayout layout_;
        Vecd current_real_;
        Vecd current_imag_;
        Real peak_amplitude_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    AphiColdCrucibleUnitBoxLayout layout_;
    Vecd current_real_;
    Vecd current_imag_;
    Real peak_amplitude_;
};

} // namespace benchmark
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_COLD_CRUCIBLE_CASE_CK_H
