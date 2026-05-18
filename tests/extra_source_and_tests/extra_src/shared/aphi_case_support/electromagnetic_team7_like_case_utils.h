#ifndef ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_H
#define ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_fields.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include "legacy_aphi_archive/baselines/electromagnetic_aphi_laplace_eigen.h"

namespace SPH
{
namespace electromagnetics
{
namespace team7_like
{

struct Team7LikeGeometry
{
    Real body_length = 1.20;
    Real body_height = 1.00;
    Real body_width = 0.30;
    Real boundary_width = 0.0;
};

struct Team7LikeRegionBoxes
{
    Real conductor_xmin = 0.0;
    Real conductor_xmax = 0.0;
    Real conductor_ymin = 0.0;
    Real conductor_ymax = 0.0;
    Real conductor_zmin = 0.0;
    Real conductor_zmax = 0.0;
    Real coil_xmin = 0.0;
    Real coil_xmax = 0.0;
    Real coil_ymin = 0.0;
    Real coil_ymax = 0.0;
    Real coil_zmin = 0.0;
    Real coil_zmax = 0.0;
};

struct Team7LikeMaterialProperties
{
    Real sigma_air = 1.0e-4;
    Real sigma_conductor = 1.0;
    Real sigma_coil = 1.0e-4;
    Real nu_air = 1.0;
    Real nu_conductor = 1.0;
    Real nu_coil = 1.0;
};

struct Team7LikeMaterialAssignment
{
    StdVec<Real> sigma;
    StdVec<Real> nu;
    StdVec<bool> is_conductor;
    StdVec<bool> is_coil;
    StdVec<bool> is_source;
    size_t conductor_particles = 0;
    size_t coil_particles = 0;
    size_t air_particles = 0;
    size_t source_particles = 0;
};

struct Team7LikeManufacturedSetup
{
    MatrixFreeAPhiFields exact_fields;
    matrix_free::MatrixFreeAPhiSources sources;
    bool has_discrete_manufactured_reference = false;
    Real reference_phi_build_residual_l2 = 0.0;
};

bool insideBoxRegion(const Vecd &position, Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax);
Real coordinateFromFraction(Real fraction, Real extent);
Real smoothBoxProfile(const Vecd &position, Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax, Real amplitude);
bool isExteriorBoundaryParticle(const Vecd &position, const Team7LikeGeometry &geometry, Real shell_thickness);

Team7LikeGeometry readTeam7LikeGeometryFromEnvironment(Real dp_0);
Team7LikeRegionBoxes readTeam7LikeRegionBoxesFromEnvironment(const Team7LikeGeometry &geometry);
Team7LikeMaterialProperties readTeam7LikeMaterialPropertiesFromEnvironment();

/** Maps compare-test tokens `source_free` / `driven` to matrix-free smoke case mode strings. */
std::string normalizeManufacturedCaseMode(const std::string &case_mode_token);

Team7LikeMaterialAssignment buildTeam7LikeMaterialAssignment(const Vecd *positions, size_t number_of_particles,
                                                             const Team7LikeGeometry &geometry,
                                                             const Team7LikeRegionBoxes &region_boxes,
                                                             const Team7LikeMaterialProperties &materials);

void enforcePhiReference(StdVec<Complex> &phi, size_t reference_index, Complex reference_value);

matrix_free::ScalarComplexHelmholtzSolverState solveDiscretePhiForSourceFreeCase(const matrix_free::MatrixFreePairwiseGraph &graph,
                                                                                 const StdVec<Real> &sigma,
                                                                                 const StdVec<Complex> &rhs,
                                                                                 StdVec<Complex> &phi);

/** Build manufactured exact fields and discrete RHS sources for TEAM7-like coulomb cases. */
bool buildTeam7LikeManufacturedSetup(const std::string &case_mode, const matrix_free::MatrixFreePairwiseGraph &graph,
                                     const Vecd *positions, size_t number_of_particles, const Team7LikeGeometry &geometry,
                                     const Team7LikeRegionBoxes &region_boxes,
                                     const Team7LikeMaterialAssignment &material_assignment,
                                     const MatrixFreeAPhiParameters &parameters,
                                     const StdVec<Real> &sigma, const StdVec<Real> &nu, Team7LikeManufacturedSetup &setup,
                                     std::string &error_message);

LaplaceStructuredAPhiBoundaryCondition buildSparseBaselineBoundaryCondition(
    const Vecd *positions, size_t number_of_particles, const Team7LikeGeometry &geometry, Real boundary_shell_thickness,
    const MatrixFreeAPhiFields &exact_fields, size_t &phi_reference_index_out);

} // namespace team7_like
} // namespace electromagnetics
} // namespace SPH

#include "aphi_case_support/electromagnetic_team7_like_case_utils.hpp"

#endif // ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_H
