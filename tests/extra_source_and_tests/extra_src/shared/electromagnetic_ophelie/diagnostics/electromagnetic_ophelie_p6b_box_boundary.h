#ifndef ELECTROMAGNETIC_OPHELIE_P6B_BOX_BOUNDARY_H
#define ELECTROMAGNETIC_OPHELIE_P6B_BOX_BOUNDARY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux_boundary_closure.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "sphinxsys.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kOphelieBoxSignedDistance = "SignedDistance";
inline constexpr const char *kOphelieBoxNormalDirection = "NormalDirection";

/** Populate SignedDistance/NormalDirection from body shape (host CK + device sync). */
inline void initBoxEdgeFluxBoundaryNormalsFromShape(SPHSystem &sph_system, SolidBody &glass_body)
{
    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(glass_body).exec();
    BaseParticles &particles = glass_body.getBaseParticles();
    syncVariableToDevice<Real>(particles, kOphelieBoxSignedDistance);
    syncVariableToDevice<Vecd>(particles, kOphelieBoxNormalDirection);
}

inline bool ophelieParamsNeedBoxEdgeFluxNormals(const OphelieParameters &params)
{
    return ophelieEdgeReconBoundaryModeUsesProductionClosure(params.edge_recon_boundary_mode_);
}

inline OpheliePhiBoundaryGeometryContext makeAnalyticBoxPhiBoundaryGeometry(const Vecd &center, const Vecd &halfsize)
{
    OpheliePhiBoundaryGeometryContext geom;
    geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticBox;
    geom.box_center = center;
    geom.box_halfsize = halfsize;
    return geom;
}

inline const char *p6bBoundarySweepLabel(const OphelieParameters &params)
{
    if (params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
    {
        return "phi-neumann";
    }
    return ophelieEdgeReconBoundaryModeName(params.edge_recon_boundary_mode_);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_P6B_BOX_BOUNDARY_H
