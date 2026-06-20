#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_NORMAL_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_NORMAL_H

#include "sphinxsys.h"

#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline constexpr const char *kTeam7NormalDirection = "NormalDirection";
inline constexpr const char *kTeam7SignedDistance = "SignedDistance";
inline constexpr const char *kTeam7EEdgeTangentLsDiag = "Team7EEdgeTangentLsDiag";
inline constexpr const char *kTeam7JEdgeTangentLsDiag = "Team7JEdgeTangentLsDiag";
inline constexpr const char *kTeam7EdgeTangentLsFallback = "Team7EdgeTangentLsFallback";

/** Run SPHinXsys NormalFromBodyShapeCK on one solid body (host exec). */
inline void computeTeam7SolidBodyNormalsFromShape(SPHSystem &sph_system, SolidBody &body)
{
    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(body).exec();
}

inline void computeTeam7CoilPlateNormalsFromShape(SPHSystem &sph_system, SolidBody &coil_body, SolidBody &plate_body)
{
    computeTeam7SolidBodyNormalsFromShape(sph_system, coil_body);
    computeTeam7SolidBodyNormalsFromShape(sph_system, plate_body);
}

inline bool team7ParticlesHaveSphinxNormals(BaseParticles &particles)
{
    return particles.getVariableDataByName<Vecd>(kTeam7NormalDirection) != nullptr &&
           particles.getVariableDataByName<Real>(kTeam7SignedDistance) != nullptr;
}

inline bool team7CoilPlateHaveSphinxNormals(SolidBody &coil_body, SolidBody &plate_body)
{
    return team7ParticlesHaveSphinxNormals(coil_body.getBaseParticles()) &&
           team7ParticlesHaveSphinxNormals(plate_body.getBaseParticles());
}

inline void registerTeam7CoilPlateNormalsForReload(ReloadParticleIO &reload_io, SolidBody &coil_body,
                                                   SolidBody &plate_body)
{
    reload_io.addToReload<Vecd>(coil_body, kTeam7NormalDirection);
    reload_io.addToReload<Real>(coil_body, kTeam7SignedDistance);
    reload_io.addToReload<Vecd>(plate_body, kTeam7NormalDirection);
    reload_io.addToReload<Real>(plate_body, kTeam7SignedDistance);
}

inline void reloadTeam7CoilPlateParticlesWithNormals(SolidBody &coil_body, SolidBody &plate_body)
{
    coil_body.generateParticles<BaseParticles, Reload>("CoilSourceBody")
        .reloadExtraVariable<Vecd>(kTeam7NormalDirection)
        .reloadExtraVariable<Real>(kTeam7SignedDistance);
    plate_body.generateParticles<BaseParticles, Reload>("PlateBody")
        .reloadExtraVariable<Vecd>(kTeam7NormalDirection)
        .reloadExtraVariable<Real>(kTeam7SignedDistance);
}

inline Real defaultTeam7BoundaryNormalWidth(Real dp) { return 2.0 * dp; }

inline bool isTeam7NearSurfaceParticle(Real signed_distance, Real boundary_width)
{
    return std::abs(signed_distance) <= boundary_width;
}

/**
 * Read SPHinXsys NormalDirection for particle i when |SignedDistance| <= boundary_width.
 * Returns false for interior particles or when reload fields are missing.
 */
inline bool getTeam7BoundaryNormal(BaseParticles &particles, size_t particle_index, Real boundary_width, Vecd &normal_out)
{
    const Real *signed_distance = particles.getVariableDataByName<Real>(kTeam7SignedDistance);
    const Vecd *normal_direction = particles.getVariableDataByName<Vecd>(kTeam7NormalDirection);
    if (signed_distance == nullptr || normal_direction == nullptr)
    {
        return false;
    }
    if (!isTeam7NearSurfaceParticle(signed_distance[particle_index], boundary_width))
    {
        return false;
    }
    const Real normal_norm = normal_direction[particle_index].norm();
    if (normal_norm <= TinyReal)
    {
        return false;
    }
    normal_out = normal_direction[particle_index] / normal_norm;
    return true;
}

inline void writeTeam7BoundaryNormalAuditCsv(SolidBody &coil_body, SolidBody &plate_body, Real dp,
                                             const std::string &csv_path, Real boundary_width_factor = 2.0)
{
    const Real boundary_width = boundary_width_factor * dp;
    std::ofstream csv(csv_path);
    if (!csv)
    {
        std::cerr << "[team7] failed to open boundary normal audit csv: " << csv_path << std::endl;
        return;
    }
    csv << "particle_id,body,region,x_m,y_m,z_m,normal_x,normal_y,normal_z,boundary_flag,signed_distance_m\n";

    const auto write_body = [&](SolidBody &body, const std::string &body_name)
    {
        BaseParticles &particles = body.getBaseParticles();
        if (!team7ParticlesHaveSphinxNormals(particles))
        {
            std::cerr << "[team7] boundary normal audit skipped for " << body_name
                      << ": missing NormalDirection/SignedDistance (rerun --relax=1)" << std::endl;
            return;
        }
        const size_t n = particles.TotalRealParticles();
        const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
        const Vecd *normal_direction = particles.getVariableDataByName<Vecd>(kTeam7NormalDirection);
        const Real *signed_distance = particles.getVariableDataByName<Real>(kTeam7SignedDistance);
        for (size_t i = 0; i < n; ++i)
        {
            const bool boundary_flag = isTeam7NearSurfaceParticle(signed_distance[i], boundary_width);
            Vecd unit_normal = Vecd::Zero();
            if (boundary_flag && normal_direction[i].norm() > TinyReal)
            {
                unit_normal = normal_direction[i].normalized();
            }
            csv << i << "," << body_name << "," << body_name << "," << pos[i][0] << "," << pos[i][1] << ","
                << pos[i][2] << "," << unit_normal[0] << "," << unit_normal[1] << "," << unit_normal[2] << ","
                << (boundary_flag ? 1 : 0) << "," << signed_distance[i] << "\n";
        }
    };

    write_body(coil_body, "CoilSourceBody");
    write_body(plate_body, "PlateBody");
    std::cout << "[team7] boundary normal audit -> " << csv_path << " boundary_width_m=" << boundary_width
              << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_NORMAL_H
