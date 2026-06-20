#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_GLASS_MESH_RELAX_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_GLASS_MESH_RELAX_H

#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_progress.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iomanip>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Log multiloop filament coil layout (axis-aligned solenoid-like rings). */
inline void logFrenchReducedCoilGeometry(const OphelieFrenchReducedCaseParams &french)
{
    std::cout << "[ophelie] coil_geometry: multiloop horizontal circular filaments (Biot I·dl), axis=+z\n"
              << "[ophelie]   stack_center=" << french.coil.stack_center.transpose()
              << " loop_radius=" << french.coil.loop_radius << " m (glass_radius=" << french.glass_radius << " m)\n"
              << "[ophelie]   num_loops=" << french.coil.num_loops << " z_range=[" << french.coil.z_min << ","
              << french.coil.z_max << "] segments_per_loop=" << french.coil.segments_per_loop
              << " cell_centered=" << (french.coil.use_cell_centered_loops ? 1 : 0) << std::endl;
    for (size_t k = 0; k < french.coil.num_loops; ++k)
    {
        const Real z_loop = multiloopZPosition(k, french.coil);
        std::cout << "[ophelie]   loop[" << k << "] center=(" << french.coil.stack_center[0] << ","
                  << french.coil.stack_center[1] << "," << z_loop << ") R=" << french.coil.loop_radius
                  << " I_per_loop=" << french.coil.current_per_loop << " A" << std::endl;
    }
}

inline BoundingBoxd frenchGlassMeshRelaxDomainBounds(const OphelieFrenchReducedCaseParams &french, Real pad)
{
    const Real extent = french.glass_radius + pad;
    const Vecd lower(french.glass_center[0] - extent, french.glass_center[1] - extent,
                     french.glass_center[2] - french.glass_half_height - pad);
    const Vecd upper(french.glass_center[0] + extent, french.glass_center[1] + extent,
                     french.glass_center[2] + french.glass_half_height + pad);
    return BoundingBoxd(lower, upper);
}

#if SPHINXSYS_USE_SYCL
/**
 * Particle relaxation for French glass cylinder using TriangleMeshShapeCylinder + SYCL-CK pipeline
 * (same pattern as test_3d_particle_relaxation_single_resolution_sycl).
 */
inline void runFrenchGlassSphinxsysStyleRelax(SPHSystem &sph_system, SolidBody &glass_body,
                                              LevelSetShape &level_set_shape, size_t relaxation_steps = 1000,
                                              size_t vtp_every = 100)
{
    const bool write_vtp = vtp_every > 0;
    const bool saved_state_recording = sph_system.StateRecording();
    if (write_vtp)
    {
        sph_system.setStateRecording(true);
    }

    OphelieProgressLogger progress("relax:GlassBody-mesh");
    progress.log("SYCL mesh-cylinder relax n_particles=" +
                 std::to_string(glass_body.getBaseParticles().TotalRealParticles()) +
                 " steps=" + std::to_string(relaxation_steps) +
                 (write_vtp ? " vtp_every=" + std::to_string(vtp_every) : " vtp=off"));

    NearShapeSurface near_body_surface(glass_body);
    Inner<> body_inner(glass_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(glass_body);
    auto &update_inner_relation = main_methods.addRelationDynamics(body_inner);
    auto &random_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(glass_body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(glass_body, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(glass_body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(glass_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);

    BodyStatesRecordingToVtpCK<execution::MainExecutionPolicy> *body_state_recorder = nullptr;
    if (write_vtp)
    {
        body_state_recorder = &main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass_body);
    }

    random_particles.exec();
    progress.log("randomize done");

    if (write_vtp && body_state_recorder != nullptr)
    {
        glass_body.setNewlyUpdated();
        body_state_recorder->writeToFile(0);
        progress.log("VTP GlassBody_ite_0 -> " + IO::getEnvironment().OutputFolder() + "/");
    }

    for (size_t step = 0; step < relaxation_steps; ++step)
    {
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        relaxation_residual.exec();
        const Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();

        const size_t step_index = step + 1;
        if (step_index % 100 == 0 || step_index == relaxation_steps)
        {
            std::cout << std::fixed << std::setprecision(0) << "[ophelie] mesh relax step " << step_index << "/"
                      << relaxation_steps << std::endl;
        }

        if (write_vtp && body_state_recorder != nullptr &&
            (step_index % vtp_every == 0 || step_index == relaxation_steps))
        {
            glass_body.setNewlyUpdated();
            body_state_recorder->writeToFile(step_index);
            progress.log("VTP GlassBody_ite_" + std::to_string(step_index) + " -> " +
                         IO::getEnvironment().OutputFolder() + "/");
        }
    }

    glass_body.updateCellLinkedList();
    if (write_vtp)
    {
        sph_system.setStateRecording(saved_state_recording);
    }
    progress.finish();
}
#endif

inline LevelSetShape &defineFrenchGlassMeshLevelSet(SolidBody &body, bool write_level_set = false)
{
#if SPHINXSYS_USE_SYCL
    if (write_level_set)
    {
        return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().writeLevelSet();
    }
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
#else
    (void)write_level_set;
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
#endif
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_GLASS_MESH_RELAX_H
