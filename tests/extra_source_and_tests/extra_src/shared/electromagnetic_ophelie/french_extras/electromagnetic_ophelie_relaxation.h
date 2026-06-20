#ifndef ELECTROMAGNETIC_OPHELIE_RELAXATION_H
#define ELECTROMAGNETIC_OPHELIE_RELAXATION_H

#include "electromagnetic_ophelie_progress.h"
#include "sphinxsys.h"

#include "io_environment.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

namespace detail
{
inline void logRelaxVtpWrite(OphelieProgressLogger &progress, const std::string &body_name, size_t step_index)
{
    const std::string folder = IO::getEnvironment().OutputFolder();
    progress.log("VTP " + body_name + "_ite_" + std::to_string(step_index) + " -> " + folder + "/");
}
} // namespace detail

#if SPHINXSYS_USE_SYCL
/**
 * SYCL particle relaxation (same pipeline as
 * tests/tests_sycl/3d_examples/test_3d_particle_relaxation_single_resolution_sycl).
 * Do NOT use legacy RelaxationStepInner on device — it can hang on complex level-set shapes.
 */
inline void relaxSolidBodyParticles(SPHSystem &sph_system, SolidBody &body, LevelSetShape &level_set_shape,
                                    const std::string &body_label, size_t relaxation_steps = 400,
                                    size_t log_every = 50, size_t relax_vtp_every = 100)
{
    const bool write_relax_vtp = relax_vtp_every > 0;
    const size_t vtp_interval = relax_vtp_every;
    const bool saved_state_recording = sph_system.StateRecording();
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(true);
    }

    OphelieProgressLogger progress("relax:" + body_label);
    progress.log("SYCL-CK relax n_particles=" + std::to_string(body.getBaseParticles().TotalRealParticles()) +
                 " steps=" + std::to_string(relaxation_steps) +
                 (write_relax_vtp ? " vtp_every=" + std::to_string(vtp_interval) : " vtp=off"));

    NearShapeSurface near_body_surface(body);
    Inner<> body_inner(body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(body);
    auto &update_inner_relation = main_methods.addRelationDynamics(body_inner);
    auto &random_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(body, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);

    BodyStatesRecordingToVtpCK<MainExecutionPolicy> *body_state_recorder = nullptr;
    if (write_relax_vtp)
    {
        body_state_recorder = &main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(body);
    }

    random_particles.exec();
    progress.log("randomize done");

    if (write_relax_vtp && body_state_recorder != nullptr)
    {
        body.setNewlyUpdated();
        body_state_recorder->writeToFile(0);
        detail::logRelaxVtpWrite(progress, body.Name(), 0);
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
        if (step_index % log_every == 0 || step_index == relaxation_steps)
        {
            progress.log("step " + std::to_string(step_index) + "/" + std::to_string(relaxation_steps));
        }

        if (write_relax_vtp && body_state_recorder != nullptr &&
            (step_index % vtp_interval == 0 || step_index == relaxation_steps))
        {
            body.setNewlyUpdated();
            body_state_recorder->writeToFile(step_index);
            detail::logRelaxVtpWrite(progress, body.Name(), step_index);
        }
    }

    body.updateCellLinkedList();
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(saved_state_recording);
    }
    progress.finish();
}
#else
/** CPU legacy relax (RelaxationStepLevelSetCorrectionInner), cf. test_3d_particle_relaxation_single_resolution. */
inline void relaxSolidBodyParticles(SPHSystem &sph_system, SolidBody &body, LevelSetShape &level_set_shape,
                                    const std::string &body_label, size_t relaxation_steps = 400,
                                    size_t log_every = 50, size_t relax_vtp_every = 100)
{
    (void)level_set_shape;
    const bool write_relax_vtp = relax_vtp_every > 0;
    const size_t vtp_interval = relax_vtp_every;
    const bool saved_state_recording = sph_system.StateRecording();
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(true);
    }

    using namespace relax_dynamics;
    OphelieProgressLogger progress("relax:" + body_label);
    progress.log("CPU legacy relax n_particles=" + std::to_string(body.getBaseParticles().TotalRealParticles()) +
                 " steps=" + std::to_string(relaxation_steps) +
                 (write_relax_vtp ? " vtp_every=" + std::to_string(vtp_interval) : " vtp=off"));

    InnerRelation body_inner(body);
    SimpleDynamics<RandomizeParticlePosition> randomize_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(body_inner);
    BodyStatesRecordingToVtp write_body_to_vtp(body);

    randomize_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    progress.log("entering relax iterations");

    if (write_relax_vtp)
    {
        body.setNewlyUpdated();
        write_body_to_vtp.writeToFile(0);
        detail::logRelaxVtpWrite(progress, body.Name(), 0);
    }

    for (size_t step = 0; step < relaxation_steps; ++step)
    {
        relaxation_step_inner.exec();
        const size_t step_index = step + 1;
        if (step_index % log_every == 0 || step_index == relaxation_steps)
        {
            progress.log("step " + std::to_string(step_index) + "/" + std::to_string(relaxation_steps));
        }

        if (write_relax_vtp && (step_index % vtp_interval == 0 || step_index == relaxation_steps))
        {
            body.setNewlyUpdated();
            write_body_to_vtp.writeToFile(step_index);
            detail::logRelaxVtpWrite(progress, body.Name(), step_index);
        }
    }
    body.updateCellLinkedList();
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(saved_state_recording);
    }
    progress.finish();
}
#endif

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_RELAXATION_H
