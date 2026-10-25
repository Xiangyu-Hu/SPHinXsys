/** @file necking_sycl.cpp
 * @brief Three-dimensional Johnson-Cook necking with thermal softening and failure.
 */
#include "necking_dynamics.h"
#include "necking_regression.h"
#include "sphinxsys.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <tbb/global_control.h>
#include <type_traits>

using namespace SPH;

template <class ExecutionPolicy>
int RunNecking(const NeckingOptions &options)
{
    tbb::global_control thread_limit(tbb::global_control::max_allowed_parallelism, options.threads);
    SPHSystem system(BoundingBoxd(Vecd(-0.008, -0.008, -0.012), Vecd(0.008, 0.008, 0.012)), options.spacing);
    SolidBody rod(system, makeShared<NeckingRod>("NeckingRod"));
    rod.defineMatterMaterial<JohnsonCookSolid>(7830, 200e9, 0.29, SteelParameters());
    BaseParticles &particles = rod.generateParticles<BaseParticles, NeckingLattice>(options.spacing);
    Vecd *position = particles.getVariableDataByName<Vecd>("Position");
    Vecd *initial_position = particles.registerStateVariableData<Vecd>("InitialPosition");
    for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
        initial_position[i] = position[i];
    particles.addEvolvingVariable<Vecd>("InitialPosition");

    auto &inner_relation = system.addInnerRelation(rod, ConfigType::Lagrangian);
    ParticleMethodContainer<ExecutionPolicy> methods(ExecutionPolicy{});
    auto &cells = methods.addCellLinkedListDynamics(rod);
    auto &configuration = methods.addRelationDynamics(inner_relation);
    auto &correct_configuration = methods.template addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(inner_relation);
    auto &first_half = methods.template addInteractionDynamicsOneLevel<solid_dynamics::JohnsonCookIntegration1stHalfCK>(inner_relation);
    auto &second_half = methods.template addInteractionDynamicsOneLevel<solid_dynamics::JohnsonCookIntegration2ndHalfCK>(inner_relation);
    auto &time_step = methods.template addReduceDynamics<solid_dynamics::JohnsonCookAcousticTimeStepCK>(rod, Real(0.3));
    using StatusReduction = QuantityReduce<SPHBody, ReduceMax<int>, SimpleEvaluation<DirectValue<int>>>;
    auto &integration_status = methods.template addReduceDynamics<StatusReduction>(rod, "JCIntegrationStatus");
    GripMotion motion;
    auto &constrain_grips = methods.template addStateDynamics<PrescribedNeckingGripsCK>(rod, motion);
    auto &reactions = methods.template addReduceDynamics<NeckingGripReactionCK>(rod);
    NeckingRegression regression(rod, options.regression_test);
#if SPHINXSYS_USE_SYCL
    if constexpr (std::is_base_of<execution::DeviceExecution<>, ExecutionPolicy>::value)
        std::cout << "Device: " << execution::execution_instance.getQueue().get_device().get_info<sycl::info::device::name>() << std::endl;
    else
#endif
        std::cout << "Device: native host computing kernels" << std::endl;
    std::unique_ptr<BodyStatesRecordingToVtpCK<ExecutionPolicy>> body_states;
    if (options.state_recording)
    {
        body_states = std::make_unique<BodyStatesRecordingToVtpCK<ExecutionPolicy>>(system);
        for (const char *name : {"JCPlasticStrain", "JCTemperature", "JCDamage", "VonMisesStress"})
            body_states->template addToWrite<Real>(rod, name);
        body_states->template addToWrite<int>(rod, "MaterialActive");
        body_states->template addToWrite<Vecd>(rod, "InitialPosition");
        body_states->template addToWrite<Vecd>(rod, "Velocity");
    }
    cells.exec();
    configuration.exec();
    correct_configuration.exec();

    const size_t particle_count = particles.TotalRealParticles();
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Real *plastic_strain = particles.getVariableDataByName<Real>("JCPlasticStrain");
    const Real *temperature = particles.getVariableDataByName<Real>("JCTemperature");
    const Real *damage = particles.getVariableDataByName<Real>("JCDamage");
    const Real *plastic_work = particles.getVariableDataByName<Real>("PlasticWorkSpecific");
    const int *active = particles.getVariableDataByName<int>("MaterialActive");
    const Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    const Matd *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real external_work = 0, previous_power = 0;
    Real peak_force = 0, maximum_damage = 0;
    size_t steps = 0;
    const auto start = std::chrono::steady_clock::now();
    std::ofstream history("history.csv");
    history << "time_s,extension_m,reaction_N,reaction_left_N,max_peeq,max_temperature_K,failed_particles,kinetic_J,plastic_work_J,external_work_J,min_J\n"
            << std::setprecision(17);
    auto record = [&](const Vec2d &reaction)
    {
        for (const char *name : {"JCPlasticStrain", "JCTemperature", "JCDamage", "PlasticWorkSpecific"})
            particles.getVariableByName<Real>(name)->prepareForOutput(ExecutionPolicy{});
        for (const char *name : {"Position", "Velocity"})
            particles.getVariableByName<Vecd>(name)->prepareForOutput(ExecutionPolicy{});
        particles.getVariableByName<Matd>("DeformationGradient")->prepareForOutput(ExecutionPolicy{});
        particles.getVariableByName<int>("MaterialActive")->prepareForOutput(ExecutionPolicy{});
        Real maximum_strain = 0, maximum_temperature = 298, kinetic_energy = 0, plastic_energy = 0, minimum_J = 1;
        size_t failed_particles = 0;
        for (size_t i = 0; i != particle_count; ++i)
        {
            const Real J = deformation[i].determinant();
            if (!position[i].allFinite() || !velocity[i].allFinite() || !std::isfinite(temperature[i]) || (active[i] && J <= 0))
                throw std::domain_error("Invalid necking state");
            maximum_strain = std::max(maximum_strain, plastic_strain[i]);
            maximum_temperature = std::max(maximum_temperature, temperature[i]);
            maximum_damage = std::max(maximum_damage, damage[i]);
            minimum_J = std::min(minimum_J, J);
            kinetic_energy += Real(0.5) * mass[i] * velocity[i].squaredNorm();
            plastic_energy += mass[i] * plastic_work[i];
            failed_particles += active[i] == 0;
        }
        peak_force = std::max(peak_force, reaction[0]);
        history << physical_time << ',' << 2 * motion.displacement << ',' << reaction[0] << ',' << reaction[1] << ','
                << maximum_strain << ',' << maximum_temperature << ',' << failed_particles << ',' << kinetic_energy << ','
                << plastic_energy << ',' << external_work << ',' << minimum_J << '\n';
    };
    auto write_frame = [&](int frame)
    {
        record(reactions.exec(motion.acceleration));
        regression.Record(frame);
        history.flush();
        if (body_states)
            body_states->writeToFile(frame);
        std::cout << "frame=" << frame << " time=" << physical_time << " steps=" << steps << std::endl;
    };
    std::ofstream metadata("case.json");
    metadata << std::setprecision(17) << "{\"spacing_m\":" << options.spacing << ",\"end_time_s\":" << options.end_time
             << ",\"frames\":" << options.frames << ",\"particles\":" << particle_count << ",\"cfl\":0.3,\"real_bytes\":" << sizeof(Real) << "}";
    metadata.close();
    write_frame(0);
    for (int frame = 1; frame <= options.frames; ++frame)
    {
        const Real output_time = options.end_time * frame / options.frames;
        while (physical_time < output_time)
        {
            const Real dt = std::min(time_step.exec(), output_time - physical_time);
            if (!(std::isfinite(dt) && dt > 1e-14))
                throw std::domain_error("Necking time step collapsed");
            motion.SetTime(physical_time, false);
            constrain_grips.exec();
            first_half.exec(dt);
            if (integration_status.exec() != 0)
                throw std::domain_error("JC device integration rejected an invalid state");
            motion.SetTime(physical_time + dt / 2, false);
            constrain_grips.exec();
            second_half.exec(dt);
            physical_time += dt;
            ++steps;
            motion.SetTime(physical_time, true);
            constrain_grips.exec();
            const Vec2d reaction = reactions.exec(motion.acceleration);
            const Real power = (reaction[0] - reaction[1]) * motion.velocity;
            if (steps % 10 == 0)
                record(reaction);
            external_work += (power + previous_power) * dt / 2;
            previous_power = power;
        }
        write_frame(frame);
    }
    if (options.end_time >= 60e-6 && maximum_damage == 0)
        throw std::runtime_error("Necking case did not exercise dynamic failure");
    regression.TestResult();
    std::ofstream summary("summary.json");
    summary << std::setprecision(17) << "{\"completed\":true,\"time_s\":" << physical_time << ",\"steps\":" << steps
            << ",\"peak_force_N\":" << peak_force << ",\"maximum_damage\":" << maximum_damage
            << ",\"elapsed_s\":" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << "}\n";
    return 0;
}
int main(int argc, char **argv)
try
{
    bool use_host = false;
    const auto options = ParseNeckingOptions(argc, argv, &use_host);
    if (use_host)
        return RunNecking<execution::ParallelPolicy>(options);
    return RunNecking<MainExecutionPolicy>(options);
}
catch (const std::exception &exception)
{
    std::cerr << exception.what() << std::endl;
    return 1;
}
