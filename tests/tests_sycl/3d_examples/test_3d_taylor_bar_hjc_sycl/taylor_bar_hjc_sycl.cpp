// SPDX-License-Identifier: Apache-2.0
/** HJC Taylor-bar impact using the shared CPU/SYCL computing kernels. */
#include "sphinxsys.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

using namespace SPH;

int main(int argc, char *argv[])
{
#if SPHINXSYS_USE_SYCL
    std::cout << "Device: " << execution::execution_instance.getQueue().get_device()
        .get_info<sycl::info::device::name>() << '\n';
#endif
    Real spacing = Real(.0005), speed = 30, cfl = Real(.2), end_time = Real(6e-5);
    bool regression_test = false;
    std::vector<char *> system_args{argv[0]};
    for (int i = 1; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "--regression-test")
            regression_test = true;
        else if (arg.rfind("--spacing=", 0) == 0)
            spacing = Real(std::stod(arg.substr(10)));
        else if (arg.rfind("--speed=", 0) == 0)
            speed = Real(std::stod(arg.substr(8)));
        else if (arg.rfind("--cfl=", 0) == 0)
            cfl = Real(std::stod(arg.substr(6)));
        else if (arg.rfind("--end-time=", 0) == 0)
            end_time = Real(std::stod(arg.substr(11)));
        else
            system_args.push_back(argv[i]);
    }
    if (!(spacing > 0 && spacing <= Real(.001) && speed > 0 && std::isfinite(speed) &&
          cfl > 0 && cfl <= Real(.2) && end_time > 0 && std::isfinite(end_time)))
        throw std::invalid_argument("Require 0 < spacing <= 0.001 m, speed > 0, 0 < CFL <= 0.2 and end-time > 0");
    if (regression_test && (spacing != Real(.001) || speed != Real(30) ||
                            cfl != Real(.2) || end_time != Real(6e-5)))
        throw std::invalid_argument("Regression requires spacing=0.001, speed=30, cfl=0.2 and end-time=0.00006");

    const Real radius = Real(.005), length = Real(.020), wall_depth = 4 * spacing;
    const Real gap = std::max(Real(.0005), spacing);
    SPHSystem system(BoundingBoxd(Vecd(-.02, -.02, -wall_depth), Vecd(.02, .02, .04)), spacing);
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(int(system_args.size()), system_args.data());
#endif
    auto &column_shape = system.addShape<TriangleMeshShapeCylinder>(
        Vec3d::UnitZ(), radius, length / 2, 48, Vecd(0, 0, length / 2 + gap), "Concrete");
    auto &wall_shape = system.addShape<TriangleMeshShapeBrick>(
        Vecd(.012, .012, wall_depth / 2), 4, Vecd(0, 0, -wall_depth / 2), "Wall");
    auto &column = system.addBody<SolidBody>(column_shape);
    HJCParameters parameters{Real(2.417e10), Real(.29), Real(2.06), Real(.0013), Real(.866),
        Real(1.19e8), Real(8.2e6), 1, Real(.01), 5, Real(4e7), Real(.00124), Real(1.2e9),
        Real(.011), Real(.04), 1, Real(1.287e10), Real(1.631e10), Real(6.495e10)};
    column.defineMatterMaterial<HJCSolid>(2700, parameters);
    auto &particles = column.generateParticles<BaseParticles, Lattice>();
    auto &wall = system.addBody<SolidBody>(wall_shape);
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();

    auto &inner = system.addInnerRelation(column, ConfigType::Lagrangian);
    auto &contact = system.addContactRelation(column, wall);
    SPHSolver solver(system);
    auto &host = solver.getHostMethodContainer();
    auto &methods = solver.getMainMethodContainer();
    auto &wall_cells = methods.addCellLinkedListDynamics(wall);
    auto &column_cells = methods.addCellLinkedListDynamics(column);
    auto &inner_configuration = methods.addRelationDynamics(inner);
    auto &contact_configuration = methods.addRelationDynamics(contact);
    host.addStateDynamics<NormalFromBodyShapeCK>(wall).exec();
    auto &correction = methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(inner);
    auto &first_half = methods.addInteractionDynamicsOneLevel<solid_dynamics::HJCIntegration1stHalfCK>(inner);
    auto &second_half = methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(inner);
    auto &contact_factor = methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(contact);
    auto &contact_force = methods.addInteractionDynamicsWithUpdate<solid_dynamics::RepulsionForceCK, Wall>(contact);
    auto &time_step = methods.addReduceDynamics<solid_dynamics::HJCAcousticTimeStepCK>(column, cfl);
    using StatusReduction = QuantityReduce<SPHBody, ReduceMax<int>, SimpleEvaluation<DirectValue<int>>>;
    auto &integration_status = methods.addReduceDynamics<StatusReduction>(column, "HJCIntegrationStatus");
    host.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(column, "Velocity", Vecd(0, 0, -speed)).exec();

    auto &states = methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(system);
    for (const char *name : {"HJCDamage", "Pressure", "VonMisesStress", "HJCPlasticStrain", "HJCPlasticVolume"})
        states.addToWrite<Real>(column, name);
    using EnergyRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<MainExecutionPolicy, TotalKineticEnergyCK>>;
    using DamageRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<MainExecutionPolicy, QuantityAverage<Real>>>;
    EnergyRegression *energy_regression = nullptr;
    DamageRegression *damage_regression = nullptr;
    if (regression_test)
    {
        energy_regression = &methods.addReduceRegression<RegressionTestEnsembleAverage, TotalKineticEnergyCK>(column);
        damage_regression = &methods.addReduceRegression<RegressionTestEnsembleAverage, QuantityAverage<Real>>(column, "HJCDamage");
    }

    wall_cells.exec();
    column_cells.exec();
    inner_configuration.exec();
    contact_configuration.exec();
    correction.exec();
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    std::ofstream metadata("case.json");
    metadata.exceptions(std::ios::failbit | std::ios::badbit);
    metadata << std::setprecision(16) << "{\"spacing\":" << spacing << ",\"cfl\":" << cfl
             << ",\"speed\":" << speed << ",\"end_time\":" << end_time << ",\"initial_gap\":" << gap
             << ",\"particles\":" << particles.TotalRealParticles() << ",\"real_bytes\":" << sizeof(Real) << "}\n";
    metadata.close();
    std::ofstream history("history.csv");
    history.exceptions(std::ios::failbit | std::ios::badbit);
    history << "time,force_z,mean_velocity_z,mean_damage,max_damage,min_J,kinetic_energy,max_pressure\n" << std::setprecision(16);
    double last_mean_damage = 0, peak_force = 0;
    auto record = [&]()
    {
        for (const char *name : {"Velocity", "RepulsionForce"})
        {
            std::string variable = std::string(name) == "RepulsionForce" ? "RepulsionForce" + contact.Name() : name;
            particles.getVariableByName<Vecd>(variable)->prepareForOutput(par_ck);
        }
        for (const char *name : {"HJCDamage", "Pressure"})
            particles.getVariableByName<Real>(name)->prepareForOutput(par_ck);
        particles.getVariableByName<Matd>("DeformationGradient")->prepareForOutput(par_ck);
        auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
        auto *force = particles.getVariableDataByName<Vecd>("RepulsionForce" + contact.Name());
        auto *mass = particles.getVariableDataByName<Real>("Mass");
        auto *damage = particles.getVariableDataByName<Real>("HJCDamage");
        auto *pressure = particles.getVariableDataByName<Real>("Pressure");
        auto *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
        double total_mass = 0, vz = 0, d = 0, max_d = 0, min_J = 1, energy = 0, fz = 0, max_p = 0;
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        {
            double J = deformation[i].determinant();
            if (!velocity[i].allFinite() || !force[i].allFinite() || !std::isfinite(J) || J <= 0 ||
                !std::isfinite(damage[i]) || damage[i] < 0 || damage[i] > 1 || !std::isfinite(pressure[i]))
                throw std::runtime_error("Nonphysical HJC impact state");
            total_mass += mass[i];
            vz += double(mass[i]) * velocity[i][2];
            d += double(mass[i]) * damage[i];
            max_d = std::max(max_d, double(damage[i]));
            min_J = std::min(min_J, J);
            energy += .5 * double(mass[i]) * velocity[i].squaredNorm();
            fz += force[i][2];
            max_p = std::max(max_p, double(pressure[i]));
        }
        d /= total_mass;
        if (d + 1e-12 < last_mean_damage)
            throw std::runtime_error("HJC damage must not heal");
        last_mean_damage = d;
        peak_force = std::max(peak_force, fz);
        history << time << ',' << fz << ',' << vz / total_mass << ',' << d << ',' << max_d << ',' << min_J << ',' << energy << ',' << max_p << '\n';
    };
    auto output = [&](int frame)
    {
        record();
        states.writeToFile(frame);
        if (regression_test)
        {
            energy_regression->writeToFile(frame);
            damage_regression->writeToFile(frame);
        }
    };
    output(0);
    size_t step = 0;
    TimeInterval compute_time;
    for (int frame = 1; frame <= 30; ++frame)
    {
        Real next_output = end_time * Real(frame) / Real(30);
        while (time < next_output)
        {
            TickCount start = TickCount::now();
            contact_factor.exec();
            contact_force.exec();
            Real dt = std::min(time_step.exec(), next_output - time);
            if (!(dt > 0 && std::isfinite(dt)) || time + dt <= time)
                throw std::runtime_error("HJC impact timestep collapsed");
            first_half.exec(dt);
            int status = integration_status.exec();
            if (status != 0)
                throw std::runtime_error("HJC constitutive update failed, status=" + std::to_string(status));
            second_half.exec(dt);
            time += dt;
            ++step;
            column_cells.exec();
            contact_configuration.exec();
            compute_time += TickCount::now() - start;
        }
        output(frame);
        history.flush();
        std::cout << "t=" << time << " steps=" << step << " mean damage=" << last_mean_damage << '\n';
    }
    std::cout << "Computation time: " << compute_time.seconds() << " s\n";
    if (regression_test && (peak_force <= 0 || last_mean_damage <= 0))
        throw std::runtime_error("Impact did not exercise contact and HJC damage");
    if (regression_test)
    {
        if (system.GenerateRegressionData())
        {
            energy_regression->generateDataBase(Real(1e-3), Real(1e-3));
            damage_regression->generateDataBase(Real(1e-3), Real(1e-3));
        }
        else
        {
            energy_regression->testResult();
            damage_regression->testResult();
        }
    }
    return 0;
}
