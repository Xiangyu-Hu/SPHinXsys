// SPDX-License-Identifier: Apache-2.0
/** HJC concrete cylinder impacting a fixed wall, following the Taylor-bar setup. */
#include "sphinxsys.h"
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

using namespace SPH;

int main(int argc, char *argv[])
{
    Real spacing = .0005, speed = 30, cfl = .2, end_time = 6e-5;
    bool regression_test = false;
    std::vector<char *> system_args{argv[0]};
    for (int i = 1; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "--regression-test")
            regression_test = true;
        else if (arg.rfind("--spacing=", 0) == 0)
            spacing = std::stod(arg.substr(10));
        else if (arg.rfind("--speed=", 0) == 0)
            speed = std::stod(arg.substr(8));
        else if (arg.rfind("--cfl=", 0) == 0)
            cfl = std::stod(arg.substr(6));
        else if (arg.rfind("--end-time=", 0) == 0)
            end_time = std::stod(arg.substr(11));
        else
            system_args.push_back(argv[i]);
    }
    if (!(spacing > 0 && spacing <= Real(.001) && speed > 0 && std::isfinite(speed) &&
          cfl > 0 && cfl <= Real(.2) && end_time > 0 && std::isfinite(end_time)))
        throw std::invalid_argument("Require 0 < spacing <= 0.001 m, speed > 0, 0 < CFL <= 0.2 and end-time > 0");
    if (regression_test && (spacing != Real(.001) || speed != Real(30) ||
                            cfl != Real(.2) || end_time != Real(6e-5)))
        throw std::invalid_argument("Regression test requires spacing=0.001, speed=30, cfl=0.2 and end-time=0.00006");
    const Real radius = .005, length = .020, wall_depth = 4 * spacing;
    // Keep the coarse lattice off the cylinder end faces.
    const Real gap = std::max(Real(.0005), spacing);
    SPHSystem system(BoundingBoxd(Vecd(-.02, -.02, -wall_depth), Vecd(.02, .02, .04)), spacing);
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(int(system_args.size()), system_args.data());
#endif
    auto shape = makeShared<ComplexShape>("Concrete");
    shape->add<TriangleMeshShapeCylinder>(Vec3d::UnitZ(), radius, length / 2, 48,
                                          Vecd(0, 0, length / 2 + gap));
    SolidBody column(system, shape);
    HJCParameters parameters{2.417e10, .29, 2.06, .0013, .866, 1.19e8, 8.2e6, 1, .01, 5,
                             4e7, .00124, 1.2e9, .011, .04, 1, 1.287e10, 1.631e10, 6.495e10};
    auto &material = column.defineMatterMaterial<HJCSolid>(2700, parameters);
    auto &particles = column.generateParticles<BaseParticles, Lattice>();
    auto wall_shape = makeShared<ComplexShape>("Wall");
    wall_shape->add<TriangleMeshShapeBrick>(Vecd(.012, .012, wall_depth / 2), 4, Vecd(0, 0, -wall_depth / 2));
    SolidBody wall(system, wall_shape);
    wall.defineMatterMaterial<SaintVenantKirchhoffSolid>(2700, material.YoungsModulus(), material.PoissonRatio());
    wall.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner(column);
    SurfaceContactRelation contact(column, {&wall});
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correction(inner);
    Dynamics1Level<solid_dynamics::HJCIntegration1stHalf> first_half(inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> second_half(inner);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> contact_factor(contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> contact_force(contact);
    ReduceDynamics<solid_dynamics::HJCAcousticTimeStep> time_step(column, cfl);
    BodyStatesRecordingToVtp states(system);
    for (const char *name : {"HJCDamage", "Pressure", "VonMisesStress", "HJCPlasticStrain", "HJCPlasticVolume"})
        states.addToWrite<Real>(column, name);

    using EnergyRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<TotalKineticEnergy>>;
    using DamageRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<Average<QuantitySummation<Real>>>>;
    std::unique_ptr<EnergyRegression> energy_regression;
    std::unique_ptr<DamageRegression> damage_regression;
    if (regression_test)
    {
        energy_regression = std::make_unique<EnergyRegression>(column);
        damage_regression = std::make_unique<DamageRegression>(column, "HJCDamage");
    }

    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    correction.exec();
    auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    auto *mass = particles.getVariableDataByName<Real>("Mass");
    auto *damage = particles.getVariableDataByName<Real>("HJCDamage");
    auto *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
    auto *force = particles.getVariableDataByName<Vecd>("RepulsionForce");
    auto *pressure = particles.getVariableDataByName<Real>("Pressure");
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        velocity[i][2] = -speed;
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    std::ofstream metadata("case.json");
    metadata << std::setprecision(16) << "{\"spacing\":" << spacing << ",\"cfl\":" << cfl
             << ",\"speed\":" << speed << ",\"end_time\":" << end_time
             << ",\"initial_gap\":" << gap << ",\"particles\":" << particles.TotalRealParticles() << "}\n";
    std::ofstream history("history.csv");
    history << "time,force_z,mean_velocity_z,mean_damage,max_damage,min_J,kinetic_energy,max_pressure\n" << std::setprecision(16);
    Real last_mean_damage = 0, peak_force = 0;
    auto record = [&]()
    {
        Real total_mass = 0, vz = 0, d = 0, max_d = 0, min_J = 1, energy = 0, fz = 0, max_p = 0;
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        {
            Real J = deformation[i].determinant();
            if (!velocity[i].allFinite() || !std::isfinite(J) || J <= 0 ||
                !std::isfinite(damage[i]) || damage[i] < 0 || damage[i] > 1)
                throw std::runtime_error("Nonphysical HJC impact state");
            total_mass += mass[i];
            vz += mass[i] * velocity[i][2];
            d += mass[i] * damage[i];
            max_d = std::max(max_d, damage[i]);
            min_J = std::min(min_J, J);
            energy += .5 * mass[i] * velocity[i].squaredNorm();
            fz += force[i][2];
            max_p = std::max(max_p, pressure[i]);
        }
        d /= total_mass;
        if (d + 1e-12 < last_mean_damage)
            throw std::runtime_error("HJC damage must not heal");
        last_mean_damage = d;
        peak_force = std::max(peak_force, fz);
        history << time << ',' << fz << ',' << vz / total_mass << ',' << d << ',' << max_d << ',' << min_J << ',' << energy << ',' << max_p << '\n';
    };
    record();
    states.writeToFile(0);
    if (regression_test)
    {
        energy_regression->writeToFile(0);
        damage_regression->writeToFile(0);
    }
    size_t step = 0;
    for (int frame = 1; frame <= 30; ++frame)
    {
        Real next_output = end_time * frame / 30;
        while (time < next_output)
        {
            contact_factor.exec();
            contact_force.exec();
            Real dt = std::min(time_step.exec(), next_output - time);
            if (!(dt > 1e-12 && std::isfinite(dt)))
                throw std::runtime_error("HJC impact timestep collapsed");
            first_half.exec(dt);
            second_half.exec(dt);
            time += dt;
            ++step;
            column.updateCellLinkedList();
            contact.updateConfiguration();
            record();
        }
        states.writeToFile(frame);
        if (regression_test)
        {
            energy_regression->writeToFile(frame);
            damage_regression->writeToFile(frame);
        }
        history.flush();
        std::cout << "t=" << time << " steps=" << step << " mean damage=" << last_mean_damage << '\n';
    }
    if (end_time >= 2e-5 && (peak_force <= 0 || last_mean_damage <= 0))
        throw std::runtime_error("Impact did not exercise contact and HJC damage");
    if (regression_test)
    {
        if (system.GenerateRegressionData())
        {
            energy_regression->generateDataBase(1e-3, 1e-3);
            damage_regression->generateDataBase(1e-3, 1e-3);
        }
        else
        {
            energy_regression->testResult();
            damage_regression->testResult();
        }
    }
    return 0;
}
