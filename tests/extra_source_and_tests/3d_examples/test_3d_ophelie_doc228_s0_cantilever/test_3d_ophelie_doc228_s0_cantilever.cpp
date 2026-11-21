/**
 * @file test_3d_ophelie_doc228_s0_cantilever.cpp
 * @brief Doc228 S0 — slender cantilever tip force vs δ=F L³/(3 E I) (≤10%).
 *
 * Load path: solid_dynamics::DistributingPointForces (same family as test_2d_plate).
 * Soft-E analytic gate. Stainless 304 preset written for paddle E0/E1.
 */
#include "electromagnetic_ophelie_client_stainless_steel.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{
Real get_physical_viscosity_general(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 0.4)
{
    return shape_constant / Real(4) * std::sqrt(rho * youngs_modulus) * length_scale;
}

struct S0Params
{
    Real length = Real(0.1);
    Real height = Real(0.0125);
    Real width = Real(0.0125);
    Real rho = Real(1265);
    Real youngs_modulus = Real(5.0e5);
    Real poisson = Real(0.30);
    Real tip_force_n = Real(0); // 0 → choose for δ_target ≈ 0.02 L
    Real end_time = Real(1.2);
    Real force_ramp_end = Real(0.2);
    Real rel_err_max = Real(0.10);
    std::string summary_json = "doc228_s0_summary.json";
    std::string out_dir = ".";
};

class FixPart : public BodyPartByParticle
{
  public:
    explicit FixPart(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&FixPart::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    }

  private:
    bool tagManually(size_t index_i) { return pos_[index_i].x() < Real(0); }
};

void applyCli(int ac, char *av[], S0Params &p)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--end-time=", 11) == 0)
            p.end_time = static_cast<Real>(std::atof(av[i] + 11));
        else if (std::strncmp(av[i], "--tip-force=", 12) == 0)
            p.tip_force_n = static_cast<Real>(std::atof(av[i] + 12));
        else if (std::strncmp(av[i], "--summary-json=", 15) == 0)
            p.summary_json = av[i] + 15;
        else if (std::strncmp(av[i], "--out-dir=", 10) == 0)
            p.out_dir = av[i] + 10;
        else if (std::strncmp(av[i], "--rel-err-max=", 14) == 0)
            p.rel_err_max = static_cast<Real>(std::atof(av[i] + 14));
        else if (std::strcmp(av[i], "--help") == 0 || std::strcmp(av[i], "-h") == 0)
        {
            std::cout << "Doc228 S0 tip-force cantilever: --end-time= --tip-force= --summary-json= "
                         "--out-dir= --rel-err-max=\n";
            std::exit(0);
        }
    }
}

void relax_solid(RealBody &body, BaseInnerRelation &inner)
{
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_particles(body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(inner);
    random_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    body.updateCellLinkedList();
    for (int ite_p = 0; ite_p < 500; ++ite_p)
    {
        relaxation_step_inner.exec();
        if ((ite_p + 1) % 100 == 0)
            std::cout << "[ophelie][doc228][S0] relax N=" << (ite_p + 1) << "\n";
    }
}
} // namespace

int main(int ac, char *av[])
{
    S0Params p;
    applyCli(ac, av, p);
    fs::create_directories(p.out_dir);

    const OphelieClientStainless304Development stainless = makeClientStainless304Development();
    writeClientStainless304Json(stainless, (fs::path(p.out_dir) / "stainless_304_dev_preset.json").string());

    const Real I = p.width * p.height * p.height * p.height / Real(12);
    if (!(p.tip_force_n > TinyReal))
    {
        const Real delta_target = Real(0.02) * p.length;
        p.tip_force_n = delta_target * (Real(3) * p.youngs_modulus * I) / std::pow(p.length, Real(3));
    }
    const Real delta_analytic =
        p.tip_force_n * std::pow(p.length, Real(3)) / (Real(3) * p.youngs_modulus * I);

    const Real physical_viscosity = get_physical_viscosity_general(p.rho, p.youngs_modulus, p.height);
    const Real dp = p.height / Real(5);
    const Real extension_length = Real(4) * dp;

    const Vec3d halfsize = Real(0.5) * (Vec3d(p.length, p.height, p.width) + extension_length * Vec3d::UnitX());
    const Vec3d translation = Real(0.5) * (p.length - extension_length) * Vec3d::UnitX();
    auto mesh = makeShared<GeometricShapeBox>(Transform(translation), halfsize, "beam");
    BoundingBoxd bb_system = mesh->getBounds();

    SPHSystem system(bb_system, dp);
    SolidBody beam_body(system, mesh);
    beam_body.defineBodyLevelSetShape();
    beam_body.defineMatterMaterial<NeoHookeanSolid>(p.rho, p.youngs_modulus, p.poisson);
    beam_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner(beam_body);
    relax_solid(beam_body, inner);
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> stress_first(inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_second(inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> damping(
        Real(0.2), inner, "Velocity", physical_viscosity);

    FixPart fix_bc_part(beam_body);
    SimpleDynamics<FixBodyPartConstraint> fix_bc(fix_bc_part);

    // Proven tip-load path (test_2d_plate): kernel-weighted point force → ForcePrior.
    StdVec<Vecd> point_forces = {Vecd(0, -p.tip_force_n, 0)};
    StdVec<Vecd> reference_positions = {Vecd(p.length, 0, 0)};
    SimpleDynamics<solid_dynamics::DistributingPointForces> tip_force(
        beam_body, point_forces, reference_positions, p.force_ramp_end, dp);

    BodyStatesRecordingToVtp vtp_output(system);
    vtp_output.addDerivedVariableRecording<SimpleDynamics<Displacement>>(beam_body);
    vtp_output.addToWrite<Vecd>(beam_body, "ForcePrior");

    const StdVec<Vecd> observation_locations = {Vecd(p.length, 0, 0)};
    ObserverBody observer(system, "TipObserver");
    observer.generateParticles<ObserverParticles>(observation_locations);
    ContactRelation observer_contact(observer, {&beam_body});
    ObservedQuantityRecording<Vecd> write_position("Position", observer_contact);

    corrected_configuration.exec();
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(beam_body);
    const Real dt_ref = computing_time_step_size.exec();

    std::cout << "[ophelie][doc228][S0] tip-force(DistributingPointForces) L=" << p.length
              << " h=w=" << p.height << " L/h=" << (p.length / p.height) << " dp=" << dp
              << " E=" << p.youngs_modulus << " F=" << p.tip_force_n
              << " delta_analytic=" << delta_analytic << " dt_ref=" << dt_ref
              << " n=" << beam_body.getBaseParticles().TotalRealParticles() << "\n";
    std::cout << "[ophelie][doc228][S0] stainless_preset E=" << stainless.youngs_modulus_pa
              << " Pa reserved for paddle E0 (this S0 uses soft E)\n";

    vtp_output.writeToFile(0);
    write_position.writeToFile(0);
    const Vecd tip_pos0 = write_position.getObservedQuantity()[0];

    int ite = 0;
    int ite_output = 0;
    const Real output_period = p.end_time / Real(20);
    Real dt = Real(0);
    const TickCount t1 = TickCount::now();

    while (physical_time < p.end_time)
    {
        Real integral_time = Real(0);
        while (integral_time < output_period)
        {
            if (ite % 500 == 0)
            {
                std::cout << "[ophelie][doc228][S0] N=" << ite << " t=" << physical_time << " dt=" << dt << "\n";
            }
            tip_force.exec(dt);
            dt = computing_time_step_size.exec();
            if (dt < dt_ref / Real(100))
            {
                std::cerr << "[ophelie][doc228][S0] FATAL: dt collapsed\n";
                return 2;
            }
            stress_first.exec(dt);
            fix_bc.exec();
            damping.exec(dt);
            fix_bc.exec();
            stress_second.exec(dt);
            ++ite;
            integral_time += dt;
            physical_time += dt;
        }
        ++ite_output;
        // Keep observer neighbors valid after large tip motion.
        observer_contact.updateConfiguration();
        vtp_output.writeToFile(ite_output);
        write_position.writeToFile(ite_output);
    }

    const Real wall_s = (TickCount::now() - t1).seconds();
    const Vecd tip_pos = write_position.getObservedQuantity()[0];
    const Real delta_num = tip_pos0[1] - tip_pos[1]; // tip drops in -Y
    const Real rel_err =
        (std::abs(delta_analytic) > TinyReal) ? std::abs(delta_num - delta_analytic) / std::abs(delta_analytic)
                                              : Real(1);
    const bool finite_ok = std::isfinite(static_cast<double>(delta_num));
    const bool pass = finite_ok && (rel_err <= p.rel_err_max) && (std::abs(delta_num) > TinyReal);

    std::cout << "[ophelie][doc228][S0] tip_pos0=(" << tip_pos0[0] << "," << tip_pos0[1] << "," << tip_pos0[2]
              << ") tip_pos=(" << tip_pos[0] << "," << tip_pos[1] << "," << tip_pos[2]
              << ") delta_num=" << delta_num << " delta_analytic=" << delta_analytic
              << " rel_err=" << rel_err << " wall_s=" << wall_s
              << " SOLID_REFERENCE_PASS=" << (pass ? 1 : 0) << "\n";

    const fs::path summary_path =
        fs::path(p.summary_json).is_absolute() ? fs::path(p.summary_json) : fs::path(p.out_dir) / p.summary_json;
    {
        std::ofstream js(summary_path);
        js << std::setprecision(10);
        js << "{\n";
        js << "  \"case\": \"doc228_s0_cantilever_tip_force\",\n";
        js << "  \"tag\": \"SOLID_REFERENCE\",\n";
        js << "  \"passed\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"SOLID_REFERENCE_PASS\": " << (pass ? 1 : 0) << ",\n";
        js << "  \"formula\": \"delta=F*L^3/(3*E*I)\",\n";
        js << "  \"operator_path\": \"NeoHookeanSolid+Integration1stHalfPK2RightCauchy+DistributingPointForces\",\n";
        js << "  \"L_m\": " << p.length << ",\n";
        js << "  \"height_m\": " << p.height << ",\n";
        js << "  \"width_m\": " << p.width << ",\n";
        js << "  \"dp_m\": " << dp << ",\n";
        js << "  \"n_particles\": " << beam_body.getBaseParticles().TotalRealParticles() << ",\n";
        js << "  \"rho\": " << p.rho << ",\n";
        js << "  \"E_Pa\": " << p.youngs_modulus << ",\n";
        js << "  \"poisson\": " << p.poisson << ",\n";
        js << "  \"tip_force_N\": " << p.tip_force_n << ",\n";
        js << "  \"I_m4\": " << I << ",\n";
        js << "  \"delta_analytic_m\": " << delta_analytic << ",\n";
        js << "  \"delta_numeric_m\": " << delta_num << ",\n";
        js << "  \"rel_err\": " << rel_err << ",\n";
        js << "  \"rel_err_max\": " << p.rel_err_max << ",\n";
        js << "  \"end_time\": " << p.end_time << ",\n";
        js << "  \"wall_seconds\": " << wall_s << ",\n";
        js << "  \"stainless_preset_json\": \"stainless_304_dev_preset.json\",\n";
        js << "  \"note\": \"S0 soft-E tip-force analytic gate; stainless E=2e11 for paddle E0/E1\"\n";
        js << "}\n";
    }
    std::cout << "[ophelie][doc228][S0] summary_json=" << summary_path.string() << "\n";
    return pass ? 0 : 1;
}
