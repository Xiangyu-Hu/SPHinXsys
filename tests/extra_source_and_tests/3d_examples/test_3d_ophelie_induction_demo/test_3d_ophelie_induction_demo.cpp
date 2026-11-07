/**
 * @file test_3d_ophelie_induction_demo.cpp
 * @brief Induction-welding demo: Case A (solid cylinder) / Case B (pipe joint marker).
 *
 * Stages: geometry | biot | em (default). VTP outputs for ParaView.
 * Source-field approximation (coil Biot only). Demo — not process / weld prediction.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

#if SPHINXSYS_USE_SYCL
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
}
#else
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
}
#endif

/** Continuous hollow aluminum pipe (Case B): single conductive domain + joint marker. */
class InductionDemoPipeShape : public ComplexShape
{
  public:
    InductionDemoPipeShape(const std::string &shape_name, const Vecd &center, Real r_outer, Real r_inner,
                           Real half_height, int mesh_resolution = 24)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_outer, half_height, mesh_resolution, center, "PipeOuter");
        // Slightly taller bore so level-set subtraction cleanly opens both ends.
        subtract<TriangleMeshShapeCylinder>(Vec3d(0, 0, 1), r_inner, half_height * Real(1.02), mesh_resolution, center,
                                            "PipeInner");
    }
};

struct InductionDemoCli
{
    char case_id = 'a';
    Real dp = 0.004;
    // Case A solid cylinder
    Real work_radius = 0.020;
    Real work_half_height = 0.050;
    // Case B pipe
    Real pipe_r_outer = 0.025;
    Real pipe_r_inner = 0.015;
    Real pipe_half_height = 0.060; // L = 0.120 m
    Real joint_half_length = 0.010; // L_joint = 0.020 m
    // Coil
    Real coil_radius = 0.035;
    size_t num_loops = 5;
    Real coil_pitch = 0.012;
    Real coil_z_center = 0.0;
    Real i_peak = 500.0;
    Real frequency_hz = 10000.0;
    Real sigma = 3.5e7;
    size_t segments_per_loop = 120;
    Real softening_length = 1.0e-6;
    std::string stage = "em"; // geometry | biot | em
};

inline void applyCaseDefaults(InductionDemoCli &cli)
{
    if (cli.case_id == 'b' || cli.case_id == 'B')
    {
        cli.case_id = 'b';
        cli.coil_radius = 0.040;
        cli.num_loops = 3;
        cli.coil_pitch = 0.010;
        cli.work_radius = cli.pipe_r_outer; // for clearance / domain
        cli.work_half_height = cli.pipe_half_height;
    }
    else
    {
        cli.case_id = 'a';
        cli.coil_radius = 0.035;
        cli.num_loops = 5;
        cli.coil_pitch = 0.012;
    }
}

inline void applyLocalCli(int ac, char *av[], InductionDemoCli &cli)
{
    char requested = 'a';
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--case=", 7) == 0)
        {
            const char c = av[i][7];
            if (c == 'a' || c == 'A' || c == 'b' || c == 'B')
            {
                requested = static_cast<char>(std::tolower(c));
            }
            else
            {
                std::cerr << "induction_demo: --case=a|b only\n";
                std::exit(2);
            }
        }
    }
    cli.case_id = requested;
    applyCaseDefaults(cli);

    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--i-peak=", 9) == 0)
        {
            cli.i_peak = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--frequency=", 12) == 0)
        {
            cli.frequency_hz = static_cast<Real>(std::atof(av[i] + 12));
        }
        else if (std::strncmp(av[i], "--sigma=", 8) == 0)
        {
            cli.sigma = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--stage=", 8) == 0)
        {
            cli.stage = av[i] + 8;
        }
        else if (std::strncmp(av[i], "--coil-z-center=", 16) == 0)
        {
            cli.coil_z_center = static_cast<Real>(std::atof(av[i] + 16));
        }
    }
}

inline OphelieMultiloopCoilSpec makeCoil(const InductionDemoCli &cli)
{
    OphelieMultiloopCoilSpec coil;
    coil.stack_center = Vecd(0.0, 0.0, cli.coil_z_center);
    coil.loop_radius = cli.coil_radius;
    const Real half_span = (cli.num_loops <= 1) ? Real(0)
                                                : 0.5 * static_cast<Real>(cli.num_loops - 1) * cli.coil_pitch;
    coil.z_min = cli.coil_z_center - half_span;
    coil.z_max = cli.coil_z_center + half_span;
    coil.num_loops = cli.num_loops;
    coil.current_per_loop = cli.i_peak;
    coil.segments_per_loop = cli.segments_per_loop;
    coil.use_cell_centered_loops = false;
    return coil;
}

inline bool geometryClearanceOk(const InductionDemoCli &cli, const OphelieMultiloopCoilSpec &coil, std::string &msg)
{
    const Real work_r = (cli.case_id == 'b') ? cli.pipe_r_outer : cli.work_radius;
    const Real radial_gap = coil.loop_radius - work_r;
    if (radial_gap <= Real(0))
    {
        msg = "coil intersects work radially (gap<=0)";
        return false;
    }
    std::ostringstream oss;
    oss << "case=" << cli.case_id << " radial_gap_m=" << radial_gap << " coil_z=[" << coil.z_min << "," << coil.z_max
        << "] work_z=[" << -cli.work_half_height << "," << cli.work_half_height << "]";
    if (cli.case_id == 'b')
    {
        oss << " pipe_r=[" << cli.pipe_r_inner << "," << cli.pipe_r_outer << "]";
    }
    msg = oss.str();
    return true;
}

inline void writeCoilPolylineVtp(const std::string &path, const OphelieMultiloopCoilSpec &coil)
{
    StdVec<OphelieCurrentMomentSample> moments;
    buildMultiloopFilamentMoments(coil, moments);
    const size_t n = moments.size();
    std::ofstream out(path);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <PolyData>\n";
    out << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << n
        << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd &p = moments[i].position;
        out << "          " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Lines>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t idx = 0;
    for (size_t k = 0; k < coil.num_loops; ++k)
    {
        for (size_t s = 0; s < coil.segments_per_loop; ++s)
        {
            const size_t a = idx + s;
            const size_t b = idx + ((s + 1) % coil.segments_per_loop);
            out << "          " << a << " " << b << "\n";
        }
        idx += coil.segments_per_loop;
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e = 1; e <= n; ++e)
    {
        out << "          " << (2 * e) << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Lines>\n";
    out << "    </Piece>\n";
    out << "  </PolyData>\n";
    out << "</VTKFile>\n";
}

/** Nominal joint ring at z=0 (marker only — not a weld interface). */
inline void writeJointMarkerVtp(const std::string &path, Real r_inner, Real r_outer, size_t n_seg = 72)
{
    const size_t n = 2 * n_seg;
    std::ofstream out(path);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <PolyData>\n";
    out << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << n
        << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < n_seg; ++i)
    {
        const Real th = 2.0 * Pi * static_cast<Real>(i) / static_cast<Real>(n_seg);
        out << "          " << (r_inner * std::cos(th)) << " " << (r_inner * std::sin(th)) << " 0\n";
    }
    for (size_t i = 0; i < n_seg; ++i)
    {
        const Real th = 2.0 * Pi * static_cast<Real>(i) / static_cast<Real>(n_seg);
        out << "          " << (r_outer * std::cos(th)) << " " << (r_outer * std::sin(th)) << " 0\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Lines>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t base : {size_t(0), n_seg})
    {
        for (size_t s = 0; s < n_seg; ++s)
        {
            out << "          " << (base + s) << " " << (base + ((s + 1) % n_seg)) << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e = 1; e <= n; ++e)
    {
        out << "          " << (2 * e) << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Lines>\n";
    out << "    </Piece>\n";
    out << "  </PolyData>\n";
    out << "</VTKFile>\n";
}

inline void fillComplexVecMag(BaseParticles &particles, const std::string &real_name, const std::string &imag_name,
                              const std::string &mag_name)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, real_name);
    syncVariableToHost<Vecd>(particles, imag_name);
    syncVariableToHost<Real>(particles, mag_name);
    const Vecd *vr = particles.getVariableDataByName<Vecd>(real_name);
    const Vecd *vi = particles.getVariableDataByName<Vecd>(imag_name);
    Real *mag = particles.getVariableDataByName<Real>(mag_name);
    for (size_t i = 0; i < n; ++i)
    {
        mag[i] = std::sqrt(vr[i].squaredNorm() + vi[i].squaredNorm());
    }
    syncVariableToDevice<Real>(particles, mag_name);
}

inline void fillJointMarkerFields(BaseParticles &particles, Real joint_half_length, const std::string &joint_name,
                                  const std::string &segment_name)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, joint_name);
    syncVariableToHost<Real>(particles, segment_name);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *joint = particles.getVariableDataByName<Real>(joint_name);
    Real *seg = particles.getVariableDataByName<Real>(segment_name);
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        joint[i] = (std::abs(z) <= joint_half_length) ? Real(1) : Real(0);
        seg[i] = (z >= Real(0)) ? Real(1) : Real(0); // nominal right/left pipe label
    }
    syncVariableToDevice<Real>(particles, joint_name);
    syncVariableToDevice<Real>(particles, segment_name);
}

inline Real hostScalarFieldMaxFinite(BaseParticles &particles, const std::string &field, size_t n, bool &all_finite)
{
    syncVariableToHost<Real>(particles, field);
    const Real *v = particles.getVariableDataByName<Real>(field);
    Real m = 0.0;
    all_finite = true;
    for (size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(v[i]))
        {
            all_finite = false;
            continue;
        }
        m = std::max(m, v[i]);
    }
    return m;
}

inline Real hostScalarFieldMin(BaseParticles &particles, const std::string &field, size_t n)
{
    syncVariableToHost<Real>(particles, field);
    const Real *v = particles.getVariableDataByName<Real>(field);
    Real m = v[0];
    for (size_t i = 1; i < n; ++i)
    {
        m = std::min(m, v[i]);
    }
    return m;
}

inline Real hostVolWeightedNormComplexVec(BaseParticles &particles, const std::string &real_name,
                                          const std::string &imag_name, size_t n)
{
    syncVariableToHost<Vecd>(particles, real_name);
    syncVariableToHost<Vecd>(particles, imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *vr = particles.getVariableDataByName<Vecd>(real_name);
    const Vecd *vi = particles.getVariableDataByName<Vecd>(imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real acc = 0.0;
    Real vol_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        acc += (vr[i].squaredNorm() + vi[i].squaredNorm()) * vol[i];
        vol_sum += vol[i];
    }
    return std::sqrt(acc / (vol_sum + TinyReal));
}

struct JointPowerMetrics
{
    Real p_total = 0.0;
    Real p_joint = 0.0;
    Real eta_joint = 0.0;
};

inline JointPowerMetrics computeJointPower(BaseParticles &particles, const std::string &q_field, Real joint_half_length)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, q_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    JointPowerMetrics m;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dP = q[i] * vol[i];
        m.p_total += dP;
        if (std::abs(pos[i][2]) <= joint_half_length)
        {
            m.p_joint += dP;
        }
    }
    m.eta_joint = m.p_joint / (m.p_total + TinyReal);
    return m;
}

inline void writePowerDensityVsZ(const std::string &path, BaseParticles &particles, const std::string &q_field,
                                 Real z_min, Real z_max, size_t n_bins)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, q_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real dz = (z_max - z_min) / static_cast<Real>(n_bins);
    std::vector<Real> p_bin(n_bins, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        const Real z = pos[i][2];
        if (z < z_min || z > z_max)
        {
            continue;
        }
        size_t b = static_cast<size_t>((z - z_min) / (dz + TinyReal));
        if (b >= n_bins)
        {
            b = n_bins - 1;
        }
        p_bin[b] += q[i] * vol[i];
    }
    std::ofstream out(path);
    out << "z_center_m,P_prime_W_per_m,P_bin_W\n";
    out.setf(std::ios::scientific);
    out.precision(10);
    for (size_t b = 0; b < n_bins; ++b)
    {
        const Real zc = z_min + (static_cast<Real>(b) + Real(0.5)) * dz;
        const Real p_prime = p_bin[b] / (dz + TinyReal);
        out << zc << "," << p_prime << "," << p_bin[b] << "\n";
    }
}

inline void writeAxisBzCsv(const std::string &path, const OphelieMultiloopCoilSpec &coil, Real mu0, Real softening)
{
    StdVec<OphelieCurrentMomentSample> moments;
    buildMultiloopFilamentMoments(coil, moments);
    std::ofstream out(path);
    out << "z_m,Bz_filament_T,Bz_analytic_T,rel_err\n";
    out.setf(std::ios::scientific);
    out.precision(10);
    for (int k = -24; k <= 24; ++k)
    {
        const Real z = 0.005 * static_cast<Real>(k);
        const Real bz_f = filamentBzOnAxisFromMoments(moments, Vecd(0, 0, z), mu0, softening);
        const Real bz_a = analyticMultiloopBzOnAxis(z, coil, mu0);
        const Real rel = std::abs(bz_f - bz_a) / (std::abs(bz_a) + TinyReal);
        out << z << "," << bz_f << "," << bz_a << "," << rel << "\n";
    }
}

inline void writeParamsTable(const std::string &path, const InductionDemoCli &cli, const OphelieMultiloopCoilSpec &coil,
                             size_t n_particles, Real max_b, Real linear_ratio_b, Real p_recon, Real max_j, Real max_q,
                             Real phi_eq_res, Real eta_joint, bool passed)
{
    std::ofstream out(path);
    out << "# OPHELIE induction demo\n";
    out << "# demo / not validated for process design; joint marker is NOT a weld\n";
    out << "# source-field approximation: coil Biot only\n";
    out << "case," << cli.case_id << "\n";
    out << "stage," << cli.stage << "\n";
    out << "material,aluminum_nonmagnetic_constant\n";
    if (cli.case_id == 'b')
    {
        out << "pipe_r_outer_m," << cli.pipe_r_outer << "\n";
        out << "pipe_r_inner_m," << cli.pipe_r_inner << "\n";
        out << "pipe_length_m," << (2.0 * cli.pipe_half_height) << "\n";
        out << "joint_half_length_m," << cli.joint_half_length << "\n";
        out << "eta_joint," << eta_joint << "\n";
    }
    else
    {
        out << "work_radius_m," << cli.work_radius << "\n";
        out << "work_length_m," << (2.0 * cli.work_half_height) << "\n";
    }
    out << "coil_radius_m," << coil.loop_radius << "\n";
    out << "num_loops," << coil.num_loops << "\n";
    out << "coil_pitch_m," << cli.coil_pitch << "\n";
    out << "coil_z_center_m," << cli.coil_z_center << "\n";
    out << "I_peak_A," << coil.current_per_loop << "\n";
    out << "frequency_Hz," << cli.frequency_hz << "\n";
    out << "sigma_S_per_m," << cli.sigma << "\n";
    out << "dp_m," << cli.dp << "\n";
    out << "n_particles," << n_particles << "\n";
    out << "max_B_src_mag_T," << max_b << "\n";
    out << "I_to_2I_Bmag_ratio," << linear_ratio_b << "\n";
    out << "P_recon_W," << p_recon << "\n";
    out << "max_J_mag_A_per_m2," << max_j << "\n";
    out << "max_JouleHeat_W_per_m3," << max_q << "\n";
    out << "phi_eq_res_vol_imag," << phi_eq_res << "\n";
    out << "passed," << (passed ? 1 : 0) << "\n";
}

} // namespace

int main(int ac, char *av[])
{
    InductionDemoCli local;
    applyLocalCli(ac, av, local);

    std::cout << "[ophelie][induction_demo] demo / not validated for process design\n"
              << "[ophelie][induction_demo] case=" << local.case_id << " stage=" << local.stage
              << " coil_z_center=" << local.coil_z_center
              << " note=source_field_coil_Biot; skin_depth_unresolved_at_dp";
    if (local.case_id == 'b')
    {
        std::cout << "; joint_marker_only_not_a_weld";
    }
    std::cout << std::endl;

    const OphelieMultiloopCoilSpec coil = makeCoil(local);
    std::string clearance_msg;
    if (!geometryClearanceOk(local, coil, clearance_msg))
    {
        std::cerr << "geometry clearance failed: " << clearance_msg << std::endl;
        return 1;
    }
    std::cout << "[ophelie][induction_demo] geometry " << clearance_msg << std::endl;

    const Vecd center = Vecd::Zero();
    const Real pad = Real(3) * local.dp;
    const Real work_r = (local.case_id == 'b') ? local.pipe_r_outer : local.work_radius;
    const Real r_domain = std::max(local.coil_radius, work_r) + pad;
    const BoundingBoxd bounds(Vecd(-r_domain, -r_domain, -local.work_half_height - pad),
                              Vecd(r_domain, r_domain, local.work_half_height + pad));

    SPHSystem sph_system(bounds, local.dp);
    SharedPtr<Shape> work_shape =
        (local.case_id == 'b')
            ? SharedPtr<Shape>(makeShared<InductionDemoPipeShape>("WorkBody", center, local.pipe_r_outer,
                                                                  local.pipe_r_inner, local.pipe_half_height))
            : SharedPtr<Shape>(makeShared<OphelieFrenchReducedGlassCylinderShape>(
                  "WorkBody", center, local.work_radius, local.work_half_height));
    SolidBody work_body(sph_system, work_shape);
    work_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    work_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(work_body);
    work_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames names;
    RegisterOphelieGlassFields register_fields(work_body, names);
    (void)register_fields;

    BaseParticles &particles = work_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const std::string b_mag_name = "BSrcMag";
    const std::string j_mag_name = "JEdgeReconMag";
    const std::string e_mag_name = "EEdgeReconMag";
    const std::string joint_name = "JointMarker";
    const std::string segment_name = "PipeSegment";
    particles.registerStateVariable<Real>(b_mag_name, Real(0));
    particles.registerStateVariable<Real>(j_mag_name, Real(0));
    particles.registerStateVariable<Real>(e_mag_name, Real(0));
    particles.registerStateVariable<Real>(joint_name, Real(0));
    particles.registerStateVariable<Real>(segment_name, Real(0));
    particles.addVariableToWrite<Real>(b_mag_name);
    particles.addVariableToWrite<Real>(j_mag_name);
    particles.addVariableToWrite<Real>(e_mag_name);
    particles.addVariableToWrite<Real>(joint_name);
    particles.addVariableToWrite<Real>(segment_name);

    if (local.case_id == 'b')
    {
        fillJointMarkerFields(particles, local.joint_half_length, joint_name, segment_name);
    }

    UniquePtr<Inner<>> work_inner = makeUnique<Inner<>>(work_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(work_body, names, local.sigma);
    assign_sigma.exec();

    const std::string out_dir = IO::getEnvironment().OutputFolder();
    const std::string case_tag = std::string("case_") + local.case_id;
    std::ostringstream coil_suffix;
    coil_suffix << case_tag;
    if (std::abs(local.coil_z_center) > TinyReal)
    {
        coil_suffix << "_zcoil_" << local.coil_z_center;
    }
    const std::string tag = coil_suffix.str();
    const std::string coil_vtp = out_dir + "/CoilFilament_" + tag + ".vtp";
    writeCoilPolylineVtp(coil_vtp, coil);
    if (local.case_id == 'b')
    {
        writeJointMarkerVtp(out_dir + "/JointMarker_" + tag + ".vtp", local.pipe_r_inner, local.pipe_r_outer);
    }

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(work_body, names.sigma);
    write_states.addToWrite<Vecd>(work_body, names.b_src_real);
    write_states.addToWrite<Vecd>(work_body, names.b_src_imag);
    write_states.addToWrite<Vecd>(work_body, names.a_src_real);
    write_states.addToWrite<Vecd>(work_body, names.a_src_imag);
    write_states.addToWrite<Real>(work_body, b_mag_name);
    write_states.addToWrite<Real>(work_body, names.phi_real);
    write_states.addToWrite<Real>(work_body, names.phi_imag);
    write_states.addToWrite<Vecd>(work_body, names.j_edge_recon_real);
    write_states.addToWrite<Vecd>(work_body, names.j_edge_recon_imag);
    write_states.addToWrite<Vecd>(work_body, names.e_edge_recon_real);
    write_states.addToWrite<Vecd>(work_body, names.e_edge_recon_imag);
    write_states.addToWrite<Real>(work_body, names.joule_heat_edge_recon_complex);
    write_states.addToWrite<Real>(work_body, j_mag_name);
    write_states.addToWrite<Real>(work_body, e_mag_name);
    write_states.addToWrite<Real>(work_body, joint_name);
    write_states.addToWrite<Real>(work_body, segment_name);

    Real max_b = 0.0;
    Real linear_ratio_b = 0.0;
    Real p_recon = 0.0;
    Real max_j = 0.0;
    Real max_q = 0.0;
    Real phi_eq_res = 0.0;
    Real eta_joint = 0.0;
    bool finite_ok = (n > 0);
    bool linear_ok = false;
    bool axis_ok = true;
    bool em_ok = (local.stage != "em");

    auto write_params = [&](bool passed) {
        writeParamsTable(out_dir + "/" + tag + "_params.csv", local, coil, n, max_b, linear_ratio_b, p_recon, max_j,
                         max_q, phi_eq_res, eta_joint, passed);
    };

    if (local.stage == "geometry")
    {
        work_body.setNewlyUpdated();
        write_states.writeToFile(0);
        write_params(true);
        std::cout << "[ophelie][induction_demo] VTP work=" << out_dir << "/WorkBody_ite_0000000000.vtp"
                  << " coil=" << coil_vtp << std::endl;
        std::cout << "test_3d_ophelie_induction_demo n=" << n << " case=" << local.case_id
                  << " stage=geometry passed=1" << std::endl;
        return 0;
    }

    OphelieParameters params;
    params.mu0_ = 4.0 * Pi * 1.0e-7;
    params.softening_length_ = local.softening_length;
    params.sigma_glass_ = local.sigma;
    params.frequency_ = local.frequency_hz;
    params.enable_power_scaling_ = false;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = true;
    params.aind_one_way_feedback_resolve_ = false;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::SolverLocal;

    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);
    logOphelieFinalParams(params, cli_options);

    applyMultiloopFilamentBiotToGlass(particles, names, coil, params.mu0_, params.softening_length_);
    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);
    syncVariableToDevice<Vecd>(particles, names.a_ind_real);
    syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
    syncVariableToDevice<Vecd>(particles, names.b_ind_real);
    syncVariableToDevice<Vecd>(particles, names.b_ind_imag);

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine(work_body, names);
    combine.exec();
    fillComplexVecMag(particles, names.b_src_real, names.b_src_imag, b_mag_name);

    bool mag_finite = true;
    max_b = hostScalarFieldMaxFinite(particles, b_mag_name, n, mag_finite);
    const Real b_norm_I = hostVolWeightedNormComplexVec(particles, names.b_src_real, names.b_src_imag, n);
    finite_ok = finite_ok && mag_finite && std::isfinite(max_b) && max_b > TinyReal;

    writeAxisBzCsv(out_dir + "/" + tag + "_axis_Bz.csv", coil, params.mu0_, params.softening_length_);
    {
        StdVec<OphelieCurrentMomentSample> moments;
        buildMultiloopFilamentMoments(coil, moments);
        Real max_axis_rel = 0.0;
        const Real z_samples[] = {coil.z_min, Real(0.5) * (coil.z_min + coil.z_max), coil.z_max, local.coil_z_center,
                                  local.work_half_height * Real(0.8)};
        for (const Real z : z_samples)
        {
            const Real bz_f = filamentBzOnAxisFromMoments(moments, Vecd(0, 0, z), params.mu0_, params.softening_length_);
            const Real bz_a = analyticMultiloopBzOnAxis(z, coil, params.mu0_);
            const Real rel = std::abs(bz_f - bz_a) / (std::abs(bz_a) + TinyReal);
            max_axis_rel = std::max(max_axis_rel, rel);
            std::cout << "[ophelie][induction_demo] axis z=" << z << " Bz_fil=" << bz_f << " Bz_an=" << bz_a
                      << " rel_err=" << rel << std::endl;
        }
        axis_ok = max_axis_rel < Real(0.02);
        std::cout << "[ophelie][induction_demo] axis max_rel_err=" << max_axis_rel << " ok=" << (axis_ok ? 1 : 0)
                  << std::endl;
    }

    OphelieMultiloopCoilSpec coil_2i = coil;
    coil_2i.current_per_loop = Real(2) * coil.current_per_loop;
    applyMultiloopFilamentBiotToGlass(particles, names, coil_2i, params.mu0_, params.softening_length_);
    combine.exec();
    const Real b_norm_2I = hostVolWeightedNormComplexVec(particles, names.b_src_real, names.b_src_imag, n);
    linear_ratio_b = b_norm_2I / (b_norm_I + TinyReal);
    linear_ok = std::isfinite(linear_ratio_b) && std::abs(linear_ratio_b - Real(2)) < Real(0.05);

    applyMultiloopFilamentBiotToGlass(particles, names, coil, params.mu0_, params.softening_length_);
    combine.exec();
    fillComplexVecMag(particles, names.b_src_real, names.b_src_imag, b_mag_name);

    if (local.stage == "biot")
    {
        work_body.setNewlyUpdated();
        write_states.writeToFile(0);
        const bool passed = finite_ok && linear_ok && axis_ok;
        write_params(passed);
        std::cout << "[ophelie][induction_demo] VTP work=" << out_dir << "/WorkBody_ite_0000000000.vtp"
                  << " coil=" << coil_vtp << std::endl;
        std::cout << "test_3d_ophelie_induction_demo n=" << n << " case=" << local.case_id
                  << " stage=biot max_B_src_mag_T=" << max_b << " I_to_2I_ratio=" << linear_ratio_b
                  << " passed=" << (passed ? 1 : 0) << std::endl;
        return passed ? 0 : 1;
    }

    OphelieComplexEdgeFluxSolveReport solve_report;
    p_recon = execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
        work_body, *work_inner, names, params, nullptr, local.dp, &solve_report);
    phi_eq_res = solve_report.phi_eq_res_vol_imag;

    fillComplexVecMag(particles, names.b_src_real, names.b_src_imag, b_mag_name);
    fillComplexVecMag(particles, names.j_edge_recon_real, names.j_edge_recon_imag, j_mag_name);
    fillComplexVecMag(particles, names.e_edge_recon_real, names.e_edge_recon_imag, e_mag_name);
    if (local.case_id == 'b')
    {
        fillJointMarkerFields(particles, local.joint_half_length, joint_name, segment_name);
    }

    bool j_finite = true;
    bool q_finite = true;
    max_j = hostScalarFieldMaxFinite(particles, j_mag_name, n, j_finite);
    max_q = hostScalarFieldMaxFinite(particles, names.joule_heat_edge_recon_complex, n, q_finite);
    const Real min_q = hostScalarFieldMin(particles, names.joule_heat_edge_recon_complex, n);
    const bool q_nonneg = min_q >= -Real(1.0e-6) * (std::abs(max_q) + Real(1));
    em_ok = j_finite && q_finite && std::isfinite(p_recon) && p_recon > TinyReal && max_j > TinyReal && q_nonneg &&
            std::isfinite(phi_eq_res);

    if (local.case_id == 'b')
    {
        const JointPowerMetrics jm =
            computeJointPower(particles, names.joule_heat_edge_recon_complex, local.joint_half_length);
        eta_joint = jm.eta_joint;
        writePowerDensityVsZ(out_dir + "/" + tag + "_Pprime_z.csv", particles, names.joule_heat_edge_recon_complex,
                             -local.pipe_half_height, local.pipe_half_height, 24);
        std::cout << "[ophelie][induction_demo] joint P_total_W=" << jm.p_total << " P_joint_W=" << jm.p_joint
                  << " eta_joint=" << jm.eta_joint << " L_joint=" << (2.0 * local.joint_half_length) << std::endl;
    }

    work_body.setNewlyUpdated();
    write_states.writeToFile(0);

    const bool passed = finite_ok && linear_ok && axis_ok && em_ok;
    write_params(passed);

    std::cout << "[ophelie][induction_demo] em P_recon_W=" << p_recon << " max_J_mag=" << max_j << " max_q=" << max_q
              << " min_q=" << min_q << " phi_eq_res_vol_imag=" << phi_eq_res
              << " phi_imag_rel_res=" << solve_report.phi_imag_solver_rel_residual
              << " em_ok=" << (em_ok ? 1 : 0) << std::endl;
    std::cout << "[ophelie][induction_demo] VTP work=" << out_dir << "/WorkBody_ite_0000000000.vtp"
              << " coil=" << coil_vtp << " params=" << out_dir << "/" << tag << "_params.csv" << std::endl;
    std::cout << "test_3d_ophelie_induction_demo n=" << n << " case=" << local.case_id << " dp=" << local.dp
              << " R_coil=" << coil.loop_radius << " N=" << coil.num_loops << " I_peak=" << coil.current_per_loop
              << " f_hz=" << local.frequency_hz << " coil_z_center=" << local.coil_z_center << " stage=em"
              << " max_B_src_mag_T=" << max_b << " I_to_2I_B_ratio=" << linear_ratio_b << " P_recon_W=" << p_recon
              << " eta_joint=" << eta_joint << " max_J_mag=" << max_j << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
