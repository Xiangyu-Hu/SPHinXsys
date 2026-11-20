/**
 * @file electromagnetic_ophelie_doc216_e1_operator_export.h
 * @brief Doc216 E1: export real cylinder SPH G coefficients (COO), C / B_pre, particles.
 *
 * Reconstructs directed edge weights a_ij identical to ComputeOphelieE2ScalarGradientCK:
 *   (G phi)_i += a_ij * (phi_i - phi_j),  a_ij = -dW_ij Vol_j * e_hat
 * with e_hat = C_i e_ij when use_linear_correction, else e_ij.
 * Also rebuilds pre-invert B_i = Σ (-r_ij)(dW Vol_j e_ij)^T for SVD diagnostics.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC216_E1_OPERATOR_EXPORT_H
#define ELECTROMAGNETIC_OPHELIE_DOC216_E1_OPERATOR_EXPORT_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "relation_ck.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * @param inner              glass Inner<> relation (neighbors after last solve update)
 * @param use_linear_correction  match E2 path (true for --doc215-e1 / --doc214-g2)
 */
inline void hostOphelieDoc216E1OperatorExport(BaseParticles &particles, Inner<> &inner,
                                              const OphelieGlassFieldNames &names, bool use_linear_correction,
                                              Real radius, Real half_height, Real dp, Real h_ref, Real omega,
                                              Real sigma_const, Real b0, const Vecd &center,
                                              const std::string &prefix)
{
    const size_t n = particles.TotalRealParticles();
#if SPHINXSYS_USE_SYCL
    {
        using execution::ParallelDevicePolicy;
        ParallelDevicePolicy par{};
        inner.dvNeighborIndex(0)->prepareForOutput(par);
        inner.dvParticleOffset(0)->prepareForOutput(par);
    }
#endif
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Matd>(particles, names.phi_grad_linear_correction);
    syncVariableToHost<Real>(particles, names.phi_real);
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_src_real);
    syncVariableToHost<Vecd>(particles, names.a_src_imag);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *sig = particles.getVariableDataByName<Real>(names.sigma);
    Matd *C = particles.getVariableDataByName<Matd>(names.phi_grad_linear_correction);
    Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
    Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
    Vecd *Ar = particles.getVariableDataByName<Vecd>(names.a_src_real);
    Vecd *Ai = particles.getVariableDataByName<Vecd>(names.a_src_imag);
    Vecd *Er = particles.getVariableDataByName<Vecd>(names.e_edge_recon_real);
    Vecd *Ei = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    Vecd *Jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    Vecd *Ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    Real *Q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);

    UnsignedInt *neighbor_index = inner.dvNeighborIndex(0)->Data();
    UnsignedInt *particle_offset = inner.dvParticleOffset(0)->Data();

    execution::SequencedPolicy seq{};
    Neighbor<SPHAdaptation, SPHAdaptation>::SmoothingKernel sk(seq, inner.getNeighborhood(0));

    const std::string part_path = prefix + ".particles.csv";
    const std::string edge_path = prefix + ".edges_coo.csv";
    const std::string meta_path = prefix + ".meta.json";

    std::ofstream pofs(part_path);
    pofs << "i,x,y,z,Vol,sigma,phi_r,phi_i,"
            "Ar_x,Ar_y,Ar_z,Ai_x,Ai_y,Ai_z,"
            "Er_x,Er_y,Er_z,Ei_x,Ei_y,Ei_z,"
            "Jr_x,Jr_y,Jr_z,Ji_x,Ji_y,Ji_z,Q,"
            "C00,C01,C02,C10,C11,C12,C20,C21,C22,C_det,"
            "B00,B01,B02,B10,B11,B12,B20,B21,B22,B_det,n_nbr\n";
    pofs << std::scientific << std::setprecision(12);

    std::ofstream eofs(edge_path);
    eofs << "i,j,a_x,a_y,a_z,dW_ij,Vol_j,r_ij,corrected\n";
    eofs << std::scientific << std::setprecision(12);

    size_t n_edges = 0;
    size_t n_nbr_min = n > 0 ? size_t(-1) : 0;
    size_t n_nbr_max = 0;
    double n_nbr_sum = 0;

    for (size_t i = 0; i < n; ++i)
    {
        Matd B_pre = Matd::Zero();
        const UnsignedInt n0 = particle_offset[i];
        const UnsignedInt n1 = particle_offset[i + 1];
        const size_t nn = static_cast<size_t>(n1 - n0);
        n_nbr_min = std::min(n_nbr_min, nn);
        n_nbr_max = std::max(n_nbr_max, nn);
        n_nbr_sum += double(nn);

        for (UnsignedInt nslot = n0; nslot != n1; ++nslot)
        {
            const UnsignedInt j = neighbor_index[nslot];
            const Vecd r_ij = sk.vec_r_ij(static_cast<UnsignedInt>(i), j);
            const Real dW = sk.dW_ij(static_cast<UnsignedInt>(i), j);
            const Vecd e = sk.e_ij(static_cast<UnsignedInt>(i), j);
            const Real Vol_j = vol[j];
            // Pre-invert B accumulation (same as ComputeOpheliePhiGradLinearCorrectionMatrixCK).
            const Vecd gradW_ijV_j = dW * Vol_j * e;
            B_pre -= r_ij * gradW_ijV_j.transpose();

            const Vecd e_hat = use_linear_correction ? (C[i] * e) : e;
            const Vecd a_ij = pairwiseGradientWeightUncorrected(dW * Vol_j, e_hat);
            eofs << i << "," << j << "," << double(a_ij[0]) << "," << double(a_ij[1]) << "," << double(a_ij[2])
                 << "," << double(dW) << "," << double(Vol_j) << "," << double(r_ij.norm()) << ","
                 << (use_linear_correction ? 1 : 0) << "\n";
            ++n_edges;
        }

        const Matd &Ci = C[i];
        const Real Cdet = Ci.determinant();
        const Real Bdet = B_pre.determinant();
        pofs << i << "," << double(pos[i][0]) << "," << double(pos[i][1]) << "," << double(pos[i][2]) << ","
             << double(vol[i]) << "," << double(sig[i]) << "," << double(phi_r[i]) << "," << double(phi_i[i]) << ","
             << double(Ar[i][0]) << "," << double(Ar[i][1]) << "," << double(Ar[i][2]) << "," << double(Ai[i][0])
             << "," << double(Ai[i][1]) << "," << double(Ai[i][2]) << "," << double(Er[i][0]) << ","
             << double(Er[i][1]) << "," << double(Er[i][2]) << "," << double(Ei[i][0]) << "," << double(Ei[i][1])
             << "," << double(Ei[i][2]) << "," << double(Jr[i][0]) << "," << double(Jr[i][1]) << ","
             << double(Jr[i][2]) << "," << double(Ji[i][0]) << "," << double(Ji[i][1]) << "," << double(Ji[i][2])
             << "," << double(Q[i]) << "," << double(Ci(0, 0)) << "," << double(Ci(0, 1)) << "," << double(Ci(0, 2))
             << "," << double(Ci(1, 0)) << "," << double(Ci(1, 1)) << "," << double(Ci(1, 2)) << ","
             << double(Ci(2, 0)) << "," << double(Ci(2, 1)) << "," << double(Ci(2, 2)) << "," << double(Cdet) << ","
             << double(B_pre(0, 0)) << "," << double(B_pre(0, 1)) << "," << double(B_pre(0, 2)) << ","
             << double(B_pre(1, 0)) << "," << double(B_pre(1, 1)) << "," << double(B_pre(1, 2)) << ","
             << double(B_pre(2, 0)) << "," << double(B_pre(2, 1)) << "," << double(B_pre(2, 2)) << ","
             << double(Bdet) << "," << nn << "\n";
    }
    pofs.close();
    eofs.close();

    std::ofstream mofs(meta_path);
    mofs << std::scientific << std::setprecision(12);
    mofs << "{\n"
         << "  \"case\": \"Doc216_E1_operator_export\",\n"
         << "  \"n_particles\": " << n << ",\n"
         << "  \"n_edges\": " << n_edges << ",\n"
         << "  \"n_nbr_min\": " << (n > 0 ? n_nbr_min : 0) << ",\n"
         << "  \"n_nbr_max\": " << n_nbr_max << ",\n"
         << "  \"n_nbr_mean\": " << (n > 0 ? n_nbr_sum / double(n) : 0.0) << ",\n"
         << "  \"use_linear_correction\": " << (use_linear_correction ? "true" : "false") << ",\n"
         << "  \"G_layout\": \"3N_x_N_particle_xyz_interleaved; edge a_ij multiplies (phi_i-phi_j) into rows 3i:3i+3\",\n"
         << "  \"radius\": " << double(radius) << ",\n"
         << "  \"half_height\": " << double(half_height) << ",\n"
         << "  \"dp\": " << double(dp) << ",\n"
         << "  \"h_ref\": " << double(h_ref) << ",\n"
         << "  \"omega\": " << double(omega) << ",\n"
         << "  \"sigma_const\": " << double(sigma_const) << ",\n"
         << "  \"b0\": " << double(b0) << ",\n"
         << "  \"center\": [" << double(center[0]) << ", " << double(center[1]) << ", " << double(center[2])
         << "],\n"
         << "  \"particles_csv\": \"" << part_path << "\",\n"
         << "  \"edges_coo_csv\": \"" << edge_path << "\"\n"
         << "}\n";
    mofs.close();

    std::cout << std::scientific << std::setprecision(6)
              << "[ophelie][doc216][e1] export prefix=" << prefix << " n=" << n << " edges=" << n_edges
              << " nbr_mean=" << (n > 0 ? n_nbr_sum / double(n) : 0.0)
              << " corrected=" << (use_linear_correction ? 1 : 0) << "\n"
              << "[ophelie][doc216][e1] particles=" << part_path << "\n"
              << "[ophelie][doc216][e1] edges=" << edge_path << "\n"
              << "[ophelie][doc216][e1] meta=" << meta_path << "\n";
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC216_E1_OPERATOR_EXPORT_H
