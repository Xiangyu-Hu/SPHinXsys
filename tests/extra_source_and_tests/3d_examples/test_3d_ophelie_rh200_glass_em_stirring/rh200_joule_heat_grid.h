/**
 * @file rh200_joule_heat_grid.h
 * @brief Eulerian JouleHeat (and EM scalar) grid: host CIC deposit, device trilinear sample, VTI output.
 */
#ifndef RH200_JOULE_HEAT_GRID_H
#define RH200_JOULE_HEAT_GRID_H

#include "base_local_dynamics.h"
#include "rh200_fake_joule_heat.h"

#if SPHINXSYS_USE_SYCL
#include "implementation_sycl.h"
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace rh200
{

inline constexpr const char *kJouleHeatGridSummaryCsv = "./output/rh200_joule_heat_grid_summary.csv";
inline constexpr const char *kJouleHeatGridVti = "./output/JouleHeatGrid.vti";

struct Rh200EmParticleDepositionHost
{
    StdVec<Vecd> position;
    StdVec<Real> volumetric_measure;
    StdVec<Real> joule_heat;
    StdVec<Real> phi_real;
    StdVec<Real> phi_imag;
    StdVec<Vecd> e_real;
    StdVec<Vecd> e_imag;
    StdVec<Vecd> j_real;
    StdVec<Vecd> j_imag;
};

struct Rh200JouleHeatGridSpec
{
    Vecd lower_bound_ = Vecd::Zero();
    Vecd upper_bound_ = Vecd::Zero();
    Real spacing_ = Real(0.008);
    size_t nx_ = 0;
    size_t ny_ = 0;
    size_t nz_ = 0;
};

struct Rh200JouleHeatGridBuildReport
{
    Real p_em_particle = Real(0);
    Real p_grid_sample_before_scale = Real(0);
    Real grid_scale_factor = Real(1);
    Real p_grid_sample_after_scale = Real(0);
    Real sample_vs_em_rel_l2 = Real(0);
    Real sample_vs_em_rel_max = Real(0);
    Real q_min = Real(0);
    Real q_max = Real(0);
    Real q_mean = Real(0);
    Real q_vol_weighted_mean = Real(0);
    Real empty_node_fraction = Real(0);
    size_t nonempty_node_count = 0;
};

inline Rh200JouleHeatGridSpec makeRh200JouleHeatGridSpecFromBounds(const BoundingBoxd &bounds, Real grid_spacing,
                                                                   Real bbox_margin_in_spacing = Real(2))
{
    Rh200JouleHeatGridSpec spec;
    const Real margin = std::max(bbox_margin_in_spacing, Real(0)) * std::max(grid_spacing, Real(1.0e-6));
    spec.lower_bound_ =
        Vecd(bounds.lower_[0], bounds.lower_[1], bounds.lower_[2]) - Vecd::Constant(margin);
    spec.upper_bound_ =
        Vecd(bounds.upper_[0], bounds.upper_[1], bounds.upper_[2]) + Vecd::Constant(margin);
    spec.spacing_ = std::max(grid_spacing, Real(1.0e-6));
    const Vecd size = spec.upper_bound_ - spec.lower_bound_;
    spec.nx_ = std::max<size_t>(2, static_cast<size_t>(std::ceil(size[0] / spec.spacing_)) + 1);
    spec.ny_ = std::max<size_t>(2, static_cast<size_t>(std::ceil(size[1] / spec.spacing_)) + 1);
    spec.nz_ = std::max<size_t>(2, static_cast<size_t>(std::ceil(size[2] / spec.spacing_)) + 1);
    return spec;
}

inline Rh200JouleHeatGridSpec makeRh200JouleHeatGridSpecFromStl(const std::string &stl_path, Real geometry_scale,
                                                                Real grid_spacing, Real bbox_margin_in_spacing = Real(2))
{
    TriangleMeshShapeSTL shape(stl_path, Vec3d::Zero(), geometry_scale, "joule_grid_domain");
    return makeRh200JouleHeatGridSpecFromBounds(shape.findBounds(), grid_spacing, bbox_margin_in_spacing);
}

class Rh200ScalarEulerianGrid
{
  public:
    Rh200JouleHeatGridSpec spec_;
    StdVec<Real> q_host_;
    StdVec<Real> weight_host_;

    size_t nodeCount() const { return spec_.nx_ * spec_.ny_ * spec_.nz_; }

    size_t index(size_t i, size_t j, size_t k) const { return i + spec_.nx_ * (j + spec_.ny_ * k); }

    void resetAccumulators()
    {
        const size_t n = nodeCount();
        q_host_.assign(n, Real(0));
        weight_host_.assign(n, Real(0));
    }

    bool inside(const Vecd &x) const
    {
        if (spec_.nx_ < 2 || spec_.ny_ < 2 || spec_.nz_ < 2)
        {
            return false;
        }
        const Vecd rel = (x - spec_.lower_bound_) / spec_.spacing_;
        return rel[0] >= Real(0) && rel[0] <= static_cast<Real>(spec_.nx_ - 1) && rel[1] >= Real(0) &&
               rel[1] <= static_cast<Real>(spec_.ny_ - 1) && rel[2] >= Real(0) &&
               rel[2] <= static_cast<Real>(spec_.nz_ - 1);
    }

    void depositScalarCloudInCell(const StdVec<Vecd> &position, const StdVec<Real> &scalar,
                                  const StdVec<Real> &volumetric_measure)
    {
        const size_t n = position.size();
        for (size_t p = 0; p < n; ++p)
        {
            const Vecd x = position[p];
            const Real value = scalar[p];
            const Real vol = volumetric_measure[p];
            if (!(vol > TinyReal) || !std::isfinite(value))
            {
                continue;
            }
            const Vecd rel = (x - spec_.lower_bound_) / spec_.spacing_;
            if (rel[0] < Real(0) || rel[0] > static_cast<Real>(spec_.nx_ - 1) || rel[1] < Real(0) ||
                rel[1] > static_cast<Real>(spec_.ny_ - 1) || rel[2] < Real(0) ||
                rel[2] > static_cast<Real>(spec_.nz_ - 1))
            {
                continue;
            }
            const int i0 = static_cast<int>(std::floor(rel[0]));
            const int j0 = static_cast<int>(std::floor(rel[1]));
            const int k0 = static_cast<int>(std::floor(rel[2]));
            const int i1 = std::min(i0 + 1, static_cast<int>(spec_.nx_) - 1);
            const int j1 = std::min(j0 + 1, static_cast<int>(spec_.ny_) - 1);
            const int k1 = std::min(k0 + 1, static_cast<int>(spec_.nz_) - 1);
            const Real fx = rel[0] - static_cast<Real>(i0);
            const Real fy = rel[1] - static_cast<Real>(j0);
            const Real fz = rel[2] - static_cast<Real>(k0);
            const Real contrib = value * vol;
            const Real w000 = (1.0 - fx) * (1.0 - fy) * (1.0 - fz) * vol;
            const Real w100 = fx * (1.0 - fy) * (1.0 - fz) * vol;
            const Real w010 = (1.0 - fx) * fy * (1.0 - fz) * vol;
            const Real w110 = fx * fy * (1.0 - fz) * vol;
            const Real w001 = (1.0 - fx) * (1.0 - fy) * fz * vol;
            const Real w101 = fx * (1.0 - fy) * fz * vol;
            const Real w011 = (1.0 - fx) * fy * fz * vol;
            const Real w111 = fx * fy * fz * vol;
            const size_t idx000 = index(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k0));
            const size_t idx100 = index(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k0));
            const size_t idx010 = index(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k0));
            const size_t idx110 = index(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k0));
            const size_t idx001 = index(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k1));
            const size_t idx101 = index(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k1));
            const size_t idx011 = index(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k1));
            const size_t idx111 = index(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k1));
            q_host_[idx000] += contrib * (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
            q_host_[idx100] += contrib * fx * (1.0 - fy) * (1.0 - fz);
            q_host_[idx010] += contrib * (1.0 - fx) * fy * (1.0 - fz);
            q_host_[idx110] += contrib * fx * fy * (1.0 - fz);
            q_host_[idx001] += contrib * (1.0 - fx) * (1.0 - fy) * fz;
            q_host_[idx101] += contrib * fx * (1.0 - fy) * fz;
            q_host_[idx011] += contrib * (1.0 - fx) * fy * fz;
            q_host_[idx111] += contrib * fx * fy * fz;
            weight_host_[idx000] += w000;
            weight_host_[idx100] += w100;
            weight_host_[idx010] += w010;
            weight_host_[idx110] += w110;
            weight_host_[idx001] += w001;
            weight_host_[idx101] += w101;
            weight_host_[idx011] += w011;
            weight_host_[idx111] += w111;
            (void)contrib;
        }
    }

    void finalizeFromAccumulators()
    {
        const size_t n = nodeCount();
        for (size_t g = 0; g < n; ++g)
        {
            if (weight_host_[g] > TinyReal)
            {
                q_host_[g] /= weight_host_[g];
            }
            else
            {
                q_host_[g] = Real(0);
            }
        }
    }

    Real sampleHost(const Vecd &x) const
    {
        if (!inside(x))
        {
            return Real(0);
        }
        const Vecd rel = (x - spec_.lower_bound_) / spec_.spacing_;
        const int i0 = static_cast<int>(std::floor(rel[0]));
        const int j0 = static_cast<int>(std::floor(rel[1]));
        const int k0 = static_cast<int>(std::floor(rel[2]));
        const int i1 = std::min(i0 + 1, static_cast<int>(spec_.nx_) - 1);
        const int j1 = std::min(j0 + 1, static_cast<int>(spec_.ny_) - 1);
        const int k1 = std::min(k0 + 1, static_cast<int>(spec_.nz_) - 1);
        const Real fx = rel[0] - static_cast<Real>(i0);
        const Real fy = rel[1] - static_cast<Real>(j0);
        const Real fz = rel[2] - static_cast<Real>(k0);
        const Real c000 = q_host_[index(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k0))];
        const Real c100 = q_host_[index(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k0))];
        const Real c010 = q_host_[index(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k0))];
        const Real c110 = q_host_[index(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k0))];
        const Real c001 = q_host_[index(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k1))];
        const Real c101 = q_host_[index(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k1))];
        const Real c011 = q_host_[index(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k1))];
        const Real c111 = q_host_[index(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k1))];
        const Real c00 = c000 * (1.0 - fx) + c100 * fx;
        const Real c10 = c010 * (1.0 - fx) + c110 * fx;
        const Real c01 = c001 * (1.0 - fx) + c101 * fx;
        const Real c11 = c011 * (1.0 - fx) + c111 * fx;
        const Real c0 = c00 * (1.0 - fy) + c10 * fy;
        const Real c1 = c01 * (1.0 - fy) + c11 * fy;
        return c0 * (1.0 - fz) + c1 * fz;
    }

    void scaleField(Real scale)
    {
        for (Real &q : q_host_)
        {
            q *= scale;
        }
    }

    void computeFieldStats(Real &q_min, Real &q_max, Real &q_mean, Real &q_vol_weighted_mean,
                           Real &empty_node_fraction, size_t &nonempty_node_count) const
    {
        const size_t n = nodeCount();
        q_min = Real(0);
        q_max = Real(0);
        q_mean = Real(0);
        q_vol_weighted_mean = Real(0);
        nonempty_node_count = 0;
        Real sum_q = Real(0);
        Real sum_qw = Real(0);
        Real sum_w = Real(0);
        bool initialized = false;
        for (size_t g = 0; g < n; ++g)
        {
            const Real q = q_host_[g];
            if (weight_host_[g] <= TinyReal)
            {
                continue;
            }
            ++nonempty_node_count;
            if (!initialized)
            {
                q_min = q_max = q;
                initialized = true;
            }
            else
            {
                q_min = std::min(q_min, q);
                q_max = std::max(q_max, q);
            }
            sum_q += q;
            sum_qw += q * weight_host_[g];
            sum_w += weight_host_[g];
        }
        if (nonempty_node_count > 0)
        {
            q_mean = sum_q / static_cast<Real>(nonempty_node_count);
            q_vol_weighted_mean = sum_qw / (sum_w + TinyReal);
        }
        empty_node_fraction = Real(1) - static_cast<Real>(nonempty_node_count) / static_cast<Real>(n + TinyReal);
    }
};

struct Rh200JouleHeatEmGridBundle
{
    Rh200ScalarEulerianGrid joule_heat;
    Rh200ScalarEulerianGrid phi_real;
    Rh200ScalarEulerianGrid phi_imag;
    Rh200ScalarEulerianGrid e_magnitude;
    Rh200ScalarEulerianGrid j_real_x;
    Rh200ScalarEulerianGrid j_real_y;
    Rh200ScalarEulerianGrid j_real_z;
    Rh200ScalarEulerianGrid j_imag_x;
    Rh200ScalarEulerianGrid j_imag_y;
    Rh200ScalarEulerianGrid j_imag_z;
};

inline Real hostSamplePowerFromScalarGrid(const Rh200ScalarEulerianGrid &grid, const StdVec<Vecd> &position,
                                          const StdVec<Real> &volumetric_measure)
{
    Real power = Real(0);
    const size_t n = position.size();
    for (size_t i = 0; i < n; ++i)
    {
        power += grid.sampleHost(position[i]) * volumetric_measure[i];
    }
    return power;
}

inline void hostCompareSampleToEmField(const Rh200ScalarEulerianGrid &grid, const Rh200EmParticleDepositionHost &em,
                                       Real &rel_l2, Real &rel_max)
{
    Real sum_sq = Real(0);
    Real sum_em_sq = Real(0);
    Real max_abs = Real(0);
    Real max_em = TinyReal;
    const size_t n = em.position.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Real q_em = em.joule_heat[i];
        const Real q_sample = grid.sampleHost(em.position[i]);
        const Real diff = q_sample - q_em;
        sum_sq += diff * diff * em.volumetric_measure[i];
        sum_em_sq += q_em * q_em * em.volumetric_measure[i];
        max_abs = std::max(max_abs, std::abs(diff));
        max_em = std::max(max_em, std::abs(q_em));
    }
    rel_l2 = std::sqrt(sum_sq / (sum_em_sq + TinyReal));
    rel_max = max_abs / (max_em + TinyReal);
}

inline Rh200JouleHeatGridBuildReport buildRh200JouleHeatEmGridBundle(
    Rh200JouleHeatEmGridBundle &bundle, const Rh200JouleHeatGridSpec &spec, const Rh200EmParticleDepositionHost &em,
    Real target_power_w)
{
    bundle.joule_heat.spec_ = spec;
    bundle.phi_real.spec_ = spec;
    bundle.phi_imag.spec_ = spec;
    bundle.e_magnitude.spec_ = spec;
    bundle.j_real_x.spec_ = spec;
    bundle.j_real_y.spec_ = spec;
    bundle.j_real_z.spec_ = spec;
    bundle.j_imag_x.spec_ = spec;
    bundle.j_imag_y.spec_ = spec;
    bundle.j_imag_z.spec_ = spec;
    bundle.joule_heat.resetAccumulators();
    bundle.phi_real.resetAccumulators();
    bundle.phi_imag.resetAccumulators();
    bundle.e_magnitude.resetAccumulators();
    bundle.j_real_x.resetAccumulators();
    bundle.j_real_y.resetAccumulators();
    bundle.j_real_z.resetAccumulators();
    bundle.j_imag_x.resetAccumulators();
    bundle.j_imag_y.resetAccumulators();
    bundle.j_imag_z.resetAccumulators();

    StdVec<Real> e_mag(em.position.size(), Real(0));
    StdVec<Real> j_real_x(em.position.size(), Real(0));
    StdVec<Real> j_real_y(em.position.size(), Real(0));
    StdVec<Real> j_real_z(em.position.size(), Real(0));
    StdVec<Real> j_imag_x(em.position.size(), Real(0));
    StdVec<Real> j_imag_y(em.position.size(), Real(0));
    StdVec<Real> j_imag_z(em.position.size(), Real(0));
    for (size_t i = 0; i < em.position.size(); ++i)
    {
        const Vecd er = em.e_real[i];
        const Vecd ei = em.e_imag[i];
        e_mag[i] = std::sqrt(er.dot(er) + ei.dot(ei));
        j_real_x[i] = em.j_real[i][0];
        j_real_y[i] = em.j_real[i][1];
        j_real_z[i] = em.j_real[i][2];
        j_imag_x[i] = em.j_imag[i][0];
        j_imag_y[i] = em.j_imag[i][1];
        j_imag_z[i] = em.j_imag[i][2];
    }

    bundle.joule_heat.depositScalarCloudInCell(em.position, em.joule_heat, em.volumetric_measure);
    bundle.phi_real.depositScalarCloudInCell(em.position, em.phi_real, em.volumetric_measure);
    bundle.phi_imag.depositScalarCloudInCell(em.position, em.phi_imag, em.volumetric_measure);
    bundle.e_magnitude.depositScalarCloudInCell(em.position, e_mag, em.volumetric_measure);
    bundle.j_real_x.depositScalarCloudInCell(em.position, j_real_x, em.volumetric_measure);
    bundle.j_real_y.depositScalarCloudInCell(em.position, j_real_y, em.volumetric_measure);
    bundle.j_real_z.depositScalarCloudInCell(em.position, j_real_z, em.volumetric_measure);
    bundle.j_imag_x.depositScalarCloudInCell(em.position, j_imag_x, em.volumetric_measure);
    bundle.j_imag_y.depositScalarCloudInCell(em.position, j_imag_y, em.volumetric_measure);
    bundle.j_imag_z.depositScalarCloudInCell(em.position, j_imag_z, em.volumetric_measure);

    bundle.joule_heat.finalizeFromAccumulators();
    bundle.phi_real.finalizeFromAccumulators();
    bundle.phi_imag.finalizeFromAccumulators();
    bundle.e_magnitude.finalizeFromAccumulators();
    bundle.j_real_x.finalizeFromAccumulators();
    bundle.j_real_y.finalizeFromAccumulators();
    bundle.j_real_z.finalizeFromAccumulators();
    bundle.j_imag_x.finalizeFromAccumulators();
    bundle.j_imag_y.finalizeFromAccumulators();
    bundle.j_imag_z.finalizeFromAccumulators();

    Rh200JouleHeatGridBuildReport report;
    Real p_em = Real(0);
    for (size_t i = 0; i < em.position.size(); ++i)
    {
        p_em += em.joule_heat[i] * em.volumetric_measure[i];
    }
    report.p_em_particle = p_em;
    report.p_grid_sample_before_scale = hostSamplePowerFromScalarGrid(bundle.joule_heat, em.position, em.volumetric_measure);
    report.grid_scale_factor = Real(1);
    if (target_power_w > TinyReal && report.p_grid_sample_before_scale > TinyReal)
    {
        report.grid_scale_factor = target_power_w / report.p_grid_sample_before_scale;
        bundle.joule_heat.scaleField(report.grid_scale_factor);
    }
    report.p_grid_sample_after_scale =
        hostSamplePowerFromScalarGrid(bundle.joule_heat, em.position, em.volumetric_measure);
    hostCompareSampleToEmField(bundle.joule_heat, em, report.sample_vs_em_rel_l2, report.sample_vs_em_rel_max);
    bundle.joule_heat.computeFieldStats(report.q_min, report.q_max, report.q_mean, report.q_vol_weighted_mean,
                                        report.empty_node_fraction, report.nonempty_node_count);
    return report;
}

inline void writeRh200JouleHeatGridSummaryCsv(const std::string &path, const Rh200JouleHeatGridSpec &spec,
                                              const Rh200JouleHeatGridBuildReport &report)
{
    std::ofstream out(path);
    out << "grid_nx,grid_ny,grid_nz,grid_spacing,bbox_min_x,bbox_min_y,bbox_min_z,bbox_max_x,bbox_max_y,bbox_max_z,"
        << "P_em_particle,P_grid_sample_initial_before_scale,grid_scale_factor,P_grid_sample_initial_after_scale,"
        << "sample_vs_em_rel_l2,sample_vs_em_rel_max,Q_grid_min,Q_grid_max,Q_grid_mean,Q_grid_volume_weighted_mean,"
        << "empty_node_fraction,nonempty_node_count\n";
    out << spec.nx_ << "," << spec.ny_ << "," << spec.nz_ << "," << spec.spacing_ << "," << spec.lower_bound_[0] << ","
        << spec.lower_bound_[1] << "," << spec.lower_bound_[2] << "," << spec.upper_bound_[0] << ","
        << spec.upper_bound_[1] << "," << spec.upper_bound_[2] << "," << report.p_em_particle << ","
        << report.p_grid_sample_before_scale << "," << report.grid_scale_factor << ","
        << report.p_grid_sample_after_scale << "," << report.sample_vs_em_rel_l2 << "," << report.sample_vs_em_rel_max
        << "," << report.q_min << "," << report.q_max << "," << report.q_mean << "," << report.q_vol_weighted_mean
        << "," << report.empty_node_fraction << "," << report.nonempty_node_count << "\n";
    std::cout << "[rh200] wrote joule heat grid summary: " << path << std::endl;
}

inline void writeRh200ScalarGridVtiAscii(const std::string &path, const Rh200JouleHeatGridSpec &spec,
                                         const StdVec<std::pair<std::string, const Rh200ScalarEulerianGrid *>> &fields)
{
    std::ofstream out(path);
    out << std::setprecision(10);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <ImageData Origin=\"" << spec.lower_bound_[0] << " " << spec.lower_bound_[1] << " "
        << spec.lower_bound_[2] << "\" Spacing=\"" << spec.spacing_ << " " << spec.spacing_ << " " << spec.spacing_
        << "\" WholeExtent=\"0 " << (spec.nx_ - 1) << " 0 " << (spec.ny_ - 1) << " 0 " << (spec.nz_ - 1) << "\">\n";
    out << "    <Piece Extent=\"0 " << (spec.nx_ - 1) << " 0 " << (spec.ny_ - 1) << " 0 " << (spec.nz_ - 1)
        << "\">\n";
    out << "      <PointData>\n";
    for (const auto &field : fields)
    {
        out << "        <DataArray type=\"Float64\" Name=\"" << field.first << "\" format=\"ascii\">\n";
        const Rh200ScalarEulerianGrid &grid = *field.second;
        for (size_t k = 0; k < spec.nz_; ++k)
        {
            for (size_t j = 0; j < spec.ny_; ++j)
            {
                for (size_t i = 0; i < spec.nx_; ++i)
                {
                    out << grid.q_host_[grid.index(i, j, k)] << " ";
                }
            }
        }
        out << "\n        </DataArray>\n";
    }
    out << "      </PointData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";
    std::cout << "[rh200] wrote Eulerian grid VTI: " << path << " (" << spec.nx_ << "x" << spec.ny_ << "x" << spec.nz_
              << ")" << std::endl;
}

inline void writeRh200JouleHeatEmGridVti(const std::string &path, const Rh200JouleHeatEmGridBundle &bundle)
{
    writeRh200ScalarGridVtiAscii(
        path, bundle.joule_heat.spec_,
        {{"JouleHeat_Wpm3", &bundle.joule_heat},
         {"PhiReal", &bundle.phi_real},
         {"PhiImag", &bundle.phi_imag},
         {"EMagnitude", &bundle.e_magnitude},
         {"JRealX", &bundle.j_real_x},
         {"JRealY", &bundle.j_real_y},
         {"JRealZ", &bundle.j_real_z},
         {"JImagX", &bundle.j_imag_x},
         {"JImagY", &bundle.j_imag_y},
         {"JImagZ", &bundle.j_imag_z}});
}

#if SPHINXSYS_USE_SYCL

class Rh200JouleHeatGridDevice
{
  public:
    ~Rh200JouleHeatGridDevice() { release(); }

    void upload(const Rh200ScalarEulerianGrid &grid)
    {
        release();
        spec_ = grid.spec_;
        n_total_ = grid.nodeCount();
        origin_ = spec_.lower_bound_;
        spacing_ = spec_.spacing_;
        nx_ = spec_.nx_;
        ny_ = spec_.ny_;
        nz_ = spec_.nz_;
        q_device_ = allocateDeviceShared<Real>(n_total_);
        copyToDevice(grid.q_host_.data(), q_device_, n_total_);
    }

    bool insideDevice(const Vecd &x) const
    {
        const Vecd rel = (x - origin_) / spacing_;
        return rel[0] >= Real(0) && rel[0] <= static_cast<Real>(nx_ - 1) && rel[1] >= Real(0) &&
               rel[1] <= static_cast<Real>(ny_ - 1) && rel[2] >= Real(0) && rel[2] <= static_cast<Real>(nz_ - 1);
    }

    Real *deviceData() const { return q_device_; }
    const Rh200JouleHeatGridSpec &spec() const { return spec_; }

  private:
    void release()
    {
        if (q_device_ != nullptr)
        {
            freeDeviceData(q_device_);
            q_device_ = nullptr;
        }
    }

    Rh200JouleHeatGridSpec spec_;
    Real *q_device_ = nullptr;
    size_t n_total_ = 0;
    Vecd origin_ = Vecd::Zero();
    Real spacing_ = Real(1);
    size_t nx_ = 0;
    size_t ny_ = 0;
    size_t nz_ = 0;
};

inline Real deviceTrilinearSampleJouleGrid(const Real *q, size_t nx, size_t ny, size_t nz, const Vecd &origin,
                                             Real spacing, const Vecd &x)
{
    const Vecd rel = (x - origin) / spacing;
    if (rel[0] < Real(0) || rel[0] > static_cast<Real>(nx - 1) || rel[1] < Real(0) ||
        rel[1] > static_cast<Real>(ny - 1) || rel[2] < Real(0) || rel[2] > static_cast<Real>(nz - 1))
    {
        return Real(0);
    }
    const int i0 = static_cast<int>(sycl::floor(rel[0]));
    const int j0 = static_cast<int>(sycl::floor(rel[1]));
    const int k0 = static_cast<int>(sycl::floor(rel[2]));
    const int i1 = sycl::min(i0 + 1, static_cast<int>(nx) - 1);
    const int j1 = sycl::min(j0 + 1, static_cast<int>(ny) - 1);
    const int k1 = sycl::min(k0 + 1, static_cast<int>(nz) - 1);
    const Real fx = rel[0] - static_cast<Real>(i0);
    const Real fy = rel[1] - static_cast<Real>(j0);
    const Real fz = rel[2] - static_cast<Real>(k0);
    const auto idx = [=](size_t i, size_t j, size_t k) { return i + nx * (j + ny * k); };
    const Real c000 = q[idx(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k0))];
    const Real c100 = q[idx(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k0))];
    const Real c010 = q[idx(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k0))];
    const Real c110 = q[idx(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k0))];
    const Real c001 = q[idx(static_cast<size_t>(i0), static_cast<size_t>(j0), static_cast<size_t>(k1))];
    const Real c101 = q[idx(static_cast<size_t>(i1), static_cast<size_t>(j0), static_cast<size_t>(k1))];
    const Real c011 = q[idx(static_cast<size_t>(i0), static_cast<size_t>(j1), static_cast<size_t>(k1))];
    const Real c111 = q[idx(static_cast<size_t>(i1), static_cast<size_t>(j1), static_cast<size_t>(k1))];
    const Real c00 = c000 * (Real(1) - fx) + c100 * fx;
    const Real c10 = c010 * (Real(1) - fx) + c110 * fx;
    const Real c01 = c001 * (Real(1) - fx) + c101 * fx;
    const Real c11 = c011 * (Real(1) - fx) + c111 * fx;
    const Real c0 = c00 * (Real(1) - fy) + c10 * fy;
    const Real c1 = c01 * (Real(1) - fy) + c11 * fy;
    return c0 * (Real(1) - fz) + c1 * fz;
}

class SampleJouleHeatFromGridCK : public LocalDynamics
{
  public:
    SampleJouleHeatFromGridCK(SPHBody &sph_body, const Rh200JouleHeatGridDevice &grid_device)
        : LocalDynamics(sph_body), grid_device_(&grid_device),
          dv_joule_(particles_->template getVariableByName<Real>(kJouleHeatField)),
          dv_pos_(particles_->template getVariableByName<Vecd>("Position"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : q_(encloser.grid_device_->deviceData()), origin_(encloser.grid_device_->spec().lower_bound_),
              spacing_(encloser.grid_device_->spec().spacing_), nx_(encloser.grid_device_->spec().nx_),
              ny_(encloser.grid_device_->spec().ny_), nz_(encloser.grid_device_->spec().nz_),
              joule_(encloser.dv_joule_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            joule_[index_i] =
                deviceTrilinearSampleJouleGrid(q_, nx_, ny_, nz_, origin_, spacing_, pos_[index_i]);
        }

      protected:
        const Real *q_;
        Vecd origin_;
        Real spacing_;
        size_t nx_;
        size_t ny_;
        size_t nz_;
        Real *joule_;
        Vecd *pos_;
    };

  protected:
    const Rh200JouleHeatGridDevice *grid_device_;
    DiscreteVariable<Real> *dv_joule_;
    DiscreteVariable<Vecd> *dv_pos_;
};

class OutOfGridJouleSampleReduceCK : public LocalDynamicsReduce<ReduceSum<Real>>
{
  public:
    OutOfGridJouleSampleReduceCK(SPHBody &sph_body, const Rh200JouleHeatGridDevice &grid_device)
        : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), grid_device_(&grid_device),
          dv_pos_(particles_->template getVariableByName<Vecd>("Position"))
    {
        quantity_name_ = "Rh200OutOfGridJouleSampleCount";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : origin_(encloser.grid_device_->spec().lower_bound_), spacing_(encloser.grid_device_->spec().spacing_),
              nx_(encloser.grid_device_->spec().nx_), ny_(encloser.grid_device_->spec().ny_),
              nz_(encloser.grid_device_->spec().nz_), pos_(encloser.dv_pos_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd rel = (pos_[index_i] - origin_) / spacing_;
            const bool inside = rel[0] >= Real(0) && rel[0] <= static_cast<Real>(nx_ - 1) && rel[1] >= Real(0) &&
                                rel[1] <= static_cast<Real>(ny_ - 1) && rel[2] >= Real(0) &&
                                rel[2] <= static_cast<Real>(nz_ - 1);
            return inside ? Real(0) : Real(1);
        }

      protected:
        Vecd origin_;
        Real spacing_;
        size_t nx_;
        size_t ny_;
        size_t nz_;
        Vecd *pos_;
    };

  protected:
    const Rh200JouleHeatGridDevice *grid_device_;
    DiscreteVariable<Vecd> *dv_pos_;
};

class SampleCurrentFromGridCK : public LocalDynamics
{
  public:
    SampleCurrentFromGridCK(SPHBody &sph_body, const Rh200JouleHeatGridDevice &j_real_x_device,
                            const Rh200JouleHeatGridDevice &j_real_y_device, const Rh200JouleHeatGridDevice &j_real_z_device,
                            const Rh200JouleHeatGridDevice &j_imag_x_device, const Rh200JouleHeatGridDevice &j_imag_y_device,
                            const Rh200JouleHeatGridDevice &j_imag_z_device, const std::string &j_real_field,
                            const std::string &j_imag_field)
        : LocalDynamics(sph_body), j_real_x_device_(&j_real_x_device), j_real_y_device_(&j_real_y_device),
          j_real_z_device_(&j_real_z_device), j_imag_x_device_(&j_imag_x_device), j_imag_y_device_(&j_imag_y_device),
          j_imag_z_device_(&j_imag_z_device), dv_j_real_(particles_->template getVariableByName<Vecd>(j_real_field)),
          dv_j_imag_(particles_->template getVariableByName<Vecd>(j_imag_field)),
          dv_pos_(particles_->template getVariableByName<Vecd>("Position"))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : jrx_(encloser.j_real_x_device_->deviceData()), jry_(encloser.j_real_y_device_->deviceData()),
              jrz_(encloser.j_real_z_device_->deviceData()), jix_(encloser.j_imag_x_device_->deviceData()),
              jiy_(encloser.j_imag_y_device_->deviceData()), jiz_(encloser.j_imag_z_device_->deviceData()),
              origin_(encloser.j_real_x_device_->spec().lower_bound_), spacing_(encloser.j_real_x_device_->spec().spacing_),
              nx_(encloser.j_real_x_device_->spec().nx_), ny_(encloser.j_real_x_device_->spec().ny_),
              nz_(encloser.j_real_x_device_->spec().nz_), j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)),
              j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)), pos_(encloser.dv_pos_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd x = pos_[index_i];
            j_real_[index_i] = Vecd(deviceTrilinearSampleJouleGrid(jrx_, nx_, ny_, nz_, origin_, spacing_, x),
                                    deviceTrilinearSampleJouleGrid(jry_, nx_, ny_, nz_, origin_, spacing_, x),
                                    deviceTrilinearSampleJouleGrid(jrz_, nx_, ny_, nz_, origin_, spacing_, x));
            j_imag_[index_i] = Vecd(deviceTrilinearSampleJouleGrid(jix_, nx_, ny_, nz_, origin_, spacing_, x),
                                    deviceTrilinearSampleJouleGrid(jiy_, nx_, ny_, nz_, origin_, spacing_, x),
                                    deviceTrilinearSampleJouleGrid(jiz_, nx_, ny_, nz_, origin_, spacing_, x));
        }

      protected:
        const Real *jrx_;
        const Real *jry_;
        const Real *jrz_;
        const Real *jix_;
        const Real *jiy_;
        const Real *jiz_;
        Vecd origin_;
        Real spacing_;
        size_t nx_;
        size_t ny_;
        size_t nz_;
        Vecd *j_real_;
        Vecd *j_imag_;
        Vecd *pos_;
    };

  protected:
    const Rh200JouleHeatGridDevice *j_real_x_device_;
    const Rh200JouleHeatGridDevice *j_real_y_device_;
    const Rh200JouleHeatGridDevice *j_real_z_device_;
    const Rh200JouleHeatGridDevice *j_imag_x_device_;
    const Rh200JouleHeatGridDevice *j_imag_y_device_;
    const Rh200JouleHeatGridDevice *j_imag_z_device_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Vecd> *dv_pos_;
};

#endif // SPHINXSYS_USE_SYCL

inline void writeEnergyBudgetGridCsvHeader(const std::string &path)
{
    std::ofstream out(path, std::ios::out);
    out << "time,P_joule_W,thermal_energy_J,E_joule_integrated_J,T_mean_K,T_max_K,P_wall_loss_W,"
        << "out_of_grid_particle_count,out_of_grid_particle_fraction\n";
}

inline void appendEnergyBudgetGridCsv(const std::string &path, Real time, Real p_joule, Real thermal_energy,
                                      Real e_joule_integrated, Real t_mean, Real t_max, Real out_of_grid_count,
                                      Real out_of_grid_fraction)
{
    std::ofstream out(path, std::ios::app);
    out << time << "," << p_joule << "," << thermal_energy << "," << e_joule_integrated << "," << t_mean << ","
        << t_max << ",0," << out_of_grid_count << "," << out_of_grid_fraction << "\n";
}

} // namespace rh200
} // namespace SPH

#endif // RH200_JOULE_HEAT_GRID_H
