#ifndef ELECTROMAGNETIC_OPHELIE_SPATIAL_SIGMA_FIELD_H
#define ELECTROMAGNETIC_OPHELIE_SPATIAL_SIGMA_FIELD_H

/**
 * Doc156: load per-particle electrical conductivity for French EM.
 * Operators already read names.sigma ("Sigma"); this only populates the field.
 * spatial-field mode must hard-fail on missing/incomplete coverage (no silent mean fallback).
 */

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieSpatialSigmaLoadReport
{
    bool ok = false;
    std::string error;
    std::string path;
    std::string match_mode; /**< original_id | coordinate | constant */
    size_t n_particles = 0;
    size_t n_assigned = 0;
    size_t n_unmatched = 0;
    Real sigma_min = 0;
    Real sigma_max = 0;
    Real sigma_mean_vol = 0;
};

inline bool ophelieParseCsvHeader(const std::string &line, std::vector<std::string> &hdr)
{
    hdr.clear();
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ','))
    {
        for (char &c : cell)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        while (!cell.empty() && (cell.front() == ' ' || cell.front() == '\t'))
            cell.erase(cell.begin());
        while (!cell.empty() && (cell.back() == ' ' || cell.back() == '\t' || cell.back() == '\r'))
            cell.pop_back();
        hdr.push_back(cell);
    }
    return !hdr.empty();
}

inline int ophelieCsvCol(const std::vector<std::string> &hdr, std::initializer_list<const char *> names)
{
    for (const char *n : names)
    {
        for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
        {
            if (hdr[static_cast<size_t>(i)] == n)
                return i;
        }
    }
    return -1;
}

/**
 * Load spatial sigma CSV and assign onto glass particles.
 * Preferred columns: OriginalID,sigma_S_m (plus optional x/y/z/vol/T/valid).
 * Falls back to nearest-coordinate match if OriginalID column absent.
 */
inline OphelieSpatialSigmaLoadReport loadAndAssignOphelieSpatialSigmaField(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const std::string &csv_path,
    Real coord_tol_m = Real(1.0e-7))
{
    OphelieSpatialSigmaLoadReport r;
    r.path = csv_path;
    std::ifstream in(csv_path);
    if (!in)
    {
        r.error = "cannot_open_sigma_field_file";
        return r;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        r.error = "empty_sigma_field_file";
        return r;
    }
    std::vector<std::string> hdr;
    if (!ophelieParseCsvHeader(line, hdr))
    {
        r.error = "bad_header";
        return r;
    }
    const int i_oid = ophelieCsvCol(hdr, {"originalid", "original_id", "id"});
    const int i_sig = ophelieCsvCol(hdr, {"sigma_s_m", "sigma", "sigma_sm"});
    const int i_x = ophelieCsvCol(hdr, {"x_m", "x"});
    const int i_y = ophelieCsvCol(hdr, {"y_m", "y"});
    const int i_z = ophelieCsvCol(hdr, {"z_m", "z"});
    const int i_valid = ophelieCsvCol(hdr, {"valid"});
    if (i_sig < 0)
    {
        r.error = "missing_sigma_S_m_column";
        return r;
    }

    struct Row
    {
        long oid = -1;
        Real x = 0, y = 0, z = 0, s = 0;
        bool has_xyz = false;
    };
    std::vector<Row> rows;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(ss, cell, ','))
            cells.push_back(cell);
        auto at = [&](int i) -> std::string {
            if (i < 0 || i >= static_cast<int>(cells.size()))
                return {};
            return cells[static_cast<size_t>(i)];
        };
        if (i_valid >= 0)
        {
            const std::string v = at(i_valid);
            if (v == "0" || v == "false" || v == "False")
                continue;
        }
        Row row;
        row.s = static_cast<Real>(std::atof(at(i_sig).c_str()));
        if (!std::isfinite(static_cast<double>(row.s)) || row.s <= TinyReal)
        {
            r.error = "nonfinite_or_nonpositive_sigma_in_file";
            return r;
        }
        if (i_oid >= 0)
            row.oid = std::atol(at(i_oid).c_str());
        if (i_x >= 0 && i_y >= 0 && i_z >= 0)
        {
            row.x = static_cast<Real>(std::atof(at(i_x).c_str()));
            row.y = static_cast<Real>(std::atof(at(i_y).c_str()));
            row.z = static_cast<Real>(std::atof(at(i_z).c_str()));
            row.has_xyz = true;
        }
        rows.push_back(row);
    }
    if (rows.empty())
    {
        r.error = "no_valid_sigma_rows";
        return r;
    }

    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, names.sigma);
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const size_t n = particles.TotalRealParticles();
    r.n_particles = n;

    std::vector<uint8_t> assigned(n, 0);
    if (i_oid >= 0)
    {
        r.match_mode = "original_id";
        std::unordered_map<long, Real> by_id;
        by_id.reserve(rows.size() * 2);
        for (const Row &row : rows)
            by_id[row.oid] = row.s;
        for (size_t i = 0; i < n; ++i)
        {
            const auto it = by_id.find(static_cast<long>(oid[i]));
            if (it == by_id.end())
            {
                ++r.n_unmatched;
                continue;
            }
            sigma[i] = it->second;
            assigned[i] = 1;
            ++r.n_assigned;
        }
    }
    else if (i_x >= 0)
    {
        r.match_mode = "coordinate";
        for (size_t i = 0; i < n; ++i)
        {
            Real best = std::numeric_limits<Real>::max();
            Real best_s = 0;
            for (const Row &row : rows)
            {
                if (!row.has_xyz)
                    continue;
                const Real d2 = (pos[i][0] - row.x) * (pos[i][0] - row.x) +
                                (pos[i][1] - row.y) * (pos[i][1] - row.y) +
                                (pos[i][2] - row.z) * (pos[i][2] - row.z);
                if (d2 < best)
                {
                    best = d2;
                    best_s = row.s;
                }
            }
            if (best > coord_tol_m * coord_tol_m)
            {
                ++r.n_unmatched;
                continue;
            }
            sigma[i] = best_s;
            assigned[i] = 1;
            ++r.n_assigned;
        }
    }
    else
    {
        r.error = "need_OriginalID_or_xyz_columns";
        return r;
    }

    if (r.n_assigned != n || r.n_unmatched > 0)
    {
        r.error = "incomplete_sigma_coverage";
        std::cerr << "[ophelie][sigma-field] FAIL incomplete coverage assigned=" << r.n_assigned << "/" << n
                  << " unmatched=" << r.n_unmatched << " path=" << csv_path << std::endl;
        return r;
    }

    Real smin = sigma[0], smax = sigma[0], acc = 0, w = 0;
    for (size_t i = 0; i < n; ++i)
    {
        smin = std::min(smin, sigma[i]);
        smax = std::max(smax, sigma[i]);
        acc += sigma[i] * vol[i];
        w += vol[i];
    }
    r.sigma_min = smin;
    r.sigma_max = smax;
    r.sigma_mean_vol = acc / (w + TinyReal);
    syncVariableToDevice<Real>(particles, names.sigma);
    r.ok = true;
    std::cout << "[ophelie][sigma-field] loaded path=" << csv_path << " mode=" << r.match_mode << " n=" << n
              << " sigma=[" << r.sigma_min << "," << r.sigma_max << "] mean_vol=" << r.sigma_mean_vol
              << std::endl;
    return r;
}

inline bool writeOphelieConstantSigmaFieldCsv(const std::string &path, BaseParticles &particles, Real sigma_value)
{
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    std::ofstream out(path);
    if (!out)
        return false;
    out.setf(std::ios::scientific);
    out.precision(10);
    out << "OriginalID,x_m,y_m,z_m,vol_m3,T_mapped_K,sigma_S_m,valid,support_radius_m\n";
    for (size_t i = 0; i < n; ++i)
    {
        out << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << vol[i] << ","
            << 1473.0 << "," << sigma_value << ",1,0\n";
    }
    return true;
}

inline bool writeOphelieEmParticlesLatticeCsv(const std::string &path, BaseParticles &particles)
{
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    std::ofstream out(path);
    if (!out)
        return false;
    out.setf(std::ios::scientific);
    out.precision(16);
    out << "OriginalID,x_m,y_m,z_m,vol_m3\n";
    for (size_t i = 0; i < n; ++i)
        out << oid[i] << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << "," << vol[i] << "\n";
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_SPATIAL_SIGMA_FIELD_H
