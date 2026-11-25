/**
 * @file spdb_temperature_program.h
 * @brief Uniform programmed temperature T(t) and demo nonelastic gate.
 * @details Extra-test only. Last 600 s hold at 950 C is a stated demo assumption.
 */
#ifndef SPDB_TEMPERATURE_PROGRAM_H
#define SPDB_TEMPERATURE_PROGRAM_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spdb
{

inline double celsius_to_kelvin(double t_c) { return t_c + 273.15; }
inline double kelvin_to_celsius(double t_k) { return t_k - 273.15; }

/** Demo gate: closed <=850 C, smoothstep to 900 C, open above. Not a phase model. */
inline double nonelastic_gate_weight_C(double temperature_C, double t_off_C = 850.0, double t_on_C = 900.0)
{
    if (!(t_on_C > t_off_C))
        return temperature_C >= t_on_C ? 1.0 : 0.0;
    if (temperature_C <= t_off_C)
        return 0.0;
    if (temperature_C >= t_on_C)
        return 1.0;
    const double s = (temperature_C - t_off_C) / (t_on_C - t_off_C);
    return s * s * (3.0 - 2.0 * s);
}

inline double nonelastic_gate_weight_K(double temperature_K, double t_off_C = 850.0, double t_on_C = 900.0)
{
    return nonelastic_gate_weight_C(kelvin_to_celsius(temperature_K), t_off_C, t_on_C);
}

class TemperatureProgram
{
  public:
    void clear() { knots_.clear(); }

    void add_knot(double time_s, double temperature_C)
    {
        knots_.emplace_back(time_s, temperature_C);
    }

    void load_csv(const std::string &path)
    {
        clear();
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("TemperatureProgram: cannot open " + path);
        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty())
                continue;
            if (line.find("time") != std::string::npos || line.find("Time") != std::string::npos)
                continue;
            for (char &c : line)
                if (c == ',')
                    c = ' ';
            std::istringstream iss(line);
            double t = 0.0;
            double tc = 0.0;
            if (!(iss >> t >> tc))
                continue;
            add_knot(t, tc);
        }
        if (knots_.size() < 2)
            throw std::runtime_error("TemperatureProgram: need >=2 knots in " + path);
        std::sort(knots_.begin(), knots_.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
    }

    double temperature_C(double time_s) const
    {
        if (knots_.empty())
            throw std::runtime_error("TemperatureProgram: empty");
        if (time_s <= knots_.front().first)
            return knots_.front().second;
        if (time_s >= knots_.back().first)
            return knots_.back().second;
        for (size_t i = 1; i < knots_.size(); ++i)
        {
            if (time_s <= knots_[i].first)
            {
                const double t0 = knots_[i - 1].first;
                const double t1 = knots_[i].first;
                const double y0 = knots_[i - 1].second;
                const double y1 = knots_[i].second;
                const double a = (time_s - t0) / (t1 - t0);
                return y0 + a * (y1 - y0);
            }
        }
        return knots_.back().second;
    }

    double temperature_K(double time_s) const { return celsius_to_kelvin(temperature_C(time_s)); }

    double end_time_s() const { return knots_.empty() ? 0.0 : knots_.back().first; }

    const std::vector<std::pair<double, double>> &knots() const { return knots_; }

    /** Times where T_C crosses level within a segment (for macro-step landing). */
    std::vector<double> crossing_times_C(double level_C) const
    {
        std::vector<double> out;
        for (size_t i = 1; i < knots_.size(); ++i)
        {
            const double y0 = knots_[i - 1].second;
            const double y1 = knots_[i].second;
            if ((y0 - level_C) * (y1 - level_C) > 0.0)
                continue;
            if (std::abs(y1 - y0) < 1.0e-15)
            {
                if (std::abs(y0 - level_C) < 1.0e-12)
                    out.push_back(knots_[i - 1].first);
                continue;
            }
            const double a = (level_C - y0) / (y1 - y0);
            if (a >= 0.0 && a <= 1.0)
                out.push_back(knots_[i - 1].first + a * (knots_[i].first - knots_[i - 1].first));
        }
        return out;
    }

  private:
    std::vector<std::pair<double, double>> knots_;
};

inline bool temperature_program_unit_ok(const TemperatureProgram &prog, double tol_C = 1.0e-9)
{
    const double checks[][2] = {{0, 20},   {150, 160}, {300, 300}, {450, 400},
                                {600, 500}, {900, 700}, {1200, 950}, {1500, 950}, {1800, 950}};
    for (const auto &c : checks)
    {
        if (std::abs(prog.temperature_C(c[0]) - c[1]) > tol_C)
            return false;
    }
    return true;
}

} // namespace spdb

#endif
