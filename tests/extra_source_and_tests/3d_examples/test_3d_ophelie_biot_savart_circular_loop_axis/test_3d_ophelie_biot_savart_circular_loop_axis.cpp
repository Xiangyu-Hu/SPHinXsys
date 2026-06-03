/**
 * @file test_3d_ophelie_biot_savart_circular_loop_axis.cpp
 * @brief Filament circular loop Bz on symmetry axis vs analytic formula.
 */
#include "electromagnetic_ophelie_multiloop_source.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;

int main(int, char *[])
{
    const Real mu0 = 4.0 * Pi * 1.0e-7;
    const Real radius = 0.10;
    const Real current = 2.0;
    const Real z_loop = 0.0;
    const std::vector<Real> z_samples = {0.0, 0.05, 0.10, 0.20, 0.35};

    OphelieCircularLoopSpec loop;
    loop.center = Vecd(0.0, 0.0, z_loop);
    loop.radius = radius;
    loop.current = current;
    loop.segments = 720;

    StdVec<OphelieCurrentMomentSample> moments;
    buildCircularLoopFilamentMoments(loop, moments);

    Real max_rel_err = 0.0;
    for (const Real z : z_samples)
    {
        const Vecd observer(0.0, 0.0, z);
        const Real bz_filament = filamentBzOnAxisFromMoments(moments, observer, mu0);
        const Real bz_analytic = analyticCircularLoopBzOnAxis(z, z_loop, radius, current, mu0);
        const Real rel_err = std::abs(bz_filament - bz_analytic) / (std::abs(bz_analytic) + TinyReal);
        max_rel_err = std::max(max_rel_err, rel_err);
        std::cout << "circular_loop_axis z=" << z << " Bz_filament=" << bz_filament << " Bz_analytic=" << bz_analytic
                  << " rel_err=" << rel_err << std::endl;
    }

    const bool passed = max_rel_err < 0.01;
    std::cout << "test_3d_ophelie_biot_savart_circular_loop_axis max_rel_err=" << max_rel_err
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
