/**
 * @file test_3d_ophelie_biot_savart_multiloop_axis.cpp
 * @brief Multi-loop circular coil Bz on axis = sum of single-loop analytic values.
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
    OphelieMultiloopCoilSpec spec;
    spec.stack_center = Vecd(0.0, 0.0, 0.0);
    spec.loop_radius = 0.12;
    spec.z_min = -0.10;
    spec.z_max = 0.10;
    spec.num_loops = 5;
    spec.current_per_loop = 1.5;
    spec.segments_per_loop = 360;

    StdVec<OphelieCurrentMomentSample> moments;
    buildMultiloopFilamentMoments(spec, moments);

    const std::vector<Real> z_samples = {-0.05, 0.0, 0.05, 0.12, 0.25};
    Real max_rel_err = 0.0;
    for (const Real z : z_samples)
    {
        const Vecd observer(0.0, 0.0, z);
        const Real bz_filament = filamentBzOnAxisFromMoments(moments, observer, mu0);
        const Real bz_analytic = analyticMultiloopBzOnAxis(z, spec, mu0);
        const Real rel_err = std::abs(bz_filament - bz_analytic) / (std::abs(bz_analytic) + TinyReal);
        max_rel_err = std::max(max_rel_err, rel_err);
        std::cout << "multiloop_axis z=" << z << " Bz_filament=" << bz_filament << " Bz_analytic=" << bz_analytic
                  << " rel_err=" << rel_err << std::endl;
    }

    const bool passed = max_rel_err < 0.01;
    std::cout << "test_3d_ophelie_biot_savart_multiloop_axis max_rel_err=" << max_rel_err
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
