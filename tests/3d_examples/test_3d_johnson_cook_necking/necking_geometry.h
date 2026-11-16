#ifndef NECKING_GEOMETRY_H
#define NECKING_GEOMETRY_H

#include "sphinxsys.h"
#include <algorithm>
#include <cmath>

namespace SPH
{
/** Axisymmetric rod with two tangent circular arcs on each side of its waist. */
class NeckingRod : public Shape
{
  public:
    static constexpr Real radius = 8e-4;
    static constexpr Real minimum_radius = 7.5e-4;
    static constexpr Real half_length = 3.2e-3;
    static constexpr Real transition_length = 7.5e-4;
    static constexpr Real grip_position = 2.8e-3;
    static constexpr Real radius_change = radius - minimum_radius;
    static constexpr Real arc_radius = (transition_length * transition_length + radius_change * radius_change) / (4 * radius_change);

    explicit NeckingRod(const std::string &name) : Shape(name) {}
    static Real Radius(Real z)
    {
        z = std::abs(z);
        if (z >= transition_length)
            return radius;
        if (z <= transition_length / 2)
            return minimum_radius + arc_radius - std::sqrt(arc_radius * arc_radius - z * z);
        return radius - arc_radius + std::sqrt(arc_radius * arc_radius - std::pow(transition_length - z, 2));
    }
    bool checkContain(const Vecd &position, bool boundary = true) override
    {
        const Real r = std::hypot(position[0], position[1]);
        return boundary ? r <= Radius(position[2]) && std::abs(position[2]) <= half_length
                        : r < Radius(position[2]) && std::abs(position[2]) < half_length;
    }
    Vecd findClosestPoint(const Vecd &position) override
    {
        const Real r = std::hypot(position[0], position[1]);
        const Real z = std::abs(position[2]);
        const Real angle = std::asin(transition_length / (2 * arc_radius));
        const Real theta_first = std::clamp(std::atan2(z, r - minimum_radius - arc_radius), Pi - angle, Pi);
        const Real theta_second = std::clamp(std::atan2(z - transition_length, r - radius + arc_radius), -angle, Real(0));
        const Vec2d candidates[] = {
            Vec2d(minimum_radius + arc_radius + arc_radius * std::cos(theta_first), arc_radius * std::sin(theta_first)),
            Vec2d(radius - arc_radius + arc_radius * std::cos(theta_second), transition_length + arc_radius * std::sin(theta_second)),
            Vec2d(radius, std::clamp(z, transition_length, half_length)),
            Vec2d(std::clamp(r, Real(0), radius), half_length)};
        Vec2d closest = candidates[0];
        Real distance_squared = (closest - Vec2d(r, z)).squaredNorm();
        for (const Vec2d &candidate : candidates)
        {
            const Real candidate_distance = (candidate - Vec2d(r, z)).squaredNorm();
            if (candidate_distance < distance_squared)
            {
                closest = candidate;
                distance_squared = candidate_distance;
            }
        }
        const Real direction_x = r > 0 ? position[0] / r : Real(1);
        const Real direction_y = r > 0 ? position[1] / r : Real(0);
        return Vecd(closest[0] * direction_x, closest[0] * direction_y, position[2] >= 0 ? closest[1] : -closest[1]);
    }

  protected:
    BoundingBoxd findBounds() override
    {
        return BoundingBoxd(Vecd(-radius, -radius, -half_length), Vecd(radius, radius, half_length));
    }
};

class NeckingLattice;
/** Equal spacing and equal volume throughout the analytic rod. */
template <>
class ParticleGenerator<BaseParticles, NeckingLattice> : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &body, BaseParticles &particles, Real spacing)
        : ParticleGenerator<BaseParticles>(body, particles), spacing_(spacing) {}
    void prepareGeometricData() override
    {
        const int transverse_count = int(std::round(2 * NeckingRod::radius / spacing_));
        const int axial_count = int(std::round(2 * NeckingRod::half_length / spacing_));
        if (std::abs(transverse_count * spacing_ - 2 * NeckingRod::radius) > 1e-10 ||
            std::abs(axial_count * spacing_ - 2 * NeckingRod::half_length) > 1e-10)
            throw std::invalid_argument("Spacing must divide the rod diameter and length to preserve symmetry");
        const Real volume = spacing_ * spacing_ * spacing_;
        for (int k = 0; k != axial_count; ++k)
        {
            const Real z = -NeckingRod::half_length + (k + Real(0.5)) * spacing_;
            const Real radius = NeckingRod::Radius(z);
            for (int j = 0; j != transverse_count; ++j)
            {
                const Real y = -NeckingRod::radius + (j + Real(0.5)) * spacing_;
                for (int i = 0; i != transverse_count; ++i)
                {
                    const Real x = -NeckingRod::radius + (i + Real(0.5)) * spacing_;
                    if (x * x + y * y < radius * radius)
                        addPositionAndVolumetricMeasure(Vecd(x, y, z), volume);
                }
            }
        }
    }

  protected:
    Real spacing_;
};
} // namespace SPH
#endif // NECKING_GEOMETRY_H
