
#include "boundary_face.h"

namespace SPH
{
    SegmentFace::SegmentFace(StdVec<Vecd> boundary_points, Vecd direction)
        : boundary_points_(boundary_points),
          direction_(direction.normalize())
    {
        center_ = Vecd(0);
        for (size_t i = 0; i < boundary_points_.size(); i++)
        {
            center_ += boundary_points_[i];
        }
        center_ /= (Real)boundary_points_.size();
    }

    SegmentFace::SegmentFace(StdVec<Vecd> boundary_points, Vecd direction, Vecd center)
        : boundary_points_(boundary_points),
          direction_(direction.normalize()),
          center_(center)
    {
    }

    SegmentFace::SegmentFace(const Vecd &center, const Vecd &direction)
        : center_(center), direction_(direction.normalize())
    {
    }

    BodyRegionWithFace::BodyRegionWithFace(RealBody &real_body, SegmentFace &segment_face, Real scale)
        : segment_face_(segment_face)
    {
        Real h_ratio_min = 1.0 / real_body.sph_adaptation_->MaximumSpacingRatio();
        region_width_ = scale * real_body.sph_adaptation_->getKernel()->CutOffRadius(h_ratio_min);

        calculate_region_bounds();
    }

    Real BodyRegionWithFace::getSignedDistance(const Vecd &point) const
    {
        Real d = dot((point - segment_face_.center()), segment_face_.direction());
        return d;
    }

    void BodyRegionWithFace::calculate_region_bounds()
    {
        const StdVec<Vecd> &face_points = segment_face_.getBoundaryPoints();
        const Vecd &face_direction = segment_face_.direction();
        StdVec<Vecd> region_bounding_points = face_points;
        if (std::abs(region_width_) > 1e-6)
            for (size_t i = 0; i < face_points.size(); ++i)
            {
                Vecd new_point = face_points[i] + region_width_ * face_direction;
                region_bounding_points.emplace_back(new_point);
            }

        region_bounds.first = Vecd(Infinity);
        region_bounds.second = Vecd(-Infinity);
        size_t dim = region_bounds.first.size();
        for (const auto &point : region_bounding_points)
        {
            for (size_t j = 0; j < dim; ++j)
            {
                region_bounds.first[j] = SMIN(region_bounds.first[j], point[j]);
                region_bounds.second[j] = SMAX(region_bounds.second[j], point[j]);
            }
        }
    }

    BodyRegionByCellsWithFace::BodyRegionByCellsWithFace(RealBody &real_body, SegmentFace &segment_face, Real scale)
        : BodyRegionWithFace(real_body, segment_face, scale),
          real_body_(real_body)
    {
    }

    PartDynamicsByCellsWithFace::PartDynamicsByCellsWithFace(RealBody &real_body, BodyRegionByCellsWithFace &body_region)
        : ParticleDynamics<void>(real_body),
          body_region_(body_region)
    {
        body_region.tagBodyDomainBoundingCells(bound_cells_);
    }

    void PartDynamicsByCellsWithFace::exec(Real dt)
    {
        setBodyUpdated();
        setupDynamics(dt);
        for (size_t i = 0; i != bound_cells_.size(); ++i)
        {
            ListDataVector &cell_list_data = bound_cells_[i]->cell_list_data_;
            for (size_t num = 0; num < cell_list_data.size(); ++num)
            {
                size_t particle_index = cell_list_data[num].first;
                particle_functor_(particle_index, dt);
            }
        }
    }

    PartSimpleDynamicsByCellsWithFace::PartSimpleDynamicsByCellsWithFace(RealBody &real_body, BodyRegionByCellsWithFace &body_region)
        : PartDynamicsByCellsWithFace(real_body, body_region)
    {
        particle_functor_ = std::bind(&PartSimpleDynamicsByCellsWithFace::Update, this, _1, _2);
    }

} // namespace SPH
