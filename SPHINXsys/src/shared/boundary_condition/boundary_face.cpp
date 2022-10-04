
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

    SegmentFace::SegmentFace(const Vecd &direction, const Vecd &center)
        : direction_(direction.normalize()), center_(center)
    {
    }
    SegmentFace::SegmentFace(const Vecd& direction, const Vecd& center, Real radius)
    : center_(center), direction_(direction.normalize()), radius_(radius)
    {
        Vecd default_normal(1.0,0.0,0.0);
        // calculate rotation matrix from defualt_normal to direction_
        Mat3d rotationMatrix(0);
        {
            default_normal = default_normal.normalize();
            direction_ = direction_.normalize();
            Real angle = acos(SimTK::dot(default_normal,direction_));
            Vecd p_rotate = cross(default_normal,direction_);
            p_rotate = p_rotate.normalize();
            
            rotationMatrix(0, 0) = cos(angle) + p_rotate[0] * p_rotate[0] * (1 - cos(angle));
            rotationMatrix(0, 1) = p_rotate[0] * p_rotate[1] * (1 - cos(angle) - p_rotate[2] * sin(angle));
            rotationMatrix(0, 2) = p_rotate[1] * sin(angle) + p_rotate[0] * p_rotate[2] * (1 - cos(angle));
        
        
            rotationMatrix(1, 0) = p_rotate[2] * sin(angle) + p_rotate[0] * p_rotate[1] * (1 - cos(angle));
            rotationMatrix(1, 1) = cos(angle) + p_rotate[1] * p_rotate[1] * (1 - cos(angle));
            rotationMatrix(1, 2) = -p_rotate[0] * sin(angle) + p_rotate[1] * p_rotate[2] * (1 - cos(angle));
        
        
            rotationMatrix(2, 0) = -p_rotate[1] * sin(angle) +p_rotate[0] * p_rotate[2] * (1 - cos(angle));
            rotationMatrix(2, 1) = p_rotate[0] * sin(angle) + p_rotate[1] * p_rotate[2] * (1 - cos(angle));
            rotationMatrix(2, 2) = cos(angle) + p_rotate[2] * p_rotate[2] * (1 - cos(angle));
        }
        // calculate rotation angle
        // Euler ZYX rotate ref:http://web.mit.edu/2.05/www/Handout/HO2.PDF
        Mat3d rotationMatrix1(0);
        {
            Real gamma, beta, alpha;
            gamma = atan(rotationMatrix(2,1)/rotationMatrix(2,2)); // gamma
            // // theta_y = atan2(direction_[0]*cos(theta_x),direction_[2]);
            beta = atan(-rotationMatrix(2,0)/sqrt(rotationMatrix(0,0)*rotationMatrix(0,0)+rotationMatrix(1,0)*rotationMatrix(1,0))); //beta
            // // theta_z = atan2(cos(theta_x),sin(theta_x)*sin(theta_y));
            alpha = atan(rotationMatrix(1,0)/rotationMatrix(0,0)); //alpha
            
            rotationMatrix1(0,0) = cos(alpha)*cos(beta);
            rotationMatrix1(0,1) = cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma);
            rotationMatrix1(0,2) = cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma);

            rotationMatrix1(1,0) = sin(alpha)*cos(beta);
            rotationMatrix1(1,1) = sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma);
            rotationMatrix1(1,2) = sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma);

            rotationMatrix1(2,0) = -sin(beta);
            rotationMatrix1(2,1) = cos(beta)*sin(gamma);
            rotationMatrix1(2,2) = cos(beta)*cos(gamma);

        }
        // defualt boundary points 
        Vecd RU1(0., 1, 0), RD1(0., -1, 0), LU1(0., 0, 1), LD1(0., 0, -1);
        // rotate and translate defualt boundary points
        RU1 = rotationMatrix1*RU1*radius_*1+ center_;
        RD1 = rotationMatrix1*RD1*radius_*1+ center_;
        LU1 = rotationMatrix1*LU1*radius_*1+ center_;
        LD1 = rotationMatrix1*LD1*radius_*1+ center_;
        StdVec<SPH::Vecd> inlet_bound_points{RU1,RD1,LU1,LD1};
        
        boundary_points_ = inlet_bound_points;
        
    }

    BodyRegionWithFace::BodyRegionWithFace(RealBody &real_body, SegmentFace &segment_face, Real scale)
        : segment_face_(segment_face),radius_(segment_face.Radius())
    {
        region_width_ = scale * real_body.sph_adaptation_->getKernel()->CutOffRadius();

        calculate_region_bounds();
    }

    Real BodyRegionWithFace::getSignedDistance(const Vecd &point) const
    {
        Real d = dot((point - segment_face_.center()), segment_face_.direction());
        return d;
    }
    /* The radius here represents the distance from the point to the line (consisting of the center and the direction.)
    *  The periodic_translation_ here represents the projected distance from the point to the line (consisting of the center and the direction.)
    *  we need d > periodic_translation_ and r < radius
    *        o(point)
    *       /|
    *      / | r
    *     /  |
    *    /_d_|_____> direction  
    *   o(center)   
    */
    bool BodyRegionWithFace::insertParticle(const Vecd &point, Real periodic_translation) const
    {
        Real d = dot((point - segment_face_.center()), segment_face_.direction());
        if(d<periodic_translation)
            return false;
        Real r2 = dot((point - segment_face_.center()), (point - segment_face_.center())) - d*d;
        if(r2>radius_*radius_)
            return false;
        return true;
    }

    bool BodyRegionWithFace::inDomain(const Vecd &point, Real periodic_translation) const
    {
        Real d = dot((point - segment_face_.center()), segment_face_.direction());
        if(d>periodic_translation)
            return false;
        Real r2 = dot((point - segment_face_.center()), (point - segment_face_.center())) - d*d;
        if(r2>radius_*radius_)
            return false;
        return true;
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
