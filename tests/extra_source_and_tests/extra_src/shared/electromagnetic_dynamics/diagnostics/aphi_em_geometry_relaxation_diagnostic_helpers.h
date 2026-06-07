#ifndef APHI_EM_GEOMETRY_RELAXATION_DIAGNOSTIC_HELPERS_H
#define APHI_EM_GEOMETRY_RELAXATION_DIAGNOSTIC_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"

#include <algorithm>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline Vecd geometryBoxRegionCenter(const benchmark::AphiBoxRegion &region)
{
    return Vecd(0.5 * (region.xmin + region.xmax), 0.5 * (region.ymin + region.ymax),
                0.5 * (region.zmin + region.zmax));
}

inline Vecd geometryBoxRegionHalfSize(const benchmark::AphiBoxRegion &region)
{
    return Vecd(0.5 * (region.xmax - region.xmin), 0.5 * (region.ymax - region.ymin),
                0.5 * (region.zmax - region.zmin));
}

class AphiTeam7GeometryAirSlabsShape : public ComplexShape
{
  public:
    AphiTeam7GeometryAirSlabsShape(const std::string &shape_name, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                   Real body_length, Real body_height, Real body_width)
        : ComplexShape(shape_name)
    {
        const benchmark::AphiBoxRegion left_air{0.0, layout.coil.xmin, 0.0, body_height, 0.0, body_width};
        const benchmark::AphiBoxRegion gap_air{layout.coil.xmax, layout.conductor.xmin, 0.0, body_height, 0.0,
                                               body_width};
        const benchmark::AphiBoxRegion right_air{layout.conductor.xmax, body_length, 0.0, body_height, 0.0, body_width};
        for (const benchmark::AphiBoxRegion &region : {left_air, gap_air, right_air})
        {
            add<GeometricShapeBox>(Transform(geometryBoxRegionCenter(region)), geometryBoxRegionHalfSize(region));
        }
    }
};

class AphiTeam7GeometryRegionBoxShape : public ComplexShape
{
  public:
    AphiTeam7GeometryRegionBoxShape(const std::string &shape_name, const benchmark::AphiBoxRegion &region)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(geometryBoxRegionCenter(region)), geometryBoxRegionHalfSize(region));
    }
};

/** Lightweight three-body lattice scaffold (no EM/CK solver headers). */
struct AphiTeam7GeometryLatticeCase
{
    SPHSystem sph_system;
    benchmark::AphiTeam7LikeUnitBoxLayout layout;
    SolidBody air_body;
    SolidBody coil_body;
    SolidBody plate_body;
    UniquePtr<Inner<>> air_inner_ck;
    UniquePtr<Inner<>> coil_inner_ck;
    UniquePtr<Inner<>> plate_inner_ck;
    UniquePtr<Contact<>> air_contact_ck;
    UniquePtr<Contact<>> coil_to_air_ck;
    UniquePtr<Contact<>> plate_to_air_ck;

    AphiTeam7GeometryLatticeCase(Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
                                 int ac, char *av[])
        : sph_system(BoundingBoxd(Vecd(-boundary_width, -boundary_width, -boundary_width),
                                  Vecd(body_length + boundary_width, body_height + boundary_width,
                                       body_width + boundary_width)),
                     dp_0),
          layout(benchmark::buildTeam7LayoutForBox(body_length, body_height, body_width)),
          air_body(sph_system,
                   makeShared<AphiTeam7GeometryAirSlabsShape>("AirBody", layout, body_length, body_height, body_width)),
          coil_body(sph_system, makeShared<AphiTeam7GeometryRegionBoxShape>("CoilBody", layout.coil)),
          plate_body(sph_system, makeShared<AphiTeam7GeometryRegionBoxShape>("PlateBody", layout.conductor))
    {
        if (ac > 0)
        {
            sph_system.handleCommandlineOptions(ac, av);
        }
        for (auto *body_ptr : {&air_body, &coil_body, &plate_body})
        {
            body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
            body_ptr->defineMaterial<Solid>();
            body_ptr->defineBodyLevelSetShape();
            body_ptr->generateParticles<BaseParticles, Lattice>();
        }
        air_inner_ck = makeUnique<Inner<>>(air_body);
        coil_inner_ck = makeUnique<Inner<>>(coil_body);
        plate_inner_ck = makeUnique<Inner<>>(plate_body);
        air_contact_ck = makeUnique<Contact<>>(air_body, StdVec<RealBody *>{&coil_body, &plate_body});
        coil_to_air_ck = makeUnique<Contact<>>(coil_body, StdVec<RealBody *>{&air_body});
        plate_to_air_ck = makeUnique<Contact<>>(plate_body, StdVec<RealBody *>{&air_body});
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
    }
};

struct AphiEmGeometryRelaxationSpec
{
    Real dp = 0.1;
    Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    Real boundary_width = 0.3;
    int relaxation_steps = 0;
};

struct AphiEmGeometryRelaxationMetrics
{
    size_t particles = 0;
    Real min_body_particle_count = 0.0;
    Real max_body_particle_count = 0.0;
    Real mean_body_particle_count = 0.0;
    bool body_particle_count_finite = true;
};

inline AphiEmGeometryRelaxationMetrics runEmGeometryRelaxationDiagnostic(int ac, char *av[],
                                                                         const AphiEmGeometryRelaxationSpec &spec)
{
    AphiEmGeometryRelaxationMetrics metrics;
    AphiTeam7GeometryLatticeCase case_setup(spec.dp, spec.body_length, spec.body_height, spec.body_width,
                                            spec.boundary_width, ac, av);

    const size_t air_particles = case_setup.air_body.getBaseParticles().TotalRealParticles();
    const size_t coil_particles = case_setup.coil_body.getBaseParticles().TotalRealParticles();
    const size_t plate_particles = case_setup.plate_body.getBaseParticles().TotalRealParticles();
    metrics.particles = air_particles + coil_particles + plate_particles;

    const Real body_counts[] = {static_cast<Real>(air_particles), static_cast<Real>(coil_particles),
                                static_cast<Real>(plate_particles)};
    metrics.min_body_particle_count = *std::min_element(body_counts, body_counts + 3);
    metrics.max_body_particle_count = *std::max_element(body_counts, body_counts + 3);
    metrics.mean_body_particle_count = static_cast<Real>(metrics.particles) / 3.0;
    metrics.body_particle_count_finite =
        std::isfinite(metrics.min_body_particle_count) && std::isfinite(metrics.max_body_particle_count);
    return metrics;
}

inline bool emGeometryRelaxationDiagnosticPassed(const AphiEmGeometryRelaxationMetrics &metrics)
{
    return metrics.particles > 0 && metrics.body_particle_count_finite && metrics.min_body_particle_count >= 8.0 &&
           metrics.max_body_particle_count <= 50000.0;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
