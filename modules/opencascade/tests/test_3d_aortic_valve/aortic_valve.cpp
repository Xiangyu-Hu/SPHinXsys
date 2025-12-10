/**
 * @file 	test_3d_surface_particle_relaxation.cpp
 * @brief 	This is the test of using OCCT to import surface model,and then generate  particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation  for a surface geometry (3D).
 */

#include "relax_dynamics_surface.h"
#include "sphinxsys.h" // SPHinXsys Library.
#include "surface_shape.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
Standard_CString full_path_to_geometry = "./input/valve.STEP";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real thickness = 1;
Real PL = 13.95;
int particle_number = 10; /** Particle number in the peripheral direction. */
Real particle_spacing_ref = PL / (Real)particle_number;
double DELTA1 = 0.1;
double DELTA2 = 0.05;

Vec3d domain_lower_bound(-15.0, 0.0, 0.0);
Vec3d domain_upper_bound(15.0, 15.0, 26.0);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

namespace SPH
{
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        return index_i >= 0 && index_i <= (2 / DELTA1 + 2 / DELTA2 - 1);
    };
};

class Leaflet;
template <>
class ParticleGenerator<SurfaceParticles, Leaflet> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles), sph_body_(sph_body) {};
    virtual void prepareGeometricData() override
    {
        SurfaceShape *a = DynamicCast<SurfaceShape>(this, &sph_body_.getInitialShape());

        Standard_Real u1 = 0;
        Standard_Real v1 = DELTA1;
        Standard_Real u2 = DELTA2;
        Standard_Real v2 = DELTA2;

        std::vector<Vecd> points;
        for (size_t k = 0; k <= 1 / DELTA1; k++)
        {
            Standard_Real u = u1 + k * DELTA1;
            points.push_back(a->getCartesianPoint(u, 0));
        }
        for (size_t k = 0; k <= 1 / DELTA1 - 1; k++)
        {
            Standard_Real v = v1 + k * DELTA1;
            points.push_back(a->getCartesianPoint(0, v));
        }

        for (size_t n = 0; n <= 1 / DELTA2 - 2; n++)
        {
            Standard_Real u = u2 + n * DELTA2;
            points.push_back(a->getCartesianPoint(u, 1));
        }

        for (size_t n = 0; n <= 1 / DELTA2 - 1; n++)
        {
            Standard_Real v = v2 + n * DELTA2;
            points.push_back(a->getCartesianPoint(1, v));
        }

        for (size_t k = 0; k <= 8; k++)
        {
            v1 = 0.5;
            u1 = 0.5;
            double V_DELTA = 0.005;
            double U_DELTA = 0.005;
            Standard_Real v = v1 + k * V_DELTA;
            for (size_t n = 0; n <= 20; n++)
            {
                Standard_Real u = u1 + n * U_DELTA;
                points.push_back(a->getCartesianPoint(u, v));
            }
        }

        for (int i = 0; i != (int)points.size(); i++)
        {
            addPositionAndVolumetricMeasure(points[i], particle_spacing_ref * particle_spacing_ref);
            Vecd n_0 = Vec3d(1.0, 1.0, 1.0);
            addSurfaceProperties(n_0, thickness);
        }
    }
    SPHBody &sph_body_;
};
} // namespace SPH

//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody leaflet(sph_system, makeShared<SurfaceShapeSTEP>(full_path_to_geometry, "Leaflet"));
    // here dummy linear elastic solid is use because no solid dynamics in particle relaxation
    leaflet.defineMaterial<Solid>();
    leaflet.generateParticles<SurfaceParticles, Leaflet>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_relaxed_particles(sph_system);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation leaflet_inner(leaflet);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    RelaxationStepInnerFirstHalf leaflet_relaxation_first_half(leaflet_inner);
    RelaxationStepInnerSecondHalf leaflet_relaxation_second_half(leaflet_inner);
    /** Constrain the boundary. */
    BoundaryGeometry boundary_geometry(leaflet);
    SimpleDynamics<SurfaceNormalDirection> surface_normal_direction(leaflet);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    write_relaxed_particles.writeToFile(0);
    leaflet.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite = 0;
    int relax_step = 1000;
    while (ite < relax_step)
    {
        leaflet_relaxation_first_half.exec();
        leaflet_relaxation_second_half.exec();
        ite += 1;
        if (ite % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
            write_relaxed_particles.writeToFile(ite);
        }
    }
    surface_normal_direction.exec();
    write_relaxed_particles.writeToFile(ite);
    std::cout << "The physics relaxation process of surface particles finish !" << std::endl;

    return 0;
}
