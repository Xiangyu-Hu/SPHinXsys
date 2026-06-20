

#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_parameters.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * French-paper-inspired reduced cold-crucible case (Jacoutot et al. 2008).
 *
 * Literature-based defaults: D=650 mm, f=300 kHz, sigma=16 S/m @ 1473 K, P_joule ~50 kW example.
 * Remaining dimensions (glass height, coil radius/z/loops, dp, wall) are reduced-case assumptions.
 */

Vec3d glass_center = Vec3d(0.0, 0.0, 0.25);
Real glass_radius = 0.325;
Real glass_half_height = 0.25;
Real dp_0 = 0.02;
int resolution(20);

BoundingBoxd system_domain_bounds(Vec3d(-0.5,-0.5,-0.5), Vec3d(1,1,1));

class GlassCylinderShape : public ComplexShape
{
public:
    GlassCylinderShape(const std::string &shape_name):
    ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0,0,1), glass_radius, glass_half_height, resolution, glass_center);
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody cyliner_glass(sph_system, makeShared<GlassCylinderShape>("CylinderGlassBody"));
    // level set shape is used for particle relaxation
    LevelSetShape &level_set_shape = cyliner_glass.defineBodyLevelSetShape(par_ck)
                                         .correctLevelSetSign()
                                         .writeLevelSet();
    cyliner_glass.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    NearShapeSurface near_body_surface(cyliner_glass);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    Inner<> cyliner_glass_inner(cyliner_glass);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &cyliner_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(cyliner_glass);
    auto &cyliner_glass_update_inner_relation = main_methods.addRelationDynamics(cyliner_glass_inner);
    auto &random_cyliner_glass_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(cyliner_glass);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(cyliner_glass_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(cyliner_glass, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(cyliner_glass);
    auto &update_particle_position =
        main_methods.addStateDynamics<PositionRelaxationCK>(cyliner_glass);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(cyliner_glass);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    random_cyliner_glass_particles.exec();

    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile(0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        cyliner_glass_cell_linked_list.exec();
        cyliner_glass_update_inner_relation.exec();

        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();

        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            body_state_recorder.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;

    return 0;
}



} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH
