/**
 * @file 	test_2d_shell_particle_relaxation.cpp
 * @brief 	This is a test for generating 2D shell particles from the middle surface of a 2D thin pipe.
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
/**
 * @brief Basic geometry parameters.
 */
Real radius = 24.5;                                 /** Inner radius of a 2D thin pipe. */
Real thickness = 1.0;                               /** Thickness of the pipe. */
Real radius_mid_surface = radius + thickness / 2.0; /** Radius of the mid circle. */
Real resolution_ref = 0.5;                          // Global reference resolution
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
Vec2d pipe_center(0.0, 0.0); /** Location of the pipe center. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-radius - thickness, -radius - thickness),
                                 Vec2d(radius + thickness, radius + thickness));
/**
 * @brief define geometry of SPH bodies
 */
class Pipe : public MultiPolygonShape
{
  public:
    explicit Pipe(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(pipe_center, radius + thickness, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(pipe_center, radius, 100, ShapeBooleanOps::sub);
    }
};
//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Build up a SPHSystem. */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);

    /** Creating body, materials and particles. */
    SolidBody pipe_body(sph_system, makeShared<Pipe>("PipeBody"));
    pipe_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    pipe_body.defineBodyLevelSetShape(level_set_refinement_ratio, UsageType::Surface)
        ->writeLevelSet(sph_system);
    pipe_body.generateParticles<SurfaceParticles, Lattice>(thickness);

    InnerRelation pipe_body_inner(pipe_body);
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_pipe_body_particles(pipe_body);
    /** A  Physics relaxation step. */
    ShellRelaxationStep relaxation_step_pipe_body_inner(pipe_body_inner);
    ShellNormalDirectionPrediction shell_normal_prediction(pipe_body_inner, thickness);
    /**
     * @brief define simple data file input and outputs functions.
     */
    BodyStatesRecordingToVtp write_real_body_states({pipe_body});
    write_real_body_states.addToWrite<Vecd>(pipe_body, "NormalDirection");
    write_real_body_states.addToWrite<int>(pipe_body, "UpdatedIndicator");
    /**
     * @brief 	Particle relaxation starts here.
     */
    random_pipe_body_particles.exec(0.25);
    relaxation_step_pipe_body_inner.MidSurfaceBounding().exec();
    write_real_body_states.writeToFile(0.0);
    pipe_body.updateCellLinkedList();

    /** relax particles of the insert body. */
    int ite_p = 0;
    while (ite_p < 2000)
    {
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
            write_real_body_states.writeToFile(ite_p);
        }
        relaxation_step_pipe_body_inner.exec();
        ite_p += 1;
    }
    shell_normal_prediction.exec();
    write_real_body_states.writeToFile(ite_p);
    std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

    return 0;
}