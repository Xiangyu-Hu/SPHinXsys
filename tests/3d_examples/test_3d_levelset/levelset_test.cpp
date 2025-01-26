/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar.
 * @author Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 * //TODO: Seems that the wall contact force should be improved.
 */
#include "levelset_test.h" /**< Case setup for this example. */

using namespace SPH;

int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    // sph_system.setRunParticleRelaxation(false);
    // sph_system.setReloadParticles(true);
    // sph_system.setGenerateRegressionData(false);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    // IOEnvironment io_environment(sph_system);

    /** create a body with corresponding material, particles and reaction model. */
    SolidBody column(sph_system, makeShared<Column>("Column"));
    column.defineAdaptationRatios(1.3, 1.0);
    BoundingBox bounding_box = column.getInitialShape().getBounds();
    LevelSetShape *target = column.defineBodyLevelSetShape();

    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<Real> dist(-1.0, 1.0);

    // Generate a random Vec3d
    Vec3d randomVec;
    randomVec << 0.1, 0.1, 0.1;
    // randomVec << dist(gen), dist(gen), dist(gen);

    target->cleanLevelSet(1.0);
    target->correctLevelSetSign(1.0);
    target->findClosestPoint(randomVec);
    target->findLevelSetGradient(randomVec);
    target->computeKernelIntegral(randomVec);
    target->computeKernelGradientIntegral(randomVec);

    return 0;
}
