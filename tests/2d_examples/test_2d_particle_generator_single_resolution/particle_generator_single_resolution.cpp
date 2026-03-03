/**
* @file 	particle_generator_single_resolution.cpp
* @brief 	This is the test of using level set to generate particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by level set for a complex geometry (2D).
*			Before particle generation, we clean the level set, then do re-initialization.

* @author 	Yongchuan Yu and Xiangyu Hu
*/

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string input_body = "./input/SPHinXsys-2d.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Real DL = 2.5;                          /**< InputBody length right part. */
Real DL1 = 2.5;                         /**< InputBody length left part. */
Real DH = 5.0;                          /**< InputBody height. */
Real global_resolution = (DL + DL1) / 120; /**< Reference resolution. */
BoundingBoxd system_domain_bounds(Vec2d(-DL1, -0.5), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Shape of the InputBody
//----------------------------------------------------------------------
class InputBody : public ComplexShape
{
  public:
    explicit InputBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon original_logo;
        original_logo.addAPolygonFromFile(input_body, ShapeBooleanOps::add);
        add<ExtrudeShape<MultiPolygonShape>>(4.0 * global_resolution, original_logo);
        subtract<MultiPolygonShape>(original_logo);
    }
};
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody input_body(sph_system, makeShared<InputBody>("SPHInXsysLogo"));
    input_body.defineBodyLevelSetShape(2.0)->writeLevelSet();
    input_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation input_body_inner(input_body);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_input_body_particles(input_body);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(input_body_inner);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp input_body_recording_to_vtp(input_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    random_input_body_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    input_body.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    input_body_recording_to_vtp.writeToFile();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            input_body_recording_to_vtp.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;

    return 0;
}
