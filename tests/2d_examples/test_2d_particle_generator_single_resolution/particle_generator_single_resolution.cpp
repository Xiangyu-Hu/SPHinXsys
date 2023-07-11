/**
* @file 	particle_generator_single_resolution.cpp
* @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
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
Real resolution_ref = (DL + DL1) / 120; /**< Reference resolution. */
BoundingBox system_domain_bounds(Vec2d(-DL1, -0.5), Vec2d(DL, DH));
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
        add<ExtrudeShape<MultiPolygonShape>>(4.0 * resolution_ref, original_logo);
        subtract<MultiPolygonShape>(original_logo);
    }
};
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody input_body(system, makeShared<InputBody>("SPHInXsysLogo"));
    input_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    input_body.defineParticlesAndMaterial();
    input_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation input_body_inner(input_body);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SimpleDynamics<RandomizeParticlePosition> random_input_body_particles(input_body);
    relax_dynamics::RelaxationStepInner relaxation_step_inner(input_body_inner, true);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp input_body_recording_to_vtp(io_environment, input_body);
    MeshRecordingToPlt cell_linked_list_recording(io_environment, input_body.getCellLinkedList());
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
    cell_linked_list_recording.writeToFile();
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
