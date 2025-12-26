/**
 * @file 	airfoil_2d.cpp
 * @brief 	This is the test of using level set to generate body fitted SPH particles.
 * @details	We use this case to test the particle generation and relaxation with a complex geometry (2D).
 *			Before the particles are generated, we clean the sharp corners and other unresolvable surfaces.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string airfoil_flap_front = "./input/airfoil_flap_front.dat";
std::string airfoil_wing = "./input/airfoil_wing.dat";
std::string airfoil_flap_rear = "./input/airfoil_flap_rear.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.25;                /**< airfoil length rear part. */
Real DL1 = 0.25;               /**< airfoil length front part. */
Real DH = 0.25;                /**< airfoil height. */
Real global_resolution = 0.02; /**< Reference resolution. */
BoundingBoxd system_domain_bounds(Vec2d(-DL1, -DH), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	import model as a complex shape
//----------------------------------------------------------------------
class ImportModel : public MultiPolygonShape
{
  public:
    explicit ImportModel(const std::string &import_model_name) : MultiPolygonShape(import_model_name)
    {
        multi_polygon_.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.setRunParticleRelaxation(true); // tag to run particle relaxation when no commandline option
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &airfoil = sph_system.addAdaptiveBody<RealBody>(
        AdaptiveNearSurface(global_resolution, 1.15, 1.0, 3), makeShared<ImportModel>("AirFoil"));
    LevelSetShape *level_set_shape = airfoil.defineBodyLevelSetShape()
                                         ->cleanLevelSet()
                                         ->addCellVariableToWrite<UnsignedInt>("CellPackageIndex")
                                         ->writeLevelSet();
    airfoil.generateParticles<BaseParticles, Lattice>();
    auto &near_body_surface = airfoil.addBodyPart<NearShapeSurface>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &airfoil_inner = sph_system.addInnerRelation(airfoil);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    host_methods.addStateDynamics<RandomizeParticlePositionCK>(airfoil).exec();
    auto &airfoil_update_cell_linked_list = main_methods.addCellLinkedListDynamics(airfoil);
    auto &airfoil_update_inner_relation = main_methods.addRelationDynamics(airfoil_inner);

    ParticleDynamicsGroup relaxation_residual;
    relaxation_residual.add(
        &main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(airfoil_inner)
             .addPostStateDynamics<LevelsetKernelGradientIntegral>(airfoil, *level_set_shape));

    ParticleDynamicsGroup update_particle_position;
    update_particle_position.add(&main_methods.addStateDynamics<PositionRelaxationCK>(airfoil));
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(near_body_surface));
    update_particle_position.add(&main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(airfoil, *level_set_shape));
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(airfoil, "SmoothingLengthRatio");
    body_state_recorder.addToWrite<UnsignedInt>(airfoil, "NeighborSize");
    BaseCellLinkedList &airfoil_cell_linked_list = airfoil.getCellLinkedList();
    airfoil_cell_linked_list.addCellVariableToWrite<UnsignedInt>("CurrentListSize");
    MeshRecordingToPlt cell_linked_list_recording(sph_system, airfoil_cell_linked_list);

    airfoil_update_cell_linked_list.exec();
    airfoil_update_inner_relation.exec();
    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    cell_linked_list_recording.writeToFile();
    return 0;
}
