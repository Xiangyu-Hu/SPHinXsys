/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/triangle_prism.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0; // Domain width.
Real DH = 2.0;  // Domain height.

//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
//** The number of fluid particles along the cross-seciton should be a int *
Real num_fluid_cross_section = 40.0;

//Real y_p_constant = dp_f/2.0;
Real y_p_constant = 0.05;

Real dp_0 = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0);
Real dp_f = dp_0;
Real dp_d = dp_0;

Real offset_distance = y_p_constant - dp_f / 2.0; //** Basically offset distance is large than or equal to 0 *

Real DL_sponge = dp_f * 20;
Real BW = dp_d * 4.0;

Vecd domain_lower_bound(-DL_sponge - 2.0 * BW, -BW);
Vecd domain_upper_bound(2.0 * DL, 2.0 * DH);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Material properties of the fluid and flow parameters.
//----------------------------------------------------------------------
Real U_inlet = 1.0;
Real U_max = 3.0 * U_inlet;
Real U_f = U_inlet; //*Characteristic velo is regarded as average velo here
Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
Real rho0_f = 1.0;                                            /**< Density. */
Real mu_f = 0.00005;
//----------------------------------------------------------------------
//	define the imported model.
//----------------------------------------------------------------------
std::vector<Vecd> createChannelWallShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
std::vector<Vecd> createComputationalDomainShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge -2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge -2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge -2.0 * BW, 0.0));

    return water_block_shape;
}
//----------------------------------------------------------------------
//	The inflow and outflow regions.
//----------------------------------------------------------------------
Real DH_dummy = DH - 2.0 * offset_distance;
Vecd emitter_halfsize = Vecd(0.5 * BW, 0.5 * DH_dummy);
Vecd emitter_translation = Vecd(-DL_sponge, 0.0) + emitter_halfsize+ Vecd(0.0, offset_distance);
//Vecd emitter_translation_local = Vecd(-DL_sponge, 0.0) + emitter_halfsize + Vecd(0.0, offset_distance);
//Vecd emitter_translation = emitter_translation_local - Vecd(0.0, offset_distance/2.0);

Vecd inlet_buffer_halfsize = Vecd(0.5 * DL_sponge, 0.5 * DH - offset_distance);
Vecd inlet_buffer_translation = Vecd(-DL_sponge, 0.0) + inlet_buffer_halfsize + Vecd(0.0, offset_distance);

Vecd disposer_halfsize = Vecd(0.5 * BW, 0.75 * DH);
Vecd disposer_translation = Vecd(DL, DH + 0.25 * DH) - disposer_halfsize;




class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createChannelWallShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance + BW, outer_boundary, "OuterBoundary");
        
        MultiPolygon computational_domain(createComputationalDomainShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(-offset_distance, computational_domain, "ComputationalDomain");
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody emitter_block(
        sph_system, makeShared<TransformShape<GeometricShapeBox>>(
            Transform(emitter_translation), emitter_halfsize, "EmitterBody"));
    emitter_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    emitter_block.generateParticles<ParticleGeneratorLattice>();
    auto level_set_shape2 = emitter_block.defineBodyLevelSetShape();
    level_set_shape2->correctLevelSetSign();
    level_set_shape2->writeLevelSet(sph_system);


    RealBody imported_model(sph_system, makeShared<WallBoundary>("Wall"));
    // level set shape is used for particle relaxation
    //imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    //imported_model.defineBodyLevelSetShape()->correctLevelSetSign();
    auto level_set_shape = imported_model.defineBodyLevelSetShape();
    level_set_shape->correctLevelSetSign();
    level_set_shape->writeLevelSet(sph_system);


    imported_model.defineParticlesAndMaterial();
    imported_model.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_imported_model_to_vtp({imported_model});
    MeshRecordingToPlt write_cell_linked_list(sph_system, imported_model.getCellLinkedList());
    BodyStatesRecordingToVtp write_emitter_body_to_vtp({ emitter_block });
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation imported_model_inner(imported_model);
    InnerRelation emitter_block_inner(emitter_block);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
    SimpleDynamics<RandomizeParticlePosition> random_emitter_block_particles(emitter_block);
    /** A  Physics relaxation step. */
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner(imported_model_inner);
    RelaxationStepLevelSetCorrectionInner relaxation_step_inner2(emitter_block_inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    //random_emitter_block_particles.exec(0.25);

    relaxation_step_inner.SurfaceBounding().exec();
    //relaxation_step_inner2.SurfaceBounding().exec();

    write_imported_model_to_vtp.writeToFile(0.0);
    write_emitter_body_to_vtp.writeToFile(0.0);

    imported_model.updateCellLinkedList();
    //emitter_block.updateCellLinkedList();
    write_cell_linked_list.writeToFile(0.0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        //relaxation_step_inner2.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
            write_imported_model_to_vtp.writeToFile(ite_p);
            //write_emitter_body_to_vtp.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process of imported model finish !" << std::endl;

    return 0;
}
