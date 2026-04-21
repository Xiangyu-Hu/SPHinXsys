/** @file 	test_body_part_boolean.cpp
 *  @brief 	Tests for BodyPartByParticle with custom variable tagging criteria.
 *          The geometry is copied from test_2d_dambreak.
 *  @author Hong Zhu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
//----------------------------------------------------------------------
//	Customized criteria
//----------------------------------------------------------------------
class UnionCriteria
{
  public:
    UnionCriteria() : x_plus_2y_(nullptr), x_multiply_y_(nullptr) {}

    void setupBaseParticles(BaseParticles &base_particles)
    {
        x_plus_2y_ = base_particles.getVariableDataByName<Real>("XPlus2Y");
        x_multiply_y_ = base_particles.getVariableDataByName<Real>("XMultiplyY");
    }

    bool operator()(size_t index_i) const
    {
        return (0.5 <= x_plus_2y_[index_i] && x_plus_2y_[index_i] <= 1.5) ||
               (0.1 <= x_multiply_y_[index_i] && x_multiply_y_[index_i] <= 0.4);
    }

  private:
    Real *x_plus_2y_;
    Real *x_multiply_y_;
};

class IntersectionCriteria
{
  public:
    IntersectionCriteria() : x_plus_2y_(nullptr), x_multiply_y_(nullptr) {}

    void setupBaseParticles(BaseParticles &base_particles)
    {
        x_plus_2y_ = base_particles.getVariableDataByName<Real>("XPlus2Y");
        x_multiply_y_ = base_particles.getVariableDataByName<Real>("XMultiplyY");
    }

    bool operator()(size_t particle_index) const
    {
        return (0.5 <= x_plus_2y_[particle_index] && x_plus_2y_[particle_index] <= 1.5) &&
               (0.1 <= x_multiply_y_[particle_index] && x_multiply_y_[particle_index] <= 0.4);
    }

  private:
    Real *x_plus_2y_;
    Real *x_multiply_y_;
};

class DifferenceCriteria
{
  public:
    DifferenceCriteria() : x_plus_2y_(nullptr), x_multiply_y_(nullptr) {}

    void setupBaseParticles(BaseParticles &base_particles)
    {
        x_plus_2y_ = base_particles.getVariableDataByName<Real>("XPlus2Y");
        x_multiply_y_ = base_particles.getVariableDataByName<Real>("XMultiplyY");
    }

    bool operator()(size_t particle_index) const
    {
        return (0.5 <= x_plus_2y_[particle_index] && x_plus_2y_[particle_index] <= 1.5) &&
               !(0.1 <= x_multiply_y_[particle_index] && x_multiply_y_[particle_index] <= 0.4);
    }

  private:
    Real *x_plus_2y_;
    Real *x_multiply_y_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    // water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    BaseParticles &base_particles = water_block.getBaseParticles();
    Real *x_plus_2y = base_particles.registerStateVariableData<Real>("XPlus2Y", Real(0));
    Real *x_multiply = base_particles.registerStateVariableData<Real>("XMultiplyY", Real(0));
    Vecd *pos = base_particles.getVariableDataByName<Vecd>("Position");
    const UnsignedInt total_real_particles = base_particles.TotalRealParticles();
    for (UnsignedInt index_i = 0; index_i < total_real_particles; ++index_i)
    {
        x_plus_2y[index_i] = pos[index_i][0] + 2.0 * pos[index_i][1];
        x_multiply[index_i] = pos[index_i][0] * pos[index_i][1];
    }
    BodyPartByRealVar part_x_plus_2y(water_block, "XPlus2Y", 0.5, 1.5);
    BodyPartByRealVar part_x_multiply(water_block, "XMultiplyY", 0.1, 0.4);

    // Demonstrate custom TagCriteria through the generalized BodyPartByParticle constructor.
    BodyPartByParticle part_union(water_block, UnionCriteria());
    BodyPartByParticle part_intersection(water_block, IntersectionCriteria());
    BodyPartByParticle part_difference(water_block, DifferenceCriteria());
    const size_t count_x_plus_2y = part_x_plus_2y.SizeOfLoopRange();
    const size_t count_x_multiply = part_x_multiply.SizeOfLoopRange();
    const size_t count_union = part_union.SizeOfLoopRange();
    const size_t count_intersection = part_intersection.SizeOfLoopRange();
    const size_t count_difference = part_difference.SizeOfLoopRange();

    const size_t expected_x_plus_2y = 800;
    const size_t expected_x_multiply = 1030;
    const size_t expected_union = 1462;
    const size_t expected_intersection = 368;
    const size_t expected_difference = 432;

    std::cout << "Total particles in PartXPlus2Y: " << count_x_plus_2y << "\n"
              << "Total particles in PartXMultiplyY: " << count_x_multiply << "\n"
              << "Total particles in Union: " << count_union << "\n"
              << "Total particles in Intersection: " << count_intersection << "\n"
              << "Total particles in Difference: " << count_difference << "\n";

    const bool pass = (count_x_plus_2y == expected_x_plus_2y) &&
                      (count_x_multiply == expected_x_multiply) &&
                      (count_union == expected_union) &&
                      (count_intersection == expected_intersection) &&
                      (count_difference == expected_difference);

    if (!pass)
    {
        std::cerr << "Particle count check failed." << std::endl;
        return 1;
    }
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.writeToFile();
    return 0;
}
