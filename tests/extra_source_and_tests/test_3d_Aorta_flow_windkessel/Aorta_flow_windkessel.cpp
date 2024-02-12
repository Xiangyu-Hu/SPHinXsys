/**
 * @file 	pulsatile_Aorta_flow.cpp
 * @brief 	3D pulsatile Aorta flow example
 * @details This is the one of the basic test cases for pressure boundary condition and bidirectional buffer.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h" 
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
#include "bidirectional_buffer.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;

std::string full_path_to_stl_file_wall = "./input/wall.stl";
std::string full_path_to_stl_file_fluid = "./input/blood.stl";

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-6.0E-2, -4.0E-2, -2.0E-2), Vecd(3.0E-2, 10.0E-2, 15.0E-2));
/**
 * @brief Material properties of the fluid.
 */
Real Inlet_pressure = 0.0;
Real Outlet_pressure = 0.0;
Real rho0_f = 1060.0;                   
Real mu_f = 0.00355;
Real U_f = 3.0;
Real c_f = 10.0*U_f;

Real resolution_ref = 0.06E-2;
Real scaling = 1.0E-2;
Vecd translation(0.0, 0.0, 0.0);
StdVec<Vecd> observer_location = {Vecd(1.57E-2, 4.36E-2, 5.7E-2)};

Vecd buffer_halfsize_inlet = Vecd(2.1 * resolution_ref, 3.0E-2, 3.0E-2);
Vecd buffer_translation_inlet = Vecd(-0.9768E-2, 4.6112E-2, 3.0052E-2);
Vecd normal_vector_inlet = Vecd(0.1000, 0.1665, 0.9810);
Vecd vector_inlet = Vecd(0, -0.9859, 0.1674);

Vecd buffer_halfsize_1 = Vecd(1.5 * resolution_ref, 1.0E-2, 1.0E-2);
Vecd buffer_translation_1 = Vecd(-1.2562E-2, 4.4252E-2, 10.0148E-2);
Vecd normal_vector_1 = Vecd(0.6420, 0.4110, 0.6472);
Vecd vector_1 = Vecd(0, -0.8442, 0.5360);

Vecd buffer_halfsize_2 = Vecd(1.6 * resolution_ref, 1.0E-2, 1.0E-2);
Vecd buffer_translation_2 = Vecd(-2.6303E-2, 3.0594E-2, 10.6919E-2);
Vecd normal_vector_2 = Vecd(-0.0988, 0.0485, 0.9939);
Vecd vector_2 = Vecd(0, -0.9988, 0.0488);

Vecd buffer_halfsize_3 = Vecd(1.6 * resolution_ref, 0.6E-2, 0.6E-2);
Vecd buffer_translation_3 = Vecd(-2.8585E-2, 1.8357E-2, 9.8034E-2);
Vecd normal_vector_3 = Vecd(-0.1471, -0.1813, 0.9724);
Vecd vector_3 = Vecd(0, -0.9831, -0.1833);

Vecd buffer_halfsize_4 = Vecd(1.1 * resolution_ref, 1.0E-2, 1.0E-2);
Vecd buffer_translation_4 = Vecd(-1.0946E-2, 1.0386E-2, 9.5016E-2);
Vecd normal_vector_4 = Vecd(0.5675, 0.4280, 0.7034);
Vecd vector_4 = Vecd(0, -0.8543, 0.5197);
  
Vecd buffer_halfsize_5 = Vecd(1.3 * resolution_ref, 2.0E-2, 2.0E-2);
Vecd buffer_translation_5 = Vecd(-1.6791E-2, -0.8069E-2, 0.5017E-2);
Vecd normal_vector_5 = Vecd(0.0327, -0.0729, 0.9968);
Vecd vector_5 = Vecd(0, -0.9973, -0.0729);
/**
 * @brief 	water body definition.
 */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_fluid, translation, scaling);
        subtract<TriangleMeshShapeSTL>(full_path_to_stl_file_wall, translation, scaling);
    }
};
/**
 * @brief 	wall body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_wall, translation, scaling);
    }
};
/**
 * @brief 	Bidirectional buffer definition.
 */
class BidirectionalBufferConditionInlet : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferConditionInlet(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                      size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Inlet_pressure;
        return pressure;
    }
};

class BidirectionalBufferCondition1 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferCondition1(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class BidirectionalBufferCondition2 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferCondition2(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class BidirectionalBufferCondition3 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferCondition3(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class BidirectionalBufferCondition4 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferCondition4(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class BidirectionalBufferCondition5 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferCondition5(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
/**
 * @brief 	inflow velocity definition.
 */
struct InflowVelocity
{
    Real u_ave;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ave(0.0){}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        u_ave = 0.3782;
        Real a[8] = {-0.1812,0.1276,-0.08981,0.04347,-0.05412,0.02642,0.008946,-0.009005};
        Real b[8] = {-0.07725,0.01466,0.004295,-0.06679,0.05679,-0.01878,0.01869,-0.01888};
        for (size_t i = 0; i < 8; i++)
            u_ave += a[i] * cos(8.302 * (i + 1) * run_time) + b[i] * sin(8.302 * (i + 1) * run_time);

        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;

        return target_velocity;
    }
};
/**
 * @brief 	outlet pressure boundary definition.
 */
class DisposerOutflowDeletionWithWindkessel_1 : public fluid_dynamics::DisposerOutflowDeletion
{
  public:
    DisposerOutflowDeletionWithWindkessel_1(BodyAlignedBoxByCell &aligned_box_part, int axis) : 
        DisposerOutflowDeletion(aligned_box_part, axis), flow_rate_(0.0), Vol_(particles_->Vol_){};
    virtual ~DisposerOutflowDeletionWithWindkessel_1(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        mutex_switch_to_buffer_.lock();
        while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
        {
            particles_->switchToBufferParticle(index_i);
            flow_rate_ += Vol_[index_i];
        }
        mutex_switch_to_buffer_.unlock();
    };
    Real flow_rate_;

  protected:
    StdLargeVec<Real> &Vol_;
};

class OutflowPressure_1 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_1(DisposerOutflowDeletionWithWindkessel_1 &outlet_windkessel, BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), flow_rate(outlet_windkessel.flow_rate_),
          Q_0(0.0), Q_n(0.0), p_0(0.0), p_n(0.0), R1(1.18E8), R2(1.84E9), C(7.7E-10), delta_t(0.006), current_flow_rate(0.0), 
        previous_flow_rate(0.0), count(1) {}
    void getFlowRate()
    {
        Real run_time = GlobalStaticVariables::physical_time_;

        if (int(run_time / delta_t) > count)
        {
            Q_0 = Q_n;
            p_0 = p_n;
            current_flow_rate = flow_rate - previous_flow_rate;
            previous_flow_rate = flow_rate;
            count += 1;
        }
    };

    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Q_n = current_flow_rate / delta_t - 0.0002*0.049;  
        p_n = ((Q_n * (1.0 + R1 / R2) + C * R1 * (Q_n - Q_0) / delta_t) * delta_t / C + p_0) / (1.0 + delta_t / (C * R2));

        return p_n;
    }
    void setupDynamics(Real dt = 0.0) override {}

  protected:
    Real &flow_rate, Q_0, Q_n, p_0, p_n, R1, R2, C, delta_t, current_flow_rate, previous_flow_rate;
    int count;
};

class DisposerOutflowDeletionWithWindkessel_2 : public fluid_dynamics::DisposerOutflowDeletion
{
  public:
    DisposerOutflowDeletionWithWindkessel_2(BodyAlignedBoxByCell &aligned_box_part, int axis) : 
        DisposerOutflowDeletion(aligned_box_part, axis), flow_rate_(0.0), Vol_(particles_->Vol_){};
    virtual ~DisposerOutflowDeletionWithWindkessel_2(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        mutex_switch_to_buffer_.lock();
        while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
        {
            particles_->switchToBufferParticle(index_i);
            flow_rate_ += Vol_[index_i];
        }
        mutex_switch_to_buffer_.unlock();
    };
    Real flow_rate_;

  protected:
    StdLargeVec<Real> &Vol_;
};

class OutflowPressure_2 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_2(DisposerOutflowDeletionWithWindkessel_2 &outlet_windkessel, BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), flow_rate(outlet_windkessel.flow_rate_),
          Q_0(0.0), Q_n(0.0), p_0(0.0), p_n(0.0), R1(1.04E8), R2(1.63E9), C(8.74E-10), delta_t(0.006), 
        current_flow_rate(0.0), previous_flow_rate(0.0), count(1) {}
    void getFlowRate()
    {
        Real run_time = GlobalStaticVariables::physical_time_;

        if (int(run_time / delta_t) > count)
        {
            Q_0 = Q_n;
            p_0 = p_n;
            current_flow_rate = flow_rate - previous_flow_rate;
            previous_flow_rate = flow_rate;
            count += 1;
        }
    };

    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Q_n = current_flow_rate / delta_t - 0.0002*0.05;  
        p_n = ((Q_n * (1.0 + R1 / R2) + C * R1 * (Q_n - Q_0) / delta_t) * delta_t / C + p_0) / (1.0 + delta_t / (C * R2));

        return p_n;
    }
    void setupDynamics(Real dt = 0.0) override {}

  protected:
    Real &flow_rate, Q_0, Q_n, p_0, p_n, R1, R2, C, delta_t, current_flow_rate, previous_flow_rate;
    int count;
};

class DisposerOutflowDeletionWithWindkessel_3 : public fluid_dynamics::DisposerOutflowDeletion
{
  public:
    DisposerOutflowDeletionWithWindkessel_3(BodyAlignedBoxByCell &aligned_box_part, int axis) : 
        DisposerOutflowDeletion(aligned_box_part, axis), flow_rate_(0.0), Vol_(particles_->Vol_){};
    virtual ~DisposerOutflowDeletionWithWindkessel_3(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        mutex_switch_to_buffer_.lock();
        while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
        {
            particles_->switchToBufferParticle(index_i);
            flow_rate_ += Vol_[index_i];
        }
        mutex_switch_to_buffer_.unlock();
    };
    Real flow_rate_;

  protected:
    StdLargeVec<Real> &Vol_;
};

class OutflowPressure_3 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_3(DisposerOutflowDeletionWithWindkessel_3 &outlet_windkessel, BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), flow_rate(outlet_windkessel.flow_rate_),
          Q_0(0.0), Q_n(0.0), p_0(0.0), p_n(0.0), R1(1.18E8), R2(1.84E9), C(7.7E-10), delta_t(0.006), 
        current_flow_rate(0.0), previous_flow_rate(0.0), count(1) {}
    void getFlowRate()
    {
        Real run_time = GlobalStaticVariables::physical_time_;

        if (int(run_time / delta_t) > count)
        {
            Q_0 = Q_n;
            p_0 = p_n;
            current_flow_rate = flow_rate - previous_flow_rate;
            previous_flow_rate = flow_rate;
            count += 1;
        }
    };

    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Q_n = current_flow_rate / delta_t - 0.0002*0.034;  
        p_n = ((Q_n * (1.0 + R1 / R2) + C * R1 * (Q_n - Q_0) / delta_t) * delta_t / C + p_0) / (1.0 + delta_t / (C * R2));

        return p_n;
    }
    void setupDynamics(Real dt = 0.0) override {}

  protected:
    Real &flow_rate, Q_0, Q_n, p_0, p_n, R1, R2, C, delta_t, current_flow_rate, previous_flow_rate;
    int count;
};

class DisposerOutflowDeletionWithWindkessel_4 : public fluid_dynamics::DisposerOutflowDeletion
{
  public:
    DisposerOutflowDeletionWithWindkessel_4(BodyAlignedBoxByCell &aligned_box_part, int axis) : 
        DisposerOutflowDeletion(aligned_box_part, axis), flow_rate_(0.0), Vol_(particles_->Vol_){};
    virtual ~DisposerOutflowDeletionWithWindkessel_4(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        mutex_switch_to_buffer_.lock();
        while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
        {
            particles_->switchToBufferParticle(index_i);
            flow_rate_ += Vol_[index_i];
        }
        mutex_switch_to_buffer_.unlock();
    };
    Real flow_rate_;

  protected:
    StdLargeVec<Real> &Vol_;
};

class OutflowPressure_4 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_4(DisposerOutflowDeletionWithWindkessel_4 &outlet_windkessel, BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), flow_rate(outlet_windkessel.flow_rate_),
          Q_0(0.0), Q_n(0.0), p_0(0.0), p_n(0.0), R1(9.7E7), R2(1.52E9), C(9.34E-10), delta_t(0.006), 
        current_flow_rate(0.0), previous_flow_rate(0.0), count(1) {}
    void getFlowRate()
    {
        Real run_time = GlobalStaticVariables::physical_time_;

        if (int(run_time / delta_t) > count)
        {
            Q_0 = Q_n;
            p_0 = p_n;
            current_flow_rate = flow_rate - previous_flow_rate;
            previous_flow_rate = flow_rate;
            count += 1;
        }
    };

    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Q_n = current_flow_rate / delta_t - 0.0002*0.059;  
        p_n = ((Q_n * (1.0 + R1 / R2) + C * R1 * (Q_n - Q_0) / delta_t) * delta_t / C + p_0) / (1.0 + delta_t / (C * R2));

        return p_n;
    }
    void setupDynamics(Real dt = 0.0) override {}

  protected:
    Real &flow_rate, Q_0, Q_n, p_0, p_n, R1, R2, C, delta_t, current_flow_rate, previous_flow_rate;
    int count;
};

class DisposerOutflowDeletionWithWindkessel_5 : public fluid_dynamics::DisposerOutflowDeletion
{
  public:
    DisposerOutflowDeletionWithWindkessel_5(BodyAlignedBoxByCell &aligned_box_part, int axis) : 
        DisposerOutflowDeletion(aligned_box_part, axis), flow_rate_(0.0), Vol_(particles_->Vol_){};
    virtual ~DisposerOutflowDeletionWithWindkessel_5(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        mutex_switch_to_buffer_.lock();
        while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
        {
            particles_->switchToBufferParticle(index_i);
            flow_rate_ += Vol_[index_i];
        }
        mutex_switch_to_buffer_.unlock();
    };
    
    Real flow_rate_;

  protected:
    StdLargeVec<Real> &Vol_;  
};

class OutflowPressure_5 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_5(DisposerOutflowDeletionWithWindkessel_5 &outlet_windkessel, BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), flow_rate(outlet_windkessel.flow_rate_),
          Q_0(0.0), Q_n(0.0), p_0(0.0), p_n(0.0), R1(1.88E7), R2(2.95E8), C(4.82E-9), delta_t(0.006), current_flow_rate(0.0), 
          previous_flow_rate(0.0), count(1) {}

    void getFlowRate()
    {
        Real run_time = GlobalStaticVariables::physical_time_;

        if (int(run_time / delta_t) > count)
        {
            Q_0 = Q_n;
            p_0 = p_n;
            current_flow_rate = flow_rate - previous_flow_rate;
            previous_flow_rate = flow_rate;
            count += 1;
        }      
    };

    Real getTargetPressure(Real dt) override
    {                 
        Real run_time = GlobalStaticVariables::physical_time_;
        Q_n = current_flow_rate / delta_t - 0.0002*0.48;                                           
        p_n = ((Q_n * (1.0 + R1 / R2) + C * R1 * (Q_n - Q_0) / delta_t) * delta_t / C + p_0) / (1.0 + delta_t / (C * R2));

        return p_n;
    }
    void setupDynamics(Real dt = 0.0) override {}

  protected:
    Real &flow_rate, Q_0, Q_n, p_0, p_n, R1, R2, C, delta_t, current_flow_rate, previous_flow_rate;
    int count;
};

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setGenerateRegressionData(false);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->correctLevelSetSign();
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<ParticleGeneratorReload>(water_block.getName())
        : water_block.generateParticles<ParticleGeneratorLattice>();

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ParticleGeneratorObserver>(observer_location);
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();
    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation water_block_inner(water_block);
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> water_block_random_inserted_body_particles(water_block);
        SimpleDynamics<RandomizeParticlePosition> wall_boundary_random_inserted_body_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp water_block_write_inserted_body_to_vtp({&water_block});
        BodyStatesRecordingToVtp wall_boundary_write_inserted_body_to_vtp({&wall_boundary});
        /** Write the particle reload files. */
        ReloadParticleIO water_block_write_particle_reload_files({&water_block});
        ReloadParticleIO wall_boundary_write_particle_reload_files({&wall_boundary});
        /** A  Physics relaxation step. */
        RelaxationStepInner water_block_relaxation_step_inner(water_block_inner);
        RelaxationStepInner wall_boundary_relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
       //----------------------------------------------------------------------
        water_block_random_inserted_body_particles.exec(0.25);
        wall_boundary_random_inserted_body_particles.exec(0.25);
        water_block_relaxation_step_inner.SurfaceBounding().exec();
        wall_boundary_relaxation_step_inner.SurfaceBounding().exec();
        water_block_write_inserted_body_to_vtp.writeToFile(0);
        wall_boundary_write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 3000)
        {
            water_block_relaxation_step_inner.exec();
            wall_boundary_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                water_block_write_inserted_body_to_vtp.writeToFile(ite_p);
                wall_boundary_write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        water_block_write_particle_reload_files.writeToFile(0);
        wall_boundary_write_particle_reload_files.writeToFile(0);
        return 0;
    }
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    /**
     * @brief 	Methods used for time stepping.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** delete outflow particles */
    BodyAlignedBoxByCell disposer_inlet(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_inlet)),
                                                                                                      -vector_inlet),
                                                                                           Vecd(buffer_translation_inlet)),
                                                                                 buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_inlet(disposer_inlet, xAxis);
    BodyAlignedBoxByCell disposer_1(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_1)),
                                                                                                  vector_1),
                                                                                       Vecd(buffer_translation_1)),
                                                                             buffer_halfsize_1));
    SimpleDynamics<DisposerOutflowDeletionWithWindkessel_1> disposer_outflow_deletion_1(disposer_1, xAxis);
    BodyAlignedBoxByCell disposer_2(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_2)),
                                                                                                  vector_2),
                                                                                       Vecd(buffer_translation_2)),
                                                                             buffer_halfsize_2));
    SimpleDynamics<DisposerOutflowDeletionWithWindkessel_2> disposer_outflow_deletion_2(disposer_2, xAxis);
    BodyAlignedBoxByCell disposer_3(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_3)),
                                                                                                  vector_3),
                                                                                       Vecd(buffer_translation_3)),
                                                                             buffer_halfsize_3));
    SimpleDynamics<DisposerOutflowDeletionWithWindkessel_3> disposer_outflow_deletion_3(disposer_3, xAxis);
    BodyAlignedBoxByCell disposer_4(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_4)),
                                                                                                  vector_4),
                                                                                       Vecd(buffer_translation_4)),
                                                                             buffer_halfsize_4));
    SimpleDynamics<DisposerOutflowDeletionWithWindkessel_4> disposer_outflow_deletion_4(disposer_4, xAxis);
    BodyAlignedBoxByCell disposer_5(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_5)),
                                                                                                  -vector_5),
                                                                                       Vecd(buffer_translation_5)),
                                                                             buffer_halfsize_5));
    SimpleDynamics<DisposerOutflowDeletionWithWindkessel_5> disposer_outflow_deletion_5(disposer_5, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        free_stream_surface_indicator(water_block_inner, water_block_contact);
    /** bidrectional buffer */
    BidirectionalBufferConditionInlet emitter_inflow_injection_inlet(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)), vector_inlet), (buffer_translation_inlet)), buffer_halfsize_inlet), 100, xAxis);
    BidirectionalBufferCondition1 emitter_inflow_injection1(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_1)), -vector_1), (buffer_translation_1)), buffer_halfsize_1), 100, xAxis);
    BidirectionalBufferCondition2 emitter_inflow_injection2(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_2)), -vector_2), (buffer_translation_2)), buffer_halfsize_2), 100, xAxis);
    BidirectionalBufferCondition3 emitter_inflow_injection3(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_3)), -vector_3), (buffer_translation_3)), buffer_halfsize_3), 100, xAxis);
    BidirectionalBufferCondition4 emitter_inflow_injection4(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_4)), -vector_4), (buffer_translation_4)), buffer_halfsize_4), 100, xAxis);
    BidirectionalBufferCondition5 emitter_inflow_injection5(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_5)), vector_5), (buffer_translation_5)), buffer_halfsize_5), 100, xAxis);
    /** output parameters */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("PositionDivergence");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /* Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** pressure boundary condition. */
    BodyAlignedBoxByCell inflow_pressure_region(water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)), vector_inlet),
                                                                                                   (buffer_translation_inlet)),
                                                                                         buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(inflow_pressure_region);
    BodyRegionByCell outflow_pressure_region1(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_1)), vector_1),
                                                                                                                   (buffer_translation_1)),
                                                                                                         buffer_halfsize_1));
    SimpleDynamics<OutflowPressure_1> outflow_pressure_condition1(disposer_outflow_deletion_1, outflow_pressure_region1, normal_vector_1);
    BodyRegionByCell outflow_pressure_region2(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_2)), vector_2),
                                                                                                                   (buffer_translation_2)),
                                                                                                         buffer_halfsize_2));
    SimpleDynamics<OutflowPressure_2> outflow_pressure_condition2(disposer_outflow_deletion_2, outflow_pressure_region2, normal_vector_2);
    BodyRegionByCell outflow_pressure_region3(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_3)), vector_3),
                                                                                                                   (buffer_translation_3)),
                                                                                                         buffer_halfsize_3));
    SimpleDynamics<OutflowPressure_3> outflow_pressure_condition3(disposer_outflow_deletion_3, outflow_pressure_region3, normal_vector_3);
    BodyRegionByCell outflow_pressure_region4(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_4)), vector_4),
                                                                                                                   (buffer_translation_4)),
                                                                                                         buffer_halfsize_4));
    SimpleDynamics<OutflowPressure_4> outflow_pressure_condition4(disposer_outflow_deletion_4, outflow_pressure_region4, normal_vector_4);
    BodyRegionByCell outflow_pressure_region5(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_5)), vector_5),
                                                                                                                   (buffer_translation_5)),
                                                                                                         buffer_halfsize_5));
    SimpleDynamics<OutflowPressure_5> outflow_pressure_condition5(disposer_outflow_deletion_5, outflow_pressure_region5, normal_vector_5);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact);

    /**
     * @brief Output.
     */
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_centerline_velocity("Velocity", velocity_observer_contact);
    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    free_stream_surface_indicator.exec();
    emitter_inflow_injection_inlet.tag_buffer_particles.exec();
    emitter_inflow_injection1.tag_buffer_particles.exec();
    emitter_inflow_injection2.tag_buffer_particles.exec();
    emitter_inflow_injection3.tag_buffer_particles.exec();
    emitter_inflow_injection4.tag_buffer_particles.exec();
    emitter_inflow_injection5.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();

    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = 0.0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 1.0;      /**< End time. */
    Real Output_Time = 0.005; /**< Time stamps for output of body states. */
    Real dt = 0.0;            /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    /** Output the start states of bodies. */
    body_states_recording.writeToFile(0);
    write_centerline_velocity.writeToFile(number_of_iterations);
    /**
    * @brief 	Main loop starts here.
    */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                emitter_buffer_inflow_condition.exec();
                outflow_pressure_condition1.getFlowRate();
                outflow_pressure_condition1.exec(dt);
                outflow_pressure_condition2.getFlowRate();
                outflow_pressure_condition2.exec(dt);
                outflow_pressure_condition3.getFlowRate();
                outflow_pressure_condition3.exec(dt);
                outflow_pressure_condition4.getFlowRate();
                outflow_pressure_condition4.exec(dt);
                outflow_pressure_condition5.getFlowRate();
                outflow_pressure_condition5.exec(dt);
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
            emitter_inflow_injection_inlet.injection.exec();
            emitter_inflow_injection1.injection.exec();
            emitter_inflow_injection2.injection.exec();
            emitter_inflow_injection3.injection.exec();
            emitter_inflow_injection4.injection.exec();
            emitter_inflow_injection5.injection.exec();
            disposer_outflow_deletion_inlet.exec();
            disposer_outflow_deletion_1.exec();
            disposer_outflow_deletion_2.exec();
            disposer_outflow_deletion_3.exec();
            disposer_outflow_deletion_4.exec();
            disposer_outflow_deletion_5.exec();

            water_block.updateCellLinkedList();
            water_block_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            free_stream_surface_indicator.exec();
            emitter_inflow_injection_inlet.tag_buffer_particles.exec();
            emitter_inflow_injection1.tag_buffer_particles.exec();
            emitter_inflow_injection2.tag_buffer_particles.exec();
            emitter_inflow_injection3.tag_buffer_particles.exec();
            emitter_inflow_injection4.tag_buffer_particles.exec();
            emitter_inflow_injection5.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        velocity_observer_contact.updateConfiguration();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

      if (sph_system.GenerateRegressionData())
       {
          write_centerline_velocity.generateDataBase(1.0e-3);
       }
       else
       {
          write_centerline_velocity.testResult();
       }

    return 0;

}
