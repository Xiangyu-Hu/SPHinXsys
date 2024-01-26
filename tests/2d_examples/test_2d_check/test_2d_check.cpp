/**
 * @file 	test_2d_fluid_around_balloon_shell.cpp
 * @brief 	Test on fluid-shell interaction when 2 shell particles are close to each other
 * @details This is a case to test fluid-shell interaction.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    Kernel *kernel_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Real> W_ijV_j_ttl_contact;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_),
          kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
          inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
        {
            contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl_contact, "TotalKernelContact");
        inner_relation.base_particles_.registerVariable(dW_ijV_je_ij_ttl, "TotalKernelGrad");
        inner_relation.base_particles_.registerVariable(number_of_inner_neighbor, "InnerNeighborNumber");
        inner_relation.base_particles_.registerVariable(number_of_contact_neighbor, "ContactNeighborNumber");
    }

    inline void exec()
    {
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                int N_inner_number = 0;
                int N_contact_number = 0;
                Real W_ijV_j_ttl_i = particles_->Vol_[index_i] * kernel_->W(0, ZeroVecd);
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                double W_ijV_j_ttl_contact_i = 0;
                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_contact_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i + W_ijV_j_ttl_contact_i;
                W_ijV_j_ttl_contact[index_i] = W_ijV_j_ttl_contact_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }

    inline Real get_W_ijV_j(size_t index_i) { return W_ijV_j_ttl[index_i]; }
    inline Vecd get_dW_ijV_je_ij(size_t index_i) { return dW_ijV_je_ij_ttl[index_i]; }
    inline int get_contact_number(size_t index_i) { return number_of_contact_neighbor[index_i]; }
};

/// The given variables must outlive the lifetime of this object !
class ParticleGeneratorDirect final : public SPH::ParticleGenerator
{
    // Particle generation fully performed by SPHBody::generateParticles(...) in base_body.h
    // Enable construction only by a SPHBody to prevent misuse
    friend class SPH::SPHBody;

    const SPH::StdLargeVec<SPH::Vecd> &positions;
    const SPH::StdLargeVec<SPH::Real> &volumes;

    ParticleGeneratorDirect(
        SPH::SPHBody &body,
        const SPH::StdLargeVec<SPH::Vecd> &positions,
        const SPH::StdLargeVec<SPH::Real> &volumes) : ParticleGenerator(body),
                                                      positions(positions),
                                                      volumes(volumes)
    {
    }

    void initializeGeometricVariables() override
    {
        for (size_t i = 0; i < positions.size(); ++i)
        {
            initializePositionAndVolumetricMeasure(positions[i], volumes[i]);
        }
    };
};

class ShellParticleGeneratorDirect : public SurfaceParticleGenerator
{
    const SPH::StdLargeVec<SPH::Vecd> &positions_;
    const SPH::StdLargeVec<SPH::Vecd> &normals_;
    Real resolution_;
    Real thickness_shell_;

  public:
    ShellParticleGeneratorDirect(SPHBody &sph_body,
                                 const SPH::StdLargeVec<SPH::Vecd> &positions,
                                 const SPH::StdLargeVec<SPH::Vecd> &normals, Real resolution, Real thickness)
        : SurfaceParticleGenerator(sph_body),
          positions_(positions), normals_(normals),
          resolution_(resolution), thickness_shell_(thickness){};
    void initializeGeometricVariables() override
    {
        for (size_t i = 0; i < positions_.size(); ++i)
        {
            initializePositionAndVolumetricMeasure(positions_[i], resolution_);
            initializeSurfaceProperties(normals_[i], thickness_shell_);
        }
    }
};

/** create a wall outer shape */
std::vector<Vecd> createPlateSolidShape(Real half_length, Real thickness)
{
    // geometry
    std::vector<Vecd> outer_shape;
    outer_shape.push_back(Vecd(-half_length, -thickness));
    outer_shape.push_back(Vecd(-half_length, 0));
    outer_shape.push_back(Vecd(half_length, 0));
    outer_shape.push_back(Vecd(half_length, -thickness));
    outer_shape.push_back(Vecd(-half_length, -thickness));

    return outer_shape;
}

class PlateSolidShape : public MultiPolygonShape
{
  public:
    PlateSolidShape(const std::string &shape_name, Real half_length, Real thickness) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createPlateSolidShape(half_length, thickness), ShapeBooleanOps::add);
    }
};

class PlateShellParticleGeneratorDirect : public SurfaceParticleGenerator
{
    Real half_length_;
    Real resolution_;
    Real thickness_shell_;

  public:
    PlateShellParticleGeneratorDirect(SPHBody &sph_body, Real half_length, Real resolution, Real thickness)
        : SurfaceParticleGenerator(sph_body),
          half_length_(half_length), resolution_(resolution), thickness_shell_(thickness){};
    void initializeGeometricVariables() override
    {
        Real x = -half_length_ + 0.5 * resolution_;
        Real y = -0.5 * resolution_;
        while (x < half_length_)
        {
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_);
            initializeSurfaceProperties(Vecd(0, 1), thickness_shell_);
            x += resolution_;
        }
    }
};

class CircleSolidShape : public MultiPolygonShape
{
  public:
    CircleSolidShape(const std::string &shape_name, Real outer_radius, Real thickness, const Vecd &center) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(center, outer_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(center, outer_radius - thickness, 100, ShapeBooleanOps::sub);
    }
};

class CircleShellParticleGeneratorDirect : public SurfaceParticleGenerator
{
    Real mid_radius_;
    Real resolution_;
    Real thickness_shell_;
    Vecd center_;

  public:
    CircleShellParticleGeneratorDirect(SPHBody &sph_body, Real mid_radius, Real resolution, Real thickness, const Vecd &center)
        : SurfaceParticleGenerator(sph_body),
          mid_radius_(mid_radius), resolution_(resolution), thickness_shell_(thickness),
          center_(center){};
    void initializeGeometricVariables() override
    {
        Real dtheta = resolution_ / mid_radius_;
        Real theta = 0.5 * dtheta;
        while (theta < Pi - 0.5 * dtheta)
        {
            Real x = mid_radius_ * sin(theta);
            Real y = mid_radius_ * cos(theta);
            initializePositionAndVolumetricMeasure(Vecd(x, y) + center_, resolution_);
            initializeSurfaceProperties(Vecd(x, y).normalized(), thickness_shell_);
            initializePositionAndVolumetricMeasure(Vecd(-x, y) + center_, resolution_);
            initializeSurfaceProperties(Vecd(-x, y).normalized(), thickness_shell_);
            theta += dtheta;
        }
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
// int main(int ac, char *av[])
// {
//     const Real scale = 0.001;
//     const Real resolution_ref = 0.5 * scale;
//     const Real resolution_shell = 0.5 * resolution_ref;
//     const Real shell_thickness = 6 * resolution_shell;
//     //----------------------------------------------------------------------
//     //	Global parameters on the fluid properties
//     //----------------------------------------------------------------------
//     const Real rho0_f = 1056.0;                           /**< Reference density of fluid. */
//     const Real mu_f = 3.5e-3;                             /**< Dynamics viscosity. */
//     const Real Re = 100.0;                                /**< Reynolds number. */
//     const Real U_f = Re * mu_f / rho0_f / (12.0 * scale); /**< Characteristic velocity. */
//     const Real U_max = 1.5 * U_f;
//     /** Reference sound speed needs to consider the flow speed in the narrow channels. */
//     const Real c_f = 10.0 * U_max;
//     //----------------------------------------------------------------------
//     //	Global parameters on the solid properties
//     //----------------------------------------------------------------------
//     const Real rho0_s = 1250.0; /**< Reference density.*/
//     const Real hardness = 50;   // Durometer hardnes: 50A
//     const Real youngs_modulus =
//         std::pow(10, 0.0235 * hardness - 0.6403) * 1e4; // actual: 1e6, ref: https://doi.org/10.5254/1.3547752 eq. 12A
//     const Real poisson_ratio = 0.495;

//     BoundingBox system_domain_bounds(Vec2d(-0.5 * resolution_ref, -shell_thickness), Vec2d(0.5 * resolution_ref, resolution_ref));
//     //----------------------------------------------------------------------
//     //	Build up the environment of a SPHSystem with global controls.
//     //----------------------------------------------------------------------
//     SPHSystem sph_system(system_domain_bounds, resolution_ref);
//     IOEnvironment io_environment(sph_system);
//     //----------------------------------------------------------------------
//     //	Creating body, materials and particles.
//     //----------------------------------------------------------------------
//     StdLargeVec<Vecd> water_positions({Vec2d(0, 0.5 * resolution_ref),
//                                        Vec2d(0, 1.5 * resolution_ref)});
//     StdLargeVec<Real> water_volumes({resolution_ref * resolution_ref, resolution_ref * resolution_ref});
//     FluidBody water_1(sph_system, makeShared<DefaultShape>("WaterBody1"));
//     water_1.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
//     water_1.generateParticles<ParticleGeneratorDirect>(water_positions, water_volumes);

//     FluidBody water_2(sph_system, makeShared<DefaultShape>("WaterBody2"));
//     water_2.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
//     water_2.generateParticles<ParticleGeneratorDirect>(water_positions, water_volumes);

//     StdLargeVec<Vecd> solid_positions;
//     StdLargeVec<Real> solid_volumes;
//     {
//         for (int i = 0; i < int(shell_thickness / resolution_shell); i++)
//         {
//             solid_positions.emplace_back(-0.5 * resolution_shell, -(0.5 + i) * resolution_shell);
//             solid_volumes.emplace_back(resolution_shell * resolution_shell);
//             solid_positions.emplace_back(0.5 * resolution_shell, -(0.5 + i) * resolution_shell);
//             solid_volumes.emplace_back(resolution_shell * resolution_shell);
//         }
//     }
//     SolidBody solid(sph_system, makeShared<DefaultShape>("Solid"));
//     solid.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
//     solid.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
//     solid.generateParticles<ParticleGeneratorDirect>(solid_positions, solid_volumes);
//     for (size_t i = 0; i < solid_positions.size(); i++)
//     {
//         auto &n = *solid.getBaseParticles().getVariableByName<Vecd>("NormalDirection");
//         auto &n0 = *solid.getBaseParticles().getVariableByName<Vecd>("InitialNormalDirection");
//         n[i] = Vec2d(0, 1);
//         n0[i] = Vec2d(0, 1);
//     }

//     StdLargeVec<Vecd> shell_positions({Vecd(-0.5 * resolution_shell, -0.5 * resolution_shell),
//                                        Vecd(0.5 * resolution_shell, -0.5 * resolution_shell)});
//     StdLargeVec<Vecd> shell_normals({Vecd(0, 1.0), Vecd(0, 1.0)});
//     SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
//     shell.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
//     shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
//     shell.generateParticles<ShellParticleGeneratorDirect>(shell_positions, shell_normals, resolution_shell, shell_thickness);
//     //----------------------------------------------------------------------
//     //	Define body relation map.
//     //	The contact map gives the topological connections between the bodies.
//     //	Basically the the range of bodies to build neighbor particle lists.
//     //   Generally, we first define all the inner relations, then the contact relations.
//     //   At last, we define the complex relaxations by combining previous defined
//     //   inner and contact relations.
//     //----------------------------------------------------------------------
//     InnerRelation water_inner_1(water_1);
//     InnerRelation water_inner_2(water_2);
//     InnerRelation solid_inner(solid);
//     InnerRelation shell_inner(shell);
//     ContactRelation water_solid_contact(water_1, {&solid});
//     ContactRelationToShell water_shell_contact(water_2, {&shell});
//     ContactRelation solid_water_contact(solid, {&water_1});
//     ContactRelationFromShell shell_water_contact(shell, {&water_2});
//     ComplexRelation water_1_complex(water_inner_1, {&water_solid_contact});
//     ComplexRelation water_2_complex(water_inner_2, {&water_shell_contact});
//     //----------------------------------------------------------------------
//     //	Algorithm
//     //----------------------------------------------------------------------
//     const auto solid_velocity = [&]()
//     {
//         auto &solid_avg_vel = *solid.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
//         auto &solid_avg_force = *solid.getBaseParticles().getVariableByName<Vecd>("AverageForce");
//         auto &shell_avg_vel = *shell.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
//         auto &shell_avg_force = *shell.getBaseParticles().getVariableByName<Vecd>("AverageForce");
//         particle_for(
//             par,
//             solid.getBaseParticles().total_real_particles_,
//             [&](size_t index_i)
//             {
//                 solid_avg_vel[index_i] = Vecd(1, 1);
//                 solid_avg_force[index_i] = Vecd(100, 100) * solid.getBaseParticles().mass_[index_i];
//             });
//         particle_for(
//             par,
//             shell.getBaseParticles().total_real_particles_,
//             [&](size_t index_i)
//             {
//                 shell_avg_vel[index_i] = Vecd(1, 1);
//                 shell_avg_force[index_i] = Vecd(100, 100) * shell.getBaseParticles().mass_[index_i];
//             });
//     };
//     solid_velocity();

//     InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_by_summation_1(water_inner_1, water_solid_contact);
//     InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_by_summation_2(water_inner_2, water_shell_contact);
//     Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> water_pressure_relaxation_1(water_inner_1, water_solid_contact);
//     Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> water_pressure_relaxation_2(water_inner_2, water_shell_contact);
//     Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> water_density_relaxation_1(water_inner_1, water_solid_contact);
//     Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> water_density_relaxation_2(water_inner_2, water_shell_contact);
//     InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_1(water_inner_1, water_solid_contact);
//     InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_2(water_inner_2, water_shell_contact);
//     /** FSI */
//     InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(solid_water_contact);
//     InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_solid_update(solid_water_contact, viscous_force_on_solid);
//     InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_water_contact);
//     InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_shell_update(shell_water_contact, viscous_force_on_shell);

//     water_1.addBodyStateForRecording<Real>("Density");
//     water_1.addBodyStateForRecording<Real>("MassiveMeasure");
//     water_1.addBodyStateForRecording<Real>("VolumetricMeasure");
//     water_1.addBodyStateForRecording<Vecd>("Force");
//     water_1.addBodyStateForRecording<Vecd>("PriorForce");
//     water_2.addBodyStateForRecording<Real>("Density");
//     water_2.addBodyStateForRecording<Real>("MassiveMeasure");
//     water_2.addBodyStateForRecording<Real>("VolumetricMeasure");
//     water_2.addBodyStateForRecording<Vecd>("Force");
//     water_2.addBodyStateForRecording<Vecd>("PriorForce");
//     solid.addBodyStateForRecording<Real>("Density");
//     solid.addBodyStateForRecording<Real>("MassiveMeasure");
//     solid.addBodyStateForRecording<Real>("VolumetricMeasure");
//     solid.addBodyStateForRecording<Vecd>("NormalDirection");
//     solid.addBodyStateForRecording<Vecd>("AllForceFromFluid");
//     shell.addBodyStateForRecording<Real>("Density");
//     shell.addBodyStateForRecording<Real>("MassiveMeasure");
//     shell.addBodyStateForRecording<Real>("VolumetricMeasure");
//     shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
//     BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
//     //----------------------------------------------------------------------
//     //	Prepare the simulation with cell linked list, configuration
//     //	and case specified initial condition if necessary.
//     //----------------------------------------------------------------------
//     sph_system.initializeSystemCellLinkedLists();
//     sph_system.initializeSystemConfigurations();
//     water_shell_contact.updateConfiguration();
//     shell_water_contact.updateConfiguration();

//     //  Check dWijVjeij
//     CheckKernelCompleteness check_kernel_completeness_1(water_inner_1, water_solid_contact);
//     CheckKernelCompleteness check_kernel_completeness_2(water_inner_2, water_shell_contact);
//     check_kernel_completeness_1.exec();
//     check_kernel_completeness_2.exec();
//     water_1.addBodyStateForRecording<Real>("TotalKernel");
//     water_1.addBodyStateForRecording<Vecd>("TotalKernelGrad");
//     water_1.addBodyStateForRecording<int>("InnerNeighborNumber");
//     water_1.addBodyStateForRecording<int>("ContactNeighborNumber");
//     water_2.addBodyStateForRecording<Real>("TotalKernel");
//     water_2.addBodyStateForRecording<Vecd>("TotalKernelGrad");
//     water_2.addBodyStateForRecording<int>("InnerNeighborNumber");
//     water_2.addBodyStateForRecording<int>("ContactNeighborNumber");

//     std::cout << "Solid: Id: 0, contact number = " << check_kernel_completeness_1.get_contact_number(0)
//               << ", Shell: " << check_kernel_completeness_2.get_contact_number(0) << std::endl;
//     std::cout << "Solid: Id: 0, W_ijV_j = " << check_kernel_completeness_1.get_W_ijV_j(0)
//               << ", Shell: " << check_kernel_completeness_2.get_W_ijV_j(0) << std::endl;
//     std::cout << "Solid: Id: 0, dW_ijV_je_ij = "
//               << check_kernel_completeness_1.get_dW_ijV_je_ij(0).x()
//               << "\t" << check_kernel_completeness_1.get_dW_ijV_je_ij(0).y()
//               << "\t,Shell: "
//               << check_kernel_completeness_2.get_dW_ijV_je_ij(0).x()
//               << "\t" << check_kernel_completeness_2.get_dW_ijV_je_ij(0).y()
//               << std::endl;

//     const Real dt = 0.1;

//     update_fluid_density_by_summation_1.exec();
//     update_fluid_density_by_summation_2.exec();

//     std::cout << "After density summation: " << std::endl;
//     std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
//               << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;

//     viscous_acceleration_1.exec(dt);
//     viscous_acceleration_2.exec(dt);
//     std::cout << "After viscous_acceleration: " << std::endl;
//     std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
//               << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
//     std::cout << "Id: 0: solid, vel = "
//               << water_1.getBaseParticles().vel_[0].x()
//               << "\t" << water_1.getBaseParticles().vel_[0].y()
//               << "\t,Shell: "
//               << water_2.getBaseParticles().vel_[0].x()
//               << "\t" << water_2.getBaseParticles().vel_[0].y()
//               << std::endl;

//     water_pressure_relaxation_1.exec(dt);
//     water_pressure_relaxation_2.exec(dt);

//     std::cout << "After pressure relaxation: " << std::endl;
//     std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
//               << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
//     std::cout << "Id: 0: solid, vel = "
//               << water_1.getBaseParticles().vel_[0].x()
//               << "\t" << water_1.getBaseParticles().vel_[0].y()
//               << "\t,Shell: "
//               << water_2.getBaseParticles().vel_[0].x()
//               << "\t" << water_2.getBaseParticles().vel_[0].y()
//               << std::endl;

//     water_density_relaxation_1.exec(dt);
//     water_density_relaxation_2.exec(dt);

//     std::cout << "After density relaxation: " << std::endl;
//     std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
//               << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
//     std::cout << "Id: 0: solid, vel = "
//               << water_1.getBaseParticles().vel_[0].x()
//               << "\t" << water_1.getBaseParticles().vel_[0].y()
//               << "\t,Shell: "
//               << water_2.getBaseParticles().vel_[0].x()
//               << "\t" << water_2.getBaseParticles().vel_[0].y()
//               << std::endl;

//     viscous_force_on_solid.exec();
//     viscous_force_on_shell.exec();

//     const auto &viscous_force_solid = *solid.getBaseParticles().getVariableByName<Vecd>("ViscousForceFromFluid");
//     const auto &viscous_force_shell = *shell.getBaseParticles().getVariableByName<Vecd>("ViscousForceFromFluid");
//     Vecd total_viscous_force_solid(0, 0);
//     Vecd total_viscous_force_shell(0, 0);
//     for (size_t i = 0; i < viscous_force_solid.size(); i++)
//     {
//         total_viscous_force_solid += viscous_force_solid[i];
//     }
//     for (size_t i = 0; i < viscous_force_shell.size(); i++)
//     {
//         total_viscous_force_shell += viscous_force_shell[i];
//     }
//     std::cout << "Total force of solid: " << total_viscous_force_solid.x() << "\t" << total_viscous_force_solid.y() << std::endl;
//     std::cout << "Total force of shell: " << total_viscous_force_shell.x() * shell_thickness << "\t" << total_viscous_force_shell.y() * shell_thickness << std::endl;
//     // for (size_t i = 0; i < solid_positions.size(); i++)
//     //     std::cout << "Acceleration of solid: " << viscous_force_solid[i].x() / solid.getBaseParticles().mass_[i] << "\t" << viscous_force_solid[i].y() / solid.getBaseParticles().mass_[i] << std::endl;
//     // std::cout << "Acceleration of shell: " << viscous_force_shell[0].x() / shell.getBaseParticles().mass_[0] << "\t" << viscous_force_shell[0].y() / shell.getBaseParticles().mass_[0] << std::endl;

//     fluid_force_on_solid_update.exec();
//     fluid_force_on_shell_update.exec();

//     const auto &force_solid = *solid.getBaseParticles().getVariableByName<Vecd>("AllForceFromFluid");
//     const auto &force_shell = *shell.getBaseParticles().getVariableByName<Vecd>("AllForceFromFluid");
//     Vecd total_force_solid(0, 0);
//     Vecd total_force_shell(0, 0);
//     for (size_t i = 0; i < force_solid.size(); i++)
//     {
//         total_force_solid += force_solid[i];
//     }
//     for (size_t i = 0; i < force_shell.size(); i++)
//     {
//         total_force_shell += force_shell[i];
//     }
//     std::cout << "Total force of solid: " << total_force_solid.x() << "\t" << total_force_solid.y() << std::endl;
//     std::cout << "Total force of shell: " << total_force_shell.x() * shell_thickness << "\t" << total_force_shell.y() * shell_thickness << std::endl;
//     // for (size_t i = 0; i < solid_positions.size(); i++)
//     //     std::cout << "Acceleration of solid: " << force_solid[i].x() / solid.getBaseParticles().mass_[i] << "\t" << force_solid[i].y() / solid.getBaseParticles().mass_[i] << std::endl;
//     // std::cout << "Acceleration of shell: " << force_shell[0].x() / shell.getBaseParticles().mass_[0] << "\t" << force_shell[0].y() / shell.getBaseParticles().mass_[0] << std::endl;

//     // water_1.setNewlyUpdated();
//     // water_2.setNewlyUpdated();

//     write_body_states.writeToFile();
//     return 0;
// }

int main(int ac, char *av[])
{
    const Real scale = 0.001;
    const Real resolution_ref = 0.5 * scale;
    const Real resolution_shell = 0.5 * resolution_ref;
    const Real shell_thickness = 6 * resolution_ref;
    const Real radius_outer = 100 * resolution_shell / (2 * Pi);
    const Real radius_mid = radius_outer - 0.5 * resolution_shell;
    const Vecd center(0, -radius_outer);
    // const Real half_length = 6 * resolution_ref;
    //----------------------------------------------------------------------
    //	Global parameters on the fluid properties
    //----------------------------------------------------------------------
    const Real rho0_f = 1056.0;                           /**< Reference density of fluid. */
    const Real mu_f = 3.5e-3;                             /**< Dynamics viscosity. */
    const Real Re = 100.0;                                /**< Reynolds number. */
    const Real U_f = Re * mu_f / rho0_f / (12.0 * scale); /**< Characteristic velocity. */
    const Real U_max = 1.5 * U_f;
    /** Reference sound speed needs to consider the flow speed in the narrow channels. */
    const Real c_f = 10.0 * U_max;
    //----------------------------------------------------------------------
    //	Global parameters on the solid properties
    //----------------------------------------------------------------------
    const Real rho0_s = 1250.0; /**< Reference density.*/
    const Real hardness = 50;   // Durometer hardnes: 50A
    const Real youngs_modulus =
        std::pow(10, 0.0235 * hardness - 0.6403) * 1e4; // actual: 1e6, ref: https://doi.org/10.5254/1.3547752 eq. 12A
    const Real poisson_ratio = 0.495;

    // BoundingBox system_domain_bounds(Vec2d(-half_length, -6 * resolution_ref), Vec2d(half_length, resolution_ref));
    BoundingBox system_domain_bounds(Vec2d(-radius_outer, -2 * radius_outer), Vec2d(radius_outer, 2 * resolution_ref));
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    StdLargeVec<Vecd> water_positions({Vec2d(0, 0.5 * resolution_ref),
                                       Vec2d(0, 1.5 * resolution_ref)});
    StdLargeVec<Real> water_volumes({resolution_ref * resolution_ref, resolution_ref * resolution_ref});
    FluidBody water_1(sph_system, makeShared<DefaultShape>("WaterBody1"));
    water_1.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_1.generateParticles<ParticleGeneratorDirect>(water_positions, water_volumes);

    FluidBody water_2(sph_system, makeShared<DefaultShape>("WaterBody2"));
    water_2.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_2.generateParticles<ParticleGeneratorDirect>(water_positions, water_volumes);

    SolidBody solid(sph_system, makeShared<CircleSolidShape>("Solid", radius_outer, shell_thickness, center));
    // SolidBody solid(sph_system, makeShared<PlateSolidShape>("Solid", half_length, shell_thickness));
    solid.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    // solid.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    solid.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    // solid.generateParticles<ParticleGeneratorLattice>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? solid.generateParticles<ParticleGeneratorReload>(io_environment, solid.getName())
        : solid.generateParticles<ParticleGeneratorLattice>();

    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    shell.generateParticles<CircleShellParticleGeneratorDirect>(radius_mid, resolution_shell, shell_thickness, center);
    // shell.generateParticles<PlateShellParticleGeneratorDirect>(half_length, resolution_shell, shell_thickness);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    //   if (sph_system.RunParticleRelaxation())
    //   {
    //       //----------------------------------------------------------------------
    //       //	Define body relation map used for particle relaxation.
    //       //----------------------------------------------------------------------
    //       InnerRelation solid_inner(solid);
    //       //----------------------------------------------------------------------
    //       //	Define the methods for particle relaxation for wall boundary.
    //       //----------------------------------------------------------------------
    //       SimpleDynamics<RandomizeParticlePosition> solid_random_particles(solid);
    //       relax_dynamics::RelaxationStepInner relaxation_step_solid_inner(solid_inner);
    //       //----------------------------------------------------------------------
    //       //	Output for particle relaxation.
    //       //----------------------------------------------------------------------
    //       // BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
    //       ReloadParticleIO write_particle_reload_files(io_environment, {&solid});
    //       //----------------------------------------------------------------------
    //       //	Particle relaxation starts here.
    //       //----------------------------------------------------------------------
    //       solid_random_particles.exec(0.25);
    //       relaxation_step_solid_inner.SurfaceBounding().exec();
    //       // write_relaxed_particles.writeToFile(0);
    //       //----------------------------------------------------------------------
    //       //	From here iteration for particle relaxation begins.
    //       //----------------------------------------------------------------------
    //       int ite = 0;
    //       int relax_step = 2000;
    //       while (ite < relax_step)
    //       {
    //           relaxation_step_solid_inner.exec();
    //           ite += 1;
    //           if (ite % 100 == 0)
    //           {
    //               std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
    //               // write_relaxed_particles.writeToFile(ite);
    //           }
    //       }
    //       std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
    //       write_particle_reload_files.writeToFile(0);
    //       return 0;
    //   }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //    Generally, we first define all the inner relations, then the contact relations.
    //    At last, we define the complex relaxations by combining previous defined
    //    inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_inner_1(water_1);
    InnerRelation water_inner_2(water_2);
    InnerRelation solid_inner(solid);
    InnerRelation shell_inner(shell);
    ContactRelation water_solid_contact(water_1, {&solid});
    ContactRelationToShell water_shell_contact(water_2, {&shell});
    ContactRelation solid_water_contact(solid, {&water_1});
    ContactRelationFromShell shell_water_contact(shell, {&water_2});
    ComplexRelation water_1_complex(water_inner_1, {&water_solid_contact});
    ComplexRelation water_2_complex(water_inner_2, {&water_shell_contact});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, water_2);
    //----------------------------------------------------------------------
    //	Algorithm
    //----------------------------------------------------------------------
    const auto solid_velocity = [&]()
    {
        auto &solid_avg_vel = *solid.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
        auto &solid_avg_force = *solid.getBaseParticles().getVariableByName<Vecd>("AverageForce");
        auto &shell_avg_vel = *shell.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
        auto &shell_avg_force = *shell.getBaseParticles().getVariableByName<Vecd>("AverageForce");
        particle_for(
            par,
            solid.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                solid_avg_vel[index_i] = Vecd(1, 1);
                solid_avg_force[index_i] = Vecd(100, 100) * solid.getBaseParticles().mass_[index_i];
            });
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                shell_avg_vel[index_i] = Vecd(1, 1);
                shell_avg_force[index_i] = Vecd(100, 100) * shell.getBaseParticles().mass_[index_i];
            });
    };
    solid_velocity();

    SimpleDynamics<NormalDirectionFromBodyShape> solid_normal_direction(solid);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);

    SimpleDynamics<TimeStepInitialization> fluid_step_initialization_1(water_1);
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization_2(water_2);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_by_summation_1(water_inner_1, water_solid_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_by_summation_2(water_inner_2, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> water_pressure_relaxation_1(water_inner_1, water_solid_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> water_pressure_relaxation_2(water_inner_2, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> water_density_relaxation_1(water_inner_1, water_solid_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> water_density_relaxation_2(water_inner_2, water_shell_contact);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_1(water_inner_1, water_solid_contact);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_2(water_inner_2, water_shell_contact);
    /** FSI */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(solid_water_contact);
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_force_on_solid_update(solid_water_contact);
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_water_contact);
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_force_on_shell_update(shell_water_contact);

    water_1.addBodyStateForRecording<Real>("Density");
    water_1.addBodyStateForRecording<Real>("MassiveMeasure");
    water_1.addBodyStateForRecording<Real>("VolumetricMeasure");
    water_1.addBodyStateForRecording<Vecd>("Force");
    water_1.addBodyStateForRecording<Vecd>("PriorForce");
    water_2.addBodyStateForRecording<Real>("Density");
    water_2.addBodyStateForRecording<Real>("MassiveMeasure");
    water_2.addBodyStateForRecording<Real>("VolumetricMeasure");
    water_2.addBodyStateForRecording<Vecd>("Force");
    water_2.addBodyStateForRecording<Vecd>("PriorForce");
    solid.addBodyStateForRecording<Real>("Density");
    solid.addBodyStateForRecording<Real>("MassiveMeasure");
    solid.addBodyStateForRecording<Real>("VolumetricMeasure");
    solid.addBodyStateForRecording<Vecd>("NormalDirection");
    solid.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    solid.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    shell.addBodyStateForRecording<Real>("Density");
    shell.addBodyStateForRecording<Real>("MassiveMeasure");
    shell.addBodyStateForRecording<Real>("VolumetricMeasure");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    shell.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    solid_normal_direction.exec();
    shell_average_curvature.exec();
    water_shell_contact.updateConfiguration();

    // const Real curvature = 1 / radius_mid;
    // const auto &curvature_ave = *shell.getBaseParticles().getVariableByName<Real>("Average1stPrincipleCurvature");
    /// std::cout << "Analytical curvature: " << curvature << ", average curvature: " << curvature_ave[0] << std::endl;

    //  Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness_1(water_inner_1, water_solid_contact);
    CheckKernelCompleteness check_kernel_completeness_2(water_inner_2, water_shell_contact);
    check_kernel_completeness_1.exec();
    check_kernel_completeness_2.exec();
    water_1.addBodyStateForRecording<Real>("TotalKernel");
    water_1.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_1.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_1.addBodyStateForRecording<int>("ContactNeighborNumber");
    water_2.addBodyStateForRecording<Real>("TotalKernel");
    water_2.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_2.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_2.addBodyStateForRecording<int>("ContactNeighborNumber");

    std::cout << "Solid: Id: 0, contact number = " << check_kernel_completeness_1.get_contact_number(0)
              << ", Shell: " << check_kernel_completeness_2.get_contact_number(0) << std::endl;
    std::cout << "Solid: Id: 1, contact number = " << check_kernel_completeness_1.get_contact_number(1)
              << ", Shell: " << check_kernel_completeness_2.get_contact_number(1) << std::endl;
    std::cout << "Solid: Id: 0, W_ijV_j = " << check_kernel_completeness_1.get_W_ijV_j(0)
              << ", Shell: " << check_kernel_completeness_2.get_W_ijV_j(0) << std::endl;
    std::cout << "Solid: Id: 1, W_ijV_j = " << check_kernel_completeness_1.get_W_ijV_j(1)
              << ", Shell: " << check_kernel_completeness_2.get_W_ijV_j(1) << std::endl;
    std::cout << "Solid: Id: 0, dW_ijV_je_ij = "
              << check_kernel_completeness_1.get_dW_ijV_je_ij(0).x()
              << "\t" << check_kernel_completeness_1.get_dW_ijV_je_ij(0).y()
              << "\t,Shell: "
              << check_kernel_completeness_2.get_dW_ijV_je_ij(0).x()
              << "\t" << check_kernel_completeness_2.get_dW_ijV_je_ij(0).y()
              << std::endl;
    std::cout << "Solid: Id: 1, dW_ijV_je_ij = "
              << check_kernel_completeness_1.get_dW_ijV_je_ij(1).x()
              << "\t" << check_kernel_completeness_1.get_dW_ijV_je_ij(1).y()
              << "\t,Shell: "
              << check_kernel_completeness_2.get_dW_ijV_je_ij(1).x()
              << "\t" << check_kernel_completeness_2.get_dW_ijV_je_ij(1).y()
              << std::endl;

    const Real dt = 0.1;

    update_fluid_density_by_summation_1.exec();
    update_fluid_density_by_summation_2.exec();

    std::cout << "After density summation: " << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
              << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[1]
              << ", Shell: " << water_2.getBaseParticles().rho_[1] << std::endl;

    viscous_acceleration_1.exec(dt);
    viscous_acceleration_2.exec(dt);
    std::cout << "After viscous_acceleration: " << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
              << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[1]
              << ", Shell: " << water_2.getBaseParticles().rho_[1] << std::endl;
    std::cout << "Id: 0: solid, vel = "
              << water_1.getBaseParticles().vel_[0].x()
              << "\t" << water_1.getBaseParticles().vel_[0].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[0].x()
              << "\t" << water_2.getBaseParticles().vel_[0].y()
              << std::endl;
    std::cout << "Id: 1: solid, vel = "
              << water_1.getBaseParticles().vel_[1].x()
              << "\t" << water_1.getBaseParticles().vel_[1].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[1].x()
              << "\t" << water_2.getBaseParticles().vel_[1].y()
              << std::endl;

    water_pressure_relaxation_1.exec(dt);
    water_pressure_relaxation_2.exec(dt);

    std::cout << "After pressure relaxation: " << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
              << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[1]
              << ", Shell: " << water_2.getBaseParticles().rho_[1] << std::endl;
    std::cout << "Id: 0: solid, vel = "
              << water_1.getBaseParticles().vel_[0].x()
              << "\t" << water_1.getBaseParticles().vel_[0].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[0].x()
              << "\t" << water_2.getBaseParticles().vel_[0].y()
              << std::endl;
    std::cout << "Id: 1: solid, vel = "
              << water_1.getBaseParticles().vel_[1].x()
              << "\t" << water_1.getBaseParticles().vel_[1].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[1].x()
              << "\t" << water_2.getBaseParticles().vel_[1].y()
              << std::endl;

    water_density_relaxation_1.exec(dt);
    water_density_relaxation_2.exec(dt);

    std::cout << "After density relaxation: " << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[0]
              << ", Shell: " << water_2.getBaseParticles().rho_[0] << std::endl;
    std::cout << "Id: 0: solid rho = " << water_1.getBaseParticles().rho_[1]
              << ", Shell: " << water_2.getBaseParticles().rho_[1] << std::endl;
    std::cout << "Id: 0: solid, vel = "
              << water_1.getBaseParticles().vel_[0].x()
              << "\t" << water_1.getBaseParticles().vel_[0].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[0].x()
              << "\t" << water_2.getBaseParticles().vel_[0].y()
              << std::endl;
    std::cout << "Id: 1: solid, vel = "
              << water_1.getBaseParticles().vel_[1].x()
              << "\t" << water_1.getBaseParticles().vel_[1].y()
              << "\t,Shell: "
              << water_2.getBaseParticles().vel_[1].x()
              << "\t" << water_2.getBaseParticles().vel_[1].y()
              << std::endl;

    viscous_force_on_solid.exec();
    viscous_force_on_shell.exec();
    auto get_total_viscous_force = [&](SolidBody &body)
    {
        Vecd total_force(0, 0);
        const auto &force_from_fluid = *body.getBaseParticles().getVariableByName<Vec2d>("ViscousForceFromFluid");
        for (size_t i = 0; i < body.getBaseParticles().total_real_particles_; i++)
            total_force += force_from_fluid[i];
        return total_force;
    };

    fluid_force_on_solid_update.exec();
    fluid_force_on_shell_update.exec();
    auto get_total_pressure_force = [&](SolidBody &body)
    {
        Vecd total_force(0, 0);
        const auto &force_from_fluid = *body.getBaseParticles().getVariableByName<Vec2d>("PressureForceFromFluid");
        for (size_t i = 0; i < body.getBaseParticles().total_real_particles_; i++)
            total_force += force_from_fluid[i];
        return total_force;
    };

    auto solid_viscous_force = get_total_viscous_force(solid);
    auto shell_viscous_force = get_total_viscous_force(shell);
    auto solid_pressure_force = get_total_pressure_force(solid);
    auto shell_pressure_force = get_total_pressure_force(shell);
    std::cout << "Total viscous force of solid: " << solid_viscous_force.x() << "\t" << solid_viscous_force.y() << std::endl;
    std::cout << "Total viscous force of shell: " << shell_viscous_force.x() * shell_thickness << "\t" << shell_viscous_force.y() * shell_thickness << std::endl;
    std::cout << "Total pressure force of solid: " << solid_pressure_force.x() << "\t" << solid_pressure_force.y() << std::endl;
    std::cout << "Total pressure force of shell: " << shell_pressure_force.x() * shell_thickness << "\t" << shell_pressure_force.y() * shell_thickness << std::endl;

    water_1.setNewlyUpdated();
    water_2.setNewlyUpdated();
    shell.setNewlyUpdated();
    solid.setNewlyUpdated();

    write_body_states.writeToFile();
    return 0;
}
