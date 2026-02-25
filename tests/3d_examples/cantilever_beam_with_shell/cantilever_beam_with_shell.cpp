// Copyright 2022 Virtonomy GmbH
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include <gtest/gtest.h>
#include "sphinxsys.h"

using namespace SPH;

constexpr Real eps = 1e-7;

const Real scaling = 1000;
// Geometry
const Real length = 0.2*scaling;
const Real width = 0.02*scaling;
const Real thickness = 0.002*scaling;
const Vecd n0 = {0,0,1};
// Material
Real rho = {};
Real Youngs_modulus = {};
Real poisson = {};
Real physical_viscosity = {};
// BC
Real ratio_fixed = {};// Percent of body fixed
bool extend_bc = {};// Extend fixed bc by kernel cutoff
Real linear_load = {};// 10 N/m
Real pressure = {}; // Pressure
// Misc
bool particles_on_boundaries = {};; // cell vs node location of particles
bool use_exact_volume = {}; // for beam generated with particles on boundaries

void reset_parameters()
{
    rho = 7800/(scaling*scaling*scaling);
    Youngs_modulus = 207e9/(scaling*scaling);
    poisson = 0.49;
    physical_viscosity = 1./ 4.0 * std::sqrt(rho * Youngs_modulus)*0.002*scaling;
    ratio_fixed = 0.1;// Percent of body fixed
    extend_bc = false;// Extend fixed bc by kernell cutoff
    linear_load = 100./scaling;// 10 N/m
    pressure = linear_load/width; // Pressure
    particles_on_boundaries = false; // cell vs node location of particles
    use_exact_volume = false; // for beam generated with particles on boundaries
}

Real compute_analytical_deflection()
{
    // https://youtu.be/eo3aUJKvnTw?t=1765 Abaqus error at 1.17% for E=207GPa and nu=0.3
    // Closer to beam theory -> rigid cross-section -> Poisson's ratio should be 0
    auto linear_load = pressure*width;
    auto moment_inertia = width*std::pow(thickness,3)/12; //rectangular cross-section
    auto deflection = linear_load*std::pow(length,4)/(8*Youngs_modulus*moment_inertia);
    return deflection;
}

//=============================================================================================//

struct ApplyLoad : public PartSimpleDynamicsByParticle
{
    ApplyLoad(SPHBody &sph_body, BodyPartByParticle &body_part, Real p, Vec3d dir) : PartSimpleDynamicsByParticle(sph_body,body_part), p(p), dir(dir) {}

    virtual void Update(size_t index_i, Real dt = 0.0) override
    {
        // auto thickness = static_cast<ShellParticles&>(*base_particles_).shell_thickness_[index_i];
        // auto normal = static_cast<ShellParticles&>(*base_particles_).n_[index_i];
        auto acceleration = -(p*base_particles_->Vol_[index_i])/(base_particles_->mass_[index_i]*thickness);
        base_particles_->acc_prior_[index_i] += -acceleration*dir;
    }
    Real p;
    Vec3d dir;
};

class BeamParticleGenerator : public SurfaceParticleGenerator
{
public:
	BeamParticleGenerator(SPHBody &sph_body, Real max_interparticle_distance, Real length, Real width, Vec3d origin) : 
        SurfaceParticleGenerator(sph_body),
        max_interparticle_distance(max_interparticle_distance),
        beam_length(length),
        beam_width(width),
        origin(origin)
	{
	}

    Real max_interparticle_distance; 
    Real beam_length; 
    Real beam_width;
    Vec3d origin;

    virtual void initializeGeometricVariables() override
    {
        auto dx = (beam_length-eps) / max_interparticle_distance;
        auto nx = int(std::ceil(dx));
        dx = beam_length/nx;
        auto dy = (beam_width-eps) / max_interparticle_distance;
        auto ny = int(std::ceil(dy));
        dy = beam_width/ny;

        if(particles_on_boundaries)
        {
            for (int i = 0; i <= nx ; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    Real x = origin.get(0) + i*dx;
                    Real y = origin.get(1) + j*dy;
                    Real z = origin.get(2);
                    Real v = dx*dy;
                    if(use_exact_volume)
                    {
                        if(i==0 || j==0 || i==nx || j ==ny) v = dx*dy/2;
                        if(i==0 && j==0) v = dx*dy/4;
                        if(i==nx && j==0) v = dx*dy/4;
                        if(i==0 && j==ny) v = dx*dy/4;
                        if(i==nx && j==ny) v = dx*dy/4;
                    }
                    initializePositionAndVolumetricMeasure(Vecd(x, y, z), v);
                    initializeSurfaceProperties(n0, thickness);
                }
            }
        }
        else
        {
            for (int i = 0; i < nx ; i++)
            {// dp/2 from border as recommended by upstream dev
                for (int j = 0; j < ny; j++)
                {
                    Real x = origin.get(0) + i*dx + dx/2;
                    Real y = origin.get(1) + j*dy + dy/2;
                    Real z = origin.get(2);
                    Real v = dx*dy;
                    initializePositionAndVolumetricMeasure(Vecd(x, y, z), v);
                    initializeSurfaceProperties(n0, thickness);
                }
            }
        }
    }
};

//=============================================================================================//

Real cantilever_beam_with_shell(Real dp_shell)
{
	// Units SI: m, Pa, N, kg
	Real end_time = 5.0;

    SPHSystem system({Vec3d{-ratio_fixed*length,0.0,-thickness}-eps,Vec3d{length,width,thickness}+eps}, dp_shell);
	InOutput in_output(system);
    
	// Shell body
	SolidBody shell_body (system, makeShared<DefaultShape>("Body"));
	shell_body.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho, Youngs_modulus, poisson);
	shell_body.generateParticles<BeamParticleGenerator>(dp_shell,(1+ratio_fixed)*length,width,Vec3d{-ratio_fixed*length,0,0});
	BodyRelationInner shell_body_inner(shell_body);

	// Shell algorithms
	TimeStepInitialization initialize_shell(shell_body);
	thin_structure_dynamics::ShellCorrectConfiguration corrected_configuration(shell_body_inner);
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(shell_body);
	thin_structure_dynamics::ShellStressRelaxationFirstHalf stress_relaxation_first_half(shell_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf stress_relaxation_second_half(shell_body_inner);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>> shell_damping(0.2,shell_body_inner,"Velocity", physical_viscosity);

    // Constrain the base
    auto shell_constraint_fixed_shape = [&]()
    {
        Real cutoff_length = (1-eps)*2.0*shell_body.sph_adaptation_->ReferenceSpacing()*shell_body.sph_adaptation_->ReferenceSmoothingLength();// number of extra particle rows fixed beyond the fixed block, to guarantee displacement(0,y,z) = 0 because of basis smoothness
        if(!extend_bc) 
            cutoff_length = 0;
        auto origin = Vec3d(-ratio_fixed*length, width/2., 0.0);
        auto half_size = Vec3d((ratio_fixed*length)+cutoff_length, width/2.+eps, thickness);
        return makeShared<TriangleMeshShapeBrick>(half_size, 20, origin);
    }();
    BodyRegionByParticle shell_constraint_fixed_part(shell_body, shell_constraint_fixed_shape);
    thin_structure_dynamics::ConstrainShellBodyRegion shell_constraint_fixed(shell_body, shell_constraint_fixed_part);
    
    // Apply pressure on the body of the beam
    auto shell_constraint_force_shape = [&]()
    {
        auto origin = Vec3d(length/2., width, -eps);
        auto half_size = Vec3d(length/2, width, thickness);
        return makeShared<TriangleMeshShapeBrick>(half_size, 20, origin);
    }();
    BodyRegionByParticle shell_constraint_force_part(shell_body, shell_constraint_force_shape);
    ApplyLoad shell_constraint_force(shell_body, shell_constraint_force_part, pressure, Vecd(0,0,-1));

	// Output
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);

	// Initialization
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	// Update the normal directions, initial conditions, and configurations
	corrected_configuration.parallel_exec();
	write_states.writeToFile(0);

	// Step parameters
	int ite = 0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;

	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	std::vector<Real> iteration, deflection;
    Real current_max_deflection = 0;
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite 
                          << " | Time: " << GlobalStaticVariables::physical_time_ 
                          << " | dt: " << dt 
                          << " | max displ:" << current_max_deflection <<"                                 \r";
			}

			initialize_shell.parallel_exec(dt);

            shell_constraint_force.parallel_exec(dt);

			stress_relaxation_first_half.parallel_exec(dt);

            shell_constraint_fixed.parallel_exec();

            shell_damping.parallel_exec(dt);

            shell_constraint_fixed.parallel_exec();
            
            stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
            
            current_max_deflection = static_cast<ShellParticles*>(shell_body.base_particles_)->getMaxDisplacement();
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;

		iteration.push_back(iteration.size());
		deflection.push_back(current_max_deflection);
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "\nTotal wall clock time: " << int(tt.seconds()) << " seconds | Deflection: "<<deflection.back()<<" mm\n";
	return deflection.back();
}

//=============================================================================================//
// Tests
//=============================================================================================//

std::vector<Real> resolutions = {thickness*2,thickness*1,thickness/2};

TEST(cantilever_beam_with_shell, cell_location)
{
    reset_parameters();

	std::vector<Real> deflections;
    auto reference = compute_analytical_deflection();
	for (auto dp: resolutions)
    {
        std::cout<<"Resolution: "<<dp<<" mm\n";
        Real deflection = cantilever_beam_with_shell(dp);
        deflections.emplace_back(deflection);
        auto relative_error = std::abs(deflection-reference)/reference*100;
        EXPECT_LE(relative_error, 5.0/*%*/);
    }
}

TEST(cantilever_beam_with_shell, node_location)
{
    reset_parameters();
    particles_on_boundaries = true;
    
    std::vector<Real> deflections;
    auto reference = compute_analytical_deflection();
	for (auto dp: resolutions)
    {
        std::cout<<"Resolution: "<<dp<<" mm\n";
        Real deflection = cantilever_beam_with_shell(dp);
        deflections.emplace_back(deflection);
        auto relative_error = std::abs(deflection-reference)/reference*100;
        EXPECT_LE(relative_error, 5.0/*%*/);
    }
}

TEST(cantilever_beam_with_shell, cell_location_with_extended_bc)
{
    reset_parameters();
    extend_bc = true;

    std::vector<Real> deflections;
    auto reference = compute_analytical_deflection();
	for (auto dp: resolutions)
    {
        std::cout<<"Resolution: "<<dp<<" mm\n";
        Real deflection = cantilever_beam_with_shell(dp);
        deflections.emplace_back(deflection);
        auto relative_error = std::abs(deflection-reference)/reference*100;
        EXPECT_LE(relative_error, 5.0/*%*/);
    }
}

TEST(cantilever_beam_with_shell, node_location_with_extended_bc)
{
    reset_parameters();
    particles_on_boundaries = true;
    extend_bc = true;
    
    std::vector<Real> deflections;
    auto reference = compute_analytical_deflection();
	for (auto dp: resolutions)
    {
        std::cout<<"Resolution: "<<dp<<" mm\n";
        Real deflection = cantilever_beam_with_shell(dp);
        deflections.emplace_back(deflection);
        auto relative_error = std::abs(deflection-reference)/reference*100;
        EXPECT_LE(relative_error, 5.0/*%*/);
    }
}

TEST(cantilever_beam_with_shell, node_location_with_extended_bc_and_exact_volume)
{
    reset_parameters();
    particles_on_boundaries = true;
    extend_bc = true;
    use_exact_volume = true;

    std::vector<Real> deflections;
    auto reference = compute_analytical_deflection();
	for (auto dp: resolutions)
    {
        std::cout<<"Resolution: "<<dp<<" mm\n";
        Real deflection = cantilever_beam_with_shell(dp);
        deflections.emplace_back(deflection);
        auto relative_error = std::abs(deflection-reference)/reference*100;
        EXPECT_LE(relative_error, 5.0/*%*/);
    }
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
